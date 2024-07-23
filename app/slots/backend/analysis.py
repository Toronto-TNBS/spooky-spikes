import scipy.signal as sps
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.mixture import GaussianMixture
from scipy.stats import median_abs_deviation
import warnings
import math
from astropy.timeseries import LombScargle
import finnpy.filters.frequency as ff
import finnpy.basic.downsampling as ds


def butterworth_filter(signal, fs, hpcutoff=None, lpcutoff=None):
    if hpcutoff == None:
        b, a = sps.iirfilter(4, lpcutoff, btype='lowpass', fs=fs)
    elif lpcutoff == None:
        b, a = sps.iirfilter(4, hpcutoff, btype='highpass', fs=fs)
    else:
        b, a = sps.iirfilter(4, [hpcutoff, lpcutoff], fs=fs)
    
    return sps.filtfilt(b, a, signal)


def get_spike_threshold(signal, factor):
    return factor * np.median(np.abs(signal)) / 0.6745


def get_spike_indices(signal, threshold, fs, inverted):
    if inverted:
        signal = -signal.copy()    # Signal inverted for spike detection, factor and threshold remain positive.
    spike_indices = sps.find_peaks(signal, height=threshold, distance=int(2/1000*fs))[0]
    window_size = int(1/1000*fs)    # Standardize this spike-matrix windowing rule so features remain consistent.
    spike_indices_adjusted = spike_indices[(spike_indices > window_size) & (spike_indices < spike_indices[-1]-window_size)]
    return spike_indices_adjusted


def get_spike_matrix(signal, spike_indices, window_size):
    spike_matrix = np.empty((1, window_size*2))
    for i in spike_indices:
        spike_matrix = np.append(spike_matrix, [signal[i-window_size: i+window_size]], axis=0)
    spike_matrix = np.delete(spike_matrix, 0, axis=0)
    return spike_matrix


def spike_sorting(spike_matrix, num_clusters):

    def clustering(spike_matrix, num_clusters):
        pca_results = PCA(n_components=2).fit_transform(spike_matrix)
        kmeans = KMeans(n_clusters=num_clusters, random_state=0).fit(pca_results)
        cluster_labels = kmeans.labels_
        return pca_results, cluster_labels

    if num_clusters is None:
        scores = []
        for c in range(2, 7):
            features_temp, label_temp = clustering(spike_matrix, c)
            scores.append(silhouette_score(features_temp, label_temp))
        num_clusters = np.array(scores).argmax() + 2

    features, labels = clustering(spike_matrix, num_clusters)
    score = silhouette_score(features, labels)

    return features, labels, score


def get_FR(spike_indices, fs):
    return fs / np.mean(np.diff(spike_indices))


def get_BI(spike_indices):
    isi = np.diff(spike_indices)
    X = np.ravel(isi).reshape(-1, 1)
    M_best = GaussianMixture(n_components=2, covariance_type="spherical", random_state=1).fit(X)
    means = M_best.means_
    return np.max(means) / np.min(means)


def get_CV(spike_indices):
    isi = np.diff(spike_indices)
    return median_abs_deviation(isi) / np.median(isi)


def get_SNR(signal, spike_indices, fs):
    b, a = sps.iirfilter(4, [300, 3000], fs=fs)
    filtered_data = sps.filtfilt(b, a, signal)
    waveforms = get_spike_matrix(filtered_data, spike_indices, window_size=int(1/1000*fs))
    mean_waveform = waveforms.mean(axis=0)
    std_waveform = waveforms.std(axis=0).mean(axis=-1)
    peak_range = mean_waveform.max(axis=-1) - mean_waveform.min(axis=-1)
    noise = 2 * std_waveform
    return peak_range / noise


def get_pISIv(spike_indices, fs):
    spiketrain = spike_indices / fs
    isi = np.diff(spiketrain)
    if (isi < 0).any():
        warnings.warn("ISI evaluated to negative values. "
                        "Please sort the input array.")
    return (sum(isi < 1 / 1000) / len(isi)) * 100


def get_spiketrain_power(spike_indices, fs):

    def autocorrelation(spike_times, bin_size, max_lag):
        bins = np.arange(0, np.max(spike_times) + bin_size, bin_size)
        spike_counts = np.histogram(spike_times, bins=bins)[0]
        autocorr = np.correlate(spike_counts, spike_counts, mode='full')[len(spike_counts) - 1:]
        # Following 2 lines are potential faster replacement for above line.
        # import scipy
        # autocorr = scipy.signal.correlate(spike_counts, spike_counts)
        autocorr_lag = np.arange(0, len(autocorr)) * bin_size
        autocorr = autocorr[autocorr_lag <= max_lag]
        autocorr_lag = autocorr_lag[autocorr_lag <= max_lag]
        return autocorr, autocorr_lag
        
    def kaneoke_oscillation_power(spike_times, low, high, bin_size, max_lag):
        autocorr, autocorr_lag = autocorrelation(spike_times, bin_size, max_lag)
        freqs = np.linspace(4, 30, num=1000)
        ls = LombScargle(autocorr_lag, autocorr)
        power = ls.power(freqs)
        # Find the maximum peak in the specified frequency band
        low_freq, high_freq = (low, high)
        idx = np.where((freqs >= low_freq) & (freqs < high_freq))[0]
        peaks, _ = sps.find_peaks(power[idx])
        if peaks.size > 0:
            max_peak_idx = idx[peaks[np.argmax(power[idx][peaks])]]
            max_peak_power = power[max_peak_idx]
            return 10 * math.log10(max_peak_power)
        else:
            return 10 * math.log10(max(power[idx]))
            
    spike_times = spike_indices / fs
    theta = kaneoke_oscillation_power(spike_times, 4, 8, 10e-3, 500e-3)
    alpha = kaneoke_oscillation_power(spike_times, 8, 12, 10e-3, 500e-3)
    lowbeta = kaneoke_oscillation_power(spike_times, 12, 21, 10e-3, 500e-3)
    highbeta = kaneoke_oscillation_power(spike_times, 21, 30, 10e-3, 500e-3)
    return theta, alpha, lowbeta, highbeta


def get_spiketrain_burstduration(spike_indices, fs):

    def spiketrain_to_binary(spiketimes, sampling_rate):
        num_timesteps = int(np.ceil(np.max(spiketimes) * sampling_rate))
        binary_sequence = np.zeros(num_timesteps)
        spike_indices = (spiketimes * sampling_rate).astype(int)
        valid_spike_indices = spike_indices[spike_indices < num_timesteps]
        binary_sequence[valid_spike_indices] = 1
        return binary_sequence

    def waveform_spiketrain_oscillation(spike_train, low_freq, high_freq):

        binary_sequence = spiketrain_to_binary(spike_train, 500)

        spiketrain_waveform = ff.butter(binary_sequence, low_freq, high_freq, 500, order=4, zero_phase=True)

        spiketrain_waveform = spiketrain_waveform - np.mean(spiketrain_waveform)
        freqs, psd = sps.welch(spiketrain_waveform, 500, nperseg=500 * 2, scaling='density')
        power = np.max(psd)
        time = np.arange(0, len(spiketrain_waveform)) / 500

        return spiketrain_waveform, time, 10 * math.log10(power)

    def get_waveform_envelope(waveform, smooth=False):
        burst = np.abs(sps.hilbert(waveform))
        if smooth > 0:
            window_size = smooth
            window = np.ones(window_size) / window_size
            burst = sps.convolve(burst, window, mode='same')
            return burst
        else:
            return burst
    
    def burst_threshold(data, fs_lfp=None, data_type='lfp'):
        if data_type == 'lfp':
            assert fs_lfp is not None, 'fs_lfp must be provided for lfp data'
            gamma_envelopes = [np.abs(sps.hilbert(ff.butter(data, f1, f2, fs_lfp, order=8, zero_phase=True)))
                               for f1, f2 in zip(range(45, 50), range(51, 56))]
        elif data_type == 'spiketrain':
            gamma_envelopes = [np.abs(sps.hilbert(waveform_spiketrain_oscillation(data, f1, f2)[0]))
                               for f1, f2 in zip(range(45, 50), range(51, 56))]
        else:
            raise ValueError('data_type must be either "lfp" or "spiketrain"')
        median_values = [np.median(data[sps.argrelextrema(data, np.less)]) for data in gamma_envelopes]
        return 4 * np.mean(median_values)
    
    def burst_features(time_ds, burst, burst_threshold):
        above_threshold = burst > burst_threshold
        cross_indices = np.where(np.diff(above_threshold))[0]
        if len(cross_indices) <= 1:
            return np.zeros(7), 0
        else:
            if not above_threshold[cross_indices[0] + 1]:
                cross_indices = cross_indices[1:]
            durations = []
            for i in range(0, len(cross_indices), 2):
                if i + 1 < len(cross_indices):
                    duration = time_ds[cross_indices[i + 1]] - time_ds[cross_indices[i]]
                    durations.append(duration)
            # Converting the list to a numpy array for future calculations
            durations = np.array(durations)
            durations = [x for x in durations if x > 0.1]
            # Define the bins
            bins = np.arange(0.1, 0.8, 0.1)
            bins = np.append(bins, np.inf)
            counts, _ = np.histogram(durations, bins)
            percentages = (counts / len(durations)) * 100  # Updated calculation
            average = np.mean(durations)
            return percentages, average
    
    events = spike_indices / fs
    spike_oscillations_theta_wave, spike_oscillations_theta_times, theta_wave_power = waveform_spiketrain_oscillation(events, 4, 8)
    spike_oscillations_alpha_wave, spike_oscillations_alpha_times, alpha_wave_power = waveform_spiketrain_oscillation(events, 8, 12)
    spike_oscillations_low_beta_wave, spike_oscillations_low_beta_times, low_beta_wave_power = waveform_spiketrain_oscillation(events,
                                                                                                    12, 21)
    spike_oscillations_high_beta_wave, spike_oscillations_high_beta_times, high_beta_wave_power = waveform_spiketrain_oscillation(
        events, 21, 30)
    
    spiketrain_burst_threshold = burst_threshold(events, data_type='spiketrain')

    theta_wave_envelope = get_waveform_envelope(spike_oscillations_theta_wave)
    alpha_wave_envelope = get_waveform_envelope(spike_oscillations_alpha_wave)
    low_beta_wave_envelope = get_waveform_envelope(spike_oscillations_low_beta_wave)
    high_beta_wave_envelope = get_waveform_envelope(spike_oscillations_high_beta_wave)
    
    theta = \
    burst_features(spike_oscillations_theta_times, theta_wave_envelope, spiketrain_burst_threshold)[1]
    alpha = \
    burst_features(spike_oscillations_alpha_times, alpha_wave_envelope, spiketrain_burst_threshold)[1]
    lowbeta = \
    burst_features(spike_oscillations_low_beta_times, low_beta_wave_envelope, spiketrain_burst_threshold)[1]
    highbeta = \
    burst_features(spike_oscillations_high_beta_times, high_beta_wave_envelope, spiketrain_burst_threshold)[1]

    return theta, alpha, lowbeta, highbeta


def get_lfp_metrics(raw_signal, fs, fs_lfp=250):

    def get_LFP_data(raw_data, fs, fs_lfp):
        return ds.run(raw_data, fs, fs_lfp)
    
    def get_LFP_power(raw_data_lfp, fs_lfp, fFrom, fTo):
        filtered_data = ff.butter(raw_data_lfp, fFrom, fTo, fs_lfp, order=8, zero_phase=True)
        # Apply Welch's method to estimate PSD
        f_lfp, Pxx_lfp = sps.welch(raw_data_lfp, fs_lfp, nperseg=256)
        # Normalize the PSD as relative power (for frequencies > 80 Hz)
        normalized_freq = (f_lfp >= 45) & (f_lfp <= 55)
        total_power_high_freq = np.sum(Pxx_lfp[normalized_freq])
        Pxx_lfp_normalized = Pxx_lfp / total_power_high_freq
        freq_range = (f_lfp >= fFrom) & (f_lfp <= fTo)
        power = 10 * math.log10(np.sum(Pxx_lfp_normalized[freq_range]))
        return filtered_data, power
    
    raw_data_lfp = get_LFP_data(raw_signal, fs, fs_lfp)
    raw_data_lfp = (raw_data_lfp - np.mean(raw_data_lfp)) / np.std(raw_data_lfp)
    
    alpha_wave, alpha_power = get_LFP_power(raw_data_lfp, fs_lfp, 8, 12)
    theta_wave, theta_power = get_LFP_power(raw_data_lfp, fs_lfp, 4, 8)
    lowbeta_wave, lowbeta_power = get_LFP_power(raw_data_lfp, fs_lfp, 12, 21)
    highbeta_wave, highbeta_power = get_LFP_power(raw_data_lfp, fs_lfp, 21, 30)
    
    psd_freqs, psd_power = sps.welch(raw_data_lfp, fs_lfp, nperseg=256)

    return (theta_power, alpha_power, lowbeta_power, highbeta_power), \
    (theta_wave, alpha_wave, lowbeta_wave, highbeta_wave), (psd_freqs, psd_power)