import numpy as np
from sonpy import lib as sp
from sklearn.mixture import GaussianMixture
from scipy.stats import median_abs_deviation
from scipy.signal import butter, filtfilt, hilbert, convolve, welch, argrelextrema, find_peaks
from astropy.timeseries import LombScargle
import math
import warnings
import finnpy.filters.frequency as ff
import finnpy.basic.downsampling as ds

warnings.filterwarnings("ignore")


class AnalyzeMER:
    def __init__(self):
        self.npz_data = None
        self.DataReadFunctions = {
            sp.DataType.Adc: sp.SonFile.ReadFloats,
            sp.DataType.EventFall: sp.SonFile.ReadEvents,
            sp.DataType.EventRise: sp.SonFile.ReadEvents,
            sp.DataType.EventBoth: sp.SonFile.ReadEvents,
            sp.DataType.Marker: sp.SonFile.ReadMarkers,
            sp.DataType.AdcMark: sp.SonFile.ReadWaveMarks,
            sp.DataType.RealMark: sp.SonFile.ReadRealMarks,
            sp.DataType.TextMark: sp.SonFile.ReadTextMarks,
            sp.DataType.RealWave: sp.SonFile.ReadFloats
        }

    @staticmethod
    def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx

    @staticmethod
    def box_pir(l):
        lengths = [i for i in map(len, l)]
        shape = (len(l), max(lengths))
        a = np.full(shape, np.nan)
        for i, r in enumerate(l):
            a[i, :lengths[i]] = r
        return a

    @staticmethod
    def get_spike_matrix(filtered_data, spiketrain, time):
        spike_matrix = []

        if spiketrain[0] <= 0.001:
            spiketrain = spiketrain[1, :]

        if np.abs(spiketrain[-1] - time[-1]) <= 0.002:
            spiketrain = spiketrain[:-1]

        for t in spiketrain:
            spike_matrix.append(
                filtered_data[AnalyzeMER.find_nearest(time, t - 0.001):AnalyzeMER.find_nearest(time, t + 0.002)])

        spike_matrix = np.squeeze(spike_matrix)

        if spike_matrix.ndim == 1:
            spike_matrix = list(spike_matrix)
            spike_matrix = AnalyzeMER.box_pir(spike_matrix)
            spike_matrix = np.delete(spike_matrix, -1, axis=1)

        return spike_matrix

    @staticmethod
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

    @staticmethod
    def spiketrain_to_binary(spiketimes, sampling_rate):
        num_timesteps = int(np.ceil(np.max(spiketimes) * sampling_rate))
        binary_sequence = np.zeros(num_timesteps)
        spike_indices = (spiketimes * sampling_rate).astype(int)
        valid_spike_indices = spike_indices[spike_indices < num_timesteps]
        binary_sequence[valid_spike_indices] = 1
        return binary_sequence

    def load_npz(self, file_path, list_files_before_load=False):
        if list_files_before_load:
            with np.load(file_path, allow_pickle=True) as data:
                print("Files in the .npz file: ", data.files)
        self.npz_data = np.load(file_path, allow_pickle=True)

    def get_npz_data(self, file_name):
        if self.npz_data is None:
            print("No .npz file loaded yet!")
        else:
            try:
                return self.npz_data[file_name]
            except KeyError:
                print(f"No such file: {file_name} in the loaded .npz file!")

    def get_spike2_file_info(self, FilePath):
        # Open file
        MyFile = sp.SonFile(FilePath, True)

        if MyFile.GetOpenError() != 0:
            print('Error opening file:', sp.GetErrorString(MyFile.GetOpenError()))

        # Initialize variable to None
        firstAdcChannel = None
        firstEventFallChannel = None

        # Loop through channels and print info
        for i in range(MyFile.MaxChannels()):
            if MyFile.ChannelType(i) != sp.DataType.Off:
                print(
                    f'Channel: {i + 1}, Type: {MyFile.ChannelType(i)}, Title: {MyFile.GetChannelTitle(i)}, Comment: {MyFile.GetChannelComment(i)}, Sampling Frequency: {round(1 / (MyFile.ChannelDivide(i) * MyFile.GetTimeBase()))} Hz')

                # Check if channel type is DataType.Adc and it has not been found yet
                if MyFile.ChannelType(i) == sp.DataType.Adc and firstAdcChannel is None:
                    firstAdcChannel = i + 1

                # Check if channel type is DataType.EventFall and it has not been found yet
                if MyFile.ChannelType(i) == sp.DataType.EventFall and firstEventFallChannel is None:
                    firstEventFallChannel = i + 1

        # Return first "Type: DataType.Adc" and first "Type: DataType.EventFall" channel indices
        return firstAdcChannel, firstEventFallChannel

    def load_spike2_data(self, FilePath, selected_channel):
        # Open file
        MyFile = sp.SonFile(FilePath, True)

        if MyFile.GetOpenError() != 0:
            print('Error opening file:', sp.GetErrorString(MyFile.GetOpenError()))

        # Check if selected channel is valid
        if MyFile.ChannelType(selected_channel - 1) == sp.DataType.Off:
            print(f'Channel {selected_channel} is turned off.')
            return None

        # Retrieve and return data for selected channel
        return {
            'channel': selected_channel,
            'fs': round(1 / (MyFile.ChannelDivide(selected_channel - 1) * MyFile.GetTimeBase())),
            'data': np.array(
                self.DataReadFunctions[MyFile.ChannelType(selected_channel - 1)](MyFile, selected_channel - 1,
                                                                                 int(MyFile.ChannelMaxTime(
                                                                                     selected_channel - 1) / MyFile.ChannelDivide(
                                                                                     selected_channel - 1)), 0,
                                                                                 MyFile.ChannelMaxTime(
                                                                                     selected_channel - 1)))
        }

    def get_LFP_data(self, raw_data, fs, fs_lfp):
        return ds.run(raw_data, fs, fs_lfp)

    def get_LFP_power(self, raw_data_lfp, fs_lfp, fFrom, fTo):

        filtered_data = ff.butter(raw_data_lfp, fFrom, fTo, fs_lfp, order=8, zero_phase=True)

        # Apply Welch's method to estimate PSD
        f_lfp, Pxx_lfp = welch(raw_data_lfp, fs_lfp, nperseg=256)

        # Normalize the PSD as relative power (for frequencies > 80 Hz)
        normalized_freq = (f_lfp >= 45) & (f_lfp <= 55)
        total_power_high_freq = np.sum(Pxx_lfp[normalized_freq])
        Pxx_lfp_normalized = Pxx_lfp / total_power_high_freq

        freq_range = (f_lfp >= fFrom) & (f_lfp <= fTo)

        power = 10 * math.log10(np.sum(Pxx_lfp_normalized[freq_range]))

        return filtered_data, power

    def get_FR(self, spiketrain):
        isi = np.diff(spiketrain)
        return 1 / isi.mean()

    def get_BI(self, spiketrain):
        isi = np.diff(spiketrain)
        X = np.ravel(isi).reshape(-1, 1)
        M_best = GaussianMixture(n_components=2, covariance_type="spherical", random_state=1).fit(X)
        means = M_best.means_
        return np.max(means) / np.min(means)

    def get_CV(self, spiketrain):
        isi = np.diff(spiketrain)
        return median_abs_deviation(isi) / np.median(isi)

    def get_segment_duration(self, raw_data, fs):
        time = np.arange(0, len(raw_data)) / fs
        return time[-1]

    def get_snr(self, raw_data, fs, spiketrain, time):

        filtered_data = ff.butter(raw_data, 300, 3000, fs, order=4, zero_phase=True)

        waveforms = AnalyzeMER.get_spike_matrix(filtered_data, spiketrain, time)
        mean_waveform = waveforms.mean(axis=0)
        std_waveform = waveforms.std(axis=0).mean(axis=-1)
        peak_range = mean_waveform.max(axis=-1) - mean_waveform.min(axis=-1)
        noise = 2 * std_waveform

        return peak_range / noise

    def kaneoke_oscillation_power(self, spike_train, low, high, bin_size, max_lag):

        autocorr, autocorr_lag = AnalyzeMER.autocorrelation(spike_train, bin_size, max_lag)

        freqs = np.linspace(4, 30, num=1000)

        ls = LombScargle(autocorr_lag, autocorr)
        power = ls.power(freqs)

        # Find the maximum peak in the specified frequency band
        low_freq, high_freq = (low, high)
        idx = np.where((freqs >= low_freq) & (freqs < high_freq))[0]
        peaks, _ = find_peaks(power[idx])

        if peaks.size > 0:
            max_peak_idx = idx[peaks[np.argmax(power[idx][peaks])]]
            max_peak_power = power[max_peak_idx]
            return 10 * math.log10(max_peak_power)
        else:
            return 10 * math.log10(max(power[idx]))

    def waveform_spiketrain_oscillation(self, spike_train, low_freq, high_freq):

        binary_sequence = AnalyzeMER.spiketrain_to_binary(spike_train, 500)

        spiketrain_waveform = ff.butter(binary_sequence, low_freq, high_freq, 500, order=4, zero_phase=True)

        spiketrain_waveform = spiketrain_waveform - np.mean(spiketrain_waveform)
        freqs, psd = welch(spiketrain_waveform, 500, nperseg=500 * 2, scaling='density')
        power = np.max(psd)
        time = np.arange(0, len(spiketrain_waveform)) / 500

        return spiketrain_waveform, time, 10 * math.log10(power)

    def get_waveform_envelope(self, waveform, smooth=False):

        burst = np.abs(hilbert(waveform))

        if smooth > 0:
            window_size = smooth
            window = np.ones(window_size) / window_size
            burst = convolve(burst, window, mode='same')
            return burst
        else:
            return burst

    def burst_threshold(self, data, fs_lfp=None, data_type='lfp'):
        if data_type == 'lfp':

            assert fs_lfp is not None, 'fs_lfp must be provided for lfp data'
            gamma_envelopes = [np.abs(hilbert(ff.butter(data, f1, f2, fs_lfp, order=8, zero_phase=True)))
                               for f1, f2 in zip(range(45, 50), range(51, 56))]
        elif data_type == 'spiketrain':
            gamma_envelopes = [np.abs(hilbert(self.waveform_spiketrain_oscillation(data, f1, f2)[0]))
                               for f1, f2 in zip(range(45, 50), range(51, 56))]
        else:
            raise ValueError('data_type must be either "lfp" or "spiketrain"')

        median_values = [np.median(data[argrelextrema(data, np.less)]) for data in gamma_envelopes]

        return 4 * np.mean(median_values)

    def burst_features(self, time_ds, burst, burst_threshold):
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

