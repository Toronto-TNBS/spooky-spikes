from matplotlib import pyplot as plt
import numpy as np
from scipy.signal import filtfilt, iirfilter, butter, sosfilt, find_peaks, welch, peak_widths
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import pandas as pd
from astropy.timeseries import LombScargle
import math
from sklearn.mixture import GaussianMixture
from scipy.stats import median_abs_deviation
import neo
import heapq as pq
import warnings

plt.style.use('ggplot')


class SpikeAnalysis:

    def __init__(self):
        self.main_times = []
        self.main_magnitudes = []
        self.main_events = []
        self.main_threshold_factor = 0
        self.main_fs = 0
        self.main_threshold = 0
        self.main_disp_threshold = False
        self.main_event_peaks_disp = False
        self.main_spiketrain_disp = False
        self.main_event_peak_times = []
        self.main_event_peak_mags = []
        self.reds = []
        self.blues = []
        self.oranges = []
        self.cyans = []
        self.purples = []
        self.browns = []
        self.sorted_events = [self.reds, self.blues, self.oranges, self.cyans, self.purples, self.browns]
        self.selected_cluster = None
        self.selected_cluster_peaks = []

        # Counts how many times threshold set button is pressed. Used to fix threshold factor display.
        self.num_threshold_sets = 0

        self.segment_inverted = False

        self.spike_matrix = [[]]
        self.peak_indices = []
        self.spikesorting_plot_disp = False
        self.number_desired_clusters = 0
        self.labels = []

        self.psd_good_file = False
        self.lfp_psd_freqs = []
        self.lfp_psd_power = []
        self.lfp_theta_wave = []
        self.lfp_alpha_wave = []
        self.lfp_low_beta_wave = []
        self.lfp_high_beta_wave = []
        self.psd_plot_xlim = [0,100]

        self.spike_oscillations_theta_wave = []
        self.spike_oscillations_alpha_wave = []
        self.spike_oscillations_low_beta_wave = []
        self.spike_oscillations_high_beta = []
        self.spike_oscillations_theta_times = []
        self.spike_oscillations_alpha_times = []
        self.spike_oscillations_low_beta_times = []
        self.spike_oscillations_high_beta_times = []

        self.lag_time = 0.5
        self.time_interval = 0.01
        self.frequencies_autocorr = []
        self.psd_autocorr = []
        self.binned_time = []
        self.autocorr = []
        self.main_threshold_set = False

        self.isi = np.array([])
        self.silhouette = 0
        self.coloured_spike_matrix = [[]]
        self.coloured_peak_indices = []
        self.coloured_frequencies_autocorr = []
        self.coloured_psd_autocorr = []
        self.coloured_binned_time = []
        self.coloured_autocorr = []

        self.features_spiketrain_indices = []    # Will be updated depending on spike-selection method (i.e. thresh vs sort)


    def main_plot(self):
        fig, ax = plt.subplots(ncols=1, nrows=2, gridspec_kw={'height_ratios': [1, 8]}, figsize=(8, 5.15), sharex=True)
        ax[0].set_ylabel('Spike Train')
        ax[0].spines[['top', 'right', 'bottom']].set_visible(False)
        ax[0].set_yticks([])
        ax[1].set_xlabel('Time (s)')
        ax[1].set_ylabel('Magnitude (V)')
        ax[1].spines['top'].set_visible(False)

        # This allows for coloured events to be saved to variables without having to check spiketrain.
        if self.spikesorting_plot_disp:
            self.get_coloured_events()

        if self.spikesorting_plot_disp:
            spike_window = int(10 ** -3 * self.main_fs)
            labels = [str(i) for i in self.labels]
            colour_list = {
                '0': 'tab:red',
                '1': 'tab:blue',
                '2': 'tab:orange',
                '3': 'tab:cyan',
                '4': 'tab:purple',
                '5': 'tab:brown'
            }

            times_matrix = np.empty([1, spike_window * 2])

            for i in self.peak_indices:
                time_window = self.main_times[int(i - spike_window):int(i + spike_window)]
                times_matrix = np.append(times_matrix, time_window.reshape(1, spike_window * 2), axis=0)
            times_matrix = times_matrix[1:len(times_matrix)]

            ax[1].plot(self.main_times, self.main_magnitudes, 'k')
            peaks = []
            red_peaks = []
            blue_peaks = []
            orange_peaks = []
            cyan_peaks = []
            purple_peaks = []
            brown_peaks = []
            for i in range(len(times_matrix)):
                if labels[i] == '0':
                    red_peaks.append(self.main_magnitudes[int(self.peak_indices[i])])
                    peaks.append(self.main_magnitudes[int(self.peak_indices[i])])
                elif labels[i] == '1':
                    blue_peaks.append(self.main_magnitudes[int(self.peak_indices[i])])
                    peaks.append(self.main_magnitudes[int(self.peak_indices[i])])
                elif labels[i] == '2':
                    orange_peaks.append(self.main_magnitudes[int(self.peak_indices[i])])
                    peaks.append(self.main_magnitudes[int(self.peak_indices[i])])
                elif labels[i] == '3':
                    cyan_peaks.append(self.main_magnitudes[int(self.peak_indices[i])])
                    peaks.append(self.main_magnitudes[int(self.peak_indices[i])])
                elif labels[i] == '4':
                    purple_peaks.append(self.main_magnitudes[int(self.peak_indices[i])])
                    peaks.append(self.main_magnitudes[int(self.peak_indices[i])])
                elif labels[i] == '5':
                    brown_peaks.append(self.main_magnitudes[int(self.peak_indices[i])])
                    peaks.append(self.main_magnitudes[int(self.peak_indices[i])])
            means = [np.mean(red_peaks), np.mean(blue_peaks), np.mean(orange_peaks), np.mean(cyan_peaks), np.mean(purple_peaks), np.mean(brown_peaks)]
            if self.selected_cluster == None:
                if self.segment_inverted:
                    plotted_peaks = means.index(max(means))
                else:
                    plotted_peaks = means.index(min(means))
                self.selected_cluster = colour_list[str(plotted_peaks)].split(':')[-1]
                print(f'Selected cluster: {self.selected_cluster}')
            else:
                print(f'Selected cluster: {self.selected_cluster}')
                plotted_peaks = [key for key, value in colour_list.items() if value == f'tab:{self.selected_cluster}']
                print(f'Plotted_peaks: {plotted_peaks}')
                plotted_peaks = int(plotted_peaks[0])
            for i in range(len(times_matrix)):
                if labels[i] == str(plotted_peaks):
                    ax[1].plot(times_matrix[i], self.spike_matrix[i], colour_list[str(plotted_peaks)])
                    ax[1].plot([self.main_times[int(self.peak_indices[i])]], [peaks[i]], colour_list[labels[i]], marker='.')
                else:
                    ax[1].plot([self.main_times[int(self.peak_indices[i])]], [peaks[i]], colour_list[labels[i]], marker='.')


        else:
            ax[1].plot(self.main_times, self.main_magnitudes, 'cornflowerblue')
        if self.main_event_peaks_disp and not self.spikesorting_plot_disp:
            events = self.get_event_peaks(
                times=self.main_times,
                magnitudes=self.main_magnitudes,
                threshold=self.main_threshold
            )
            self.main_event_peak_times = events[0]
            self.main_event_peak_mags = events[1]
            ax[1].plot(self.main_event_peak_times, self.main_event_peak_mags, '.', linewidth=3)
        if self.main_disp_threshold:
            threshold_bar = [self.main_threshold for i in range(len(self.main_times))]
            ax[1].plot(self.main_times, threshold_bar, 'k--', linewidth=0.5)
        if self.main_spiketrain_disp and self.spikesorting_plot_disp:
            colour_list = ['tab:red', 'tab:blue', 'tab:orange', 'tab:cyan', 'tab:purple', 'tab:brown']
            events, indices = self.get_coloured_events()
            events = events[colour_list.index(f'tab:{self.selected_cluster}')]
            ax[0].eventplot(events, colors=f'tab:{self.selected_cluster}')
        elif self.main_spiketrain_disp and self.main_threshold != 0:
            events = self.get_event_peaks(
                times=self.main_times,
                magnitudes=self.main_magnitudes,
                threshold=self.main_threshold
            )
            self.main_event_peak_times = events[0]
            ax[0].eventplot(self.main_event_peak_times)

        return fig


    def spikesorting_plot(self):
        fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(10, 6.5))
        ax.set_xlabel('Principal Component 1')
        ax.set_ylabel('Principal Component 2')

        if self.spikesorting_plot_disp and self.main_threshold_set:

            ## ADDED AUTOMATIC FUNCTIONALITY
            ss = []
            for c in range(2, 7):
                features_temp, label_temp = self.spike_sorting(self.spike_matrix, c)
                ss.append(silhouette_score(features_temp, label_temp))

            clusters_auto = np.array(ss).argmax() + 2
            if self.number_desired_clusters == 0:
                clusters = clusters_auto
                self.number_desired_clusters = clusters
            else:
                clusters = self.number_desired_clusters

            features, labels = self.spike_sorting(self.spike_matrix, clusters)
            self.silhouette = silhouette_score(features, labels)

            labels = [str(i) for i in labels]
            df = pd.DataFrame({'colours': labels})
            colour_list = {
                '0': 'tab:red',
                '1': 'tab:blue',
                '2': 'tab:orange',
                '3': 'tab:cyan',
                '4': 'tab:purple',
                '5': 'tab:brown'
            }

            ax.scatter(features[:, 0], features[:, 1], c=df['colours'].map(colour_list))

        else:
            ax.plot()

        return fig


    def psd_plot(self):
        fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(10, 7.5))
        ax.set_xlabel('Frequency')
        ax.set_ylabel('Power')
        if self.psd_good_file:
            ax.plot(self.psd_frequencies, self.psd_power, 'k')
        else:
            ax.plot()
        ax.set_xlim(self.psd_plot_xlim)
        return fig


    def oscillations_plot(self):
        fig, ax = plt.subplots(ncols=1, nrows=2, figsize=(10, 8.5))
        plt.subplots_adjust(hspace=0.3)
        # ax[0].stem(self.binned_time, self.autocorr, markerfmt=" ")
        ax[0].set_title('Spiketrain Autocorrelation Function')
        ax[0].set_xlabel('Lag')
        ax[0].set_ylabel('Autocorrelation')

        ax[1].set_xlim([0, 50])
        ax[1].set_title('Lomb-Scargle Periodogram')
        ax[1].set_xlabel('Frequency')
        ax[1].set_ylabel('Power')

        if self.main_threshold_set and self.lag_time != 0 and self.time_interval != 0:
            ax[0].set_ylim([min(self.autocorr) - 0.1, max(self.autocorr) + 0.1])
            ax[0].plot(self.binned_time, self.autocorr)
            ax[1].plot(self.frequencies_autocorr, self.psd_autocorr)
        else:
            ax[0].plot()
            ax[1].plot()

        return fig


    def get_coloured_events(self):
        events = [self.main_times[int(i)] for i in self.peak_indices]
        index = 0
        self.reds = []
        self.reds_indices = []
        self.blues = []
        self.blues_indices = []
        self.oranges = []
        self.oranges_indices = []
        self.cyans = []
        self.cyans_indices = []
        self.purples = []
        self.purples_indices = []
        self.browns = []
        self.browns_indices = []
        for i in self.labels:
            if i == 0:
                self.reds_indices.append(self.peak_indices[index])
                self.reds.append(events[index])
            elif i == 1:
                self.blues_indices.append(self.peak_indices[index])
                self.blues.append(events[index])
            elif i == 2:
                self.oranges_indices.append(self.peak_indices[index])
                self.oranges.append(events[index])
            elif i == 3:
                self.cyans_indices.append(self.peak_indices[index])
                self.cyans.append(events[index])
            elif i == 4:
                self.purples_indices.append(self.peak_indices[index])
                self.purples.append(events[index])
            elif i == 5:
                self.browns_indices.append(self.peak_indices[index])
                self.browns.append(events[index])
            index += 1
        events = [self.reds, self.blues, self.oranges, self.cyans, self.purples, self.browns]
        events_indices = [self.reds_indices, self.blues_indices, self.oranges_indices, self.cyans_indices, self.purples_indices, self.browns_indices]
        return events, events_indices


    def filter_wave(self, wave, mincutoff, maxcutoff, fs):
        if mincutoff == None and maxcutoff == None:
            return wave
        elif mincutoff == None:
            b, a = iirfilter(4, int(maxcutoff), btype='lowpass', ftype='butter', fs=float(fs))
        elif maxcutoff == None:
            b, a = iirfilter(4, int(maxcutoff), btype='highpass', ftype='butter', fs=float(fs))
        else:
            b, a = iirfilter(4, [int(mincutoff), int(maxcutoff)], btype='bandpass', ftype='butter', fs=float(fs))
        filt_signal = filtfilt(b, a, wave)
        return np.array(filt_signal)


    def get_main_threshold(self, magnitudes, factor):
        if self.segment_inverted:
            threshold = factor * np.median(np.abs(magnitudes))/0.6745
        else:
            threshold = -factor * np.median(np.abs(magnitudes))/0.6745
        return threshold


    def get_event_peaks(self, times, magnitudes, threshold):
        min_dist = round(1.5/1000 * self.main_fs)
        if self.segment_inverted:
            peaks = find_peaks(magnitudes, height=threshold, distance=min_dist)
        else:
            peaks = find_peaks(-np.array(magnitudes), height=-threshold, distance=min_dist)
        peaks_index = peaks[0]
        return np.array([times[peaks_index], magnitudes[peaks_index], peaks_index])


    def get_spike_windows(self, magnitudes, event_indices, spike_window):
        event_indices = event_indices[
            (event_indices > spike_window) & (event_indices < event_indices[-1] - spike_window)]
        spikes_matrix = np.empty([1, spike_window * 2])
        for k in event_indices:
            temp_spike = magnitudes[int(k - spike_window):int(k + spike_window)]
            spikes_matrix = np.append(spikes_matrix, temp_spike.reshape(1, spike_window * 2), axis=0)

        spikes_matrix = np.delete(spikes_matrix, 0, axis=0)

        return event_indices, spikes_matrix


    def spike_sorting(self, spike_matrix, num_clusters):
        pca = PCA(n_components=2)
        features = pca.fit_transform(spike_matrix)

        kmeans = KMeans(n_clusters=num_clusters, random_state=0).fit(features)
        labels = kmeans.labels_

        self.labels = labels

        return features, labels


    def get_psd(self, magnitudes, fs):
        # F = 10
        # freqs, psd = welch(magnitudes, nfft=F * fs, fs=fs, nperseg=fs)
        freqs, power = welch(magnitudes, fs, nperseg=256)
        return freqs, power


    def get_spike_oscillations(self, magnitudes, event_indices, fs, lag_time, time_interval):
        event_indices = np.array(event_indices, dtype='int')
        spiketrain_index = np.array(event_indices, dtype='int')
        spike_data = np.zeros(magnitudes.shape, dtype='int')
        spike_data[event_indices] = 1

        event_indices = event_indices / fs

        t = np.arange(0, len(magnitudes)) / fs
        t_stop = t[-1]

        spiketrain_trimmed = event_indices[(event_indices > lag_time) & (event_indices < t_stop - lag_time)]
        spiketrain_index_trimmed = spiketrain_index[(event_indices > lag_time) & (event_indices < t_stop - lag_time)]
        lag = math.floor(lag_time * fs)
        lead = lag
        time_interval = math.floor(time_interval * fs)
        print(f'Lag: {lag}')
        print(f'Interval: {time_interval}')
        bins = np.arange(0, lag, time_interval)

        # LAGGING

        spike_data_matrix_lag = []

        for s in spiketrain_index_trimmed:
            temp = spike_data[s - lag:s]
            spike_data_matrix_lag.append(temp)

        spike_data_matrix_lag = np.stack(spike_data_matrix_lag, axis=0)

        binned_spike_data_lag = []
        binned_spike_data_matrix_lag = []
        for r in range(len(spike_data_matrix_lag)):
            for b in bins:
                binned_spike_data_lag.append(sum(spike_data_matrix_lag[r][b:b + time_interval]))

            binned_spike_data_matrix_lag.append(np.array(binned_spike_data_lag))
            binned_spike_data_lag = []

        binned_spike_data_matrix_lag = np.stack(binned_spike_data_matrix_lag, axis=0)

        autocorr_lag = np.sum(binned_spike_data_matrix_lag, axis=0)

        # LEADING

        spike_data_matrix_lead = []

        for s in spiketrain_index_trimmed:
            temp = spike_data[s:s + lead]
            spike_data_matrix_lead.append(temp)

        spike_data_matrix_lead = np.stack(spike_data_matrix_lead, axis=0)

        binned_spike_data_lead = []
        binned_spike_data_matrix_lead = []
        for r in range(len(spike_data_matrix_lead)):
            for b in bins:
                binned_spike_data_lead.append(sum(spike_data_matrix_lead[r][b:b + time_interval]))

            binned_spike_data_matrix_lead.append(np.array(binned_spike_data_lead))
            binned_spike_data_lead = []

        binned_spike_data_matrix_lead = np.stack(binned_spike_data_matrix_lead, axis=0)

        autocorr_lead = np.sum(binned_spike_data_matrix_lead, axis=0)

        if float(fs).is_integer() == False:
            autocorr_lag = np.delete(autocorr_lag, -1)
            autocorr_lead = np.delete(autocorr_lead, -1)

        autocorr_lead[0] = np.mean([autocorr_lag[-1], autocorr_lead[1]])

        autocorr = np.concatenate([autocorr_lag, autocorr_lead])
        autocorr = (autocorr / len(spike_data_matrix_lead)) / time_interval
        autocorr = autocorr / np.mean(autocorr)
        autocorr = autocorr - np.mean(autocorr)

        binned_time = np.linspace(-lag_time, lag_time, len(autocorr))

        # freqs_autocorr, psd_autocorr = welch(autocorr, fs = 1/0.01)
        freqs_autocorr, psd_autocorr = LombScargle(binned_time, autocorr).autopower()
        print('done getting spike oscillations')
        return freqs_autocorr, psd_autocorr, binned_time, autocorr


    def get_snr(self):
        if self.selected_cluster == None:
            waveforms = self.spike_matrix
        else:
            waveforms = self.coloured_spike_matrix
        # el.waveform_features.waveform_snr(spike_matrix)
        # ----- Replacement from source code:

        if isinstance(waveforms, neo.SpikeTrain):
            warnings.warn("spiketrain input is deprecated; pass "
                          "'spiketrain.waveforms' directly.", DeprecationWarning)
            waveforms = waveforms.waveforms
        # asarray removes quantities, if present
        waveforms = np.squeeze(np.asarray(waveforms))
        # average over all waveforms for each bin
        mean_waveform = waveforms.mean(axis=0)
        # standard deviation over all waveforms over all bins
        std_waveform = waveforms.std(axis=0).mean(axis=-1)

        # peak to trough voltage signal
        peak_range = mean_waveform.max(axis=-1) - mean_waveform.min(axis=-1)
        # noise
        noise = 2 * std_waveform

        snr = peak_range / noise
        if np.isnan(snr).any():
            warnings.warn('The waveforms noise was evaluated to 0. Returning NaN')

        return snr


    def get_percent_isi_violations(self):
        if self.selected_cluster == None:
            peak_indices = self.peak_indices
        else:
            peak_indices = self.coloured_peak_indices
        spiketrain = peak_indices / self.main_fs
        axis = -1
        # self.isi = el.statistics.isi(spiketrain)
        # ---- Replacement from source code:
        if isinstance(spiketrain, neo.SpikeTrain):
            intervals = np.diff(spiketrain.magnitude, axis=axis)
            # np.diff makes a copy
            intervals = pq.Quantity(intervals, units=spiketrain.units, copy=False)
        else:
            intervals = np.diff(spiketrain, axis=axis)
        if (intervals < 0).any():
            warnings.warn("ISI evaluated to negative values. "
                          "Please sort the input array.")

        self.isi = intervals

        return (sum(self.isi < 1 / 1000) / len(self.isi)) * 100


    def get_firing_rate(self):
        return 1 / self.isi.mean()


    def get_burst_index(self):
        X = np.ravel(self.isi).reshape(-1, 1)
        M_best = GaussianMixture(n_components=2, covariance_type="spherical").fit(X)
        means = M_best.means_
        return np.max(means) / np.min(means)


    def get_variation_coefficient(self):
        return median_abs_deviation(self.isi) / np.median(self.isi)


    def find_nearest(self, array, value):
        idx = np.searchsorted(array, value, side="left")
        if idx > 0 and (idx == len(array) or math.fabs(value - array[idx - 1]) < math.fabs(value - array[idx])):
            return idx
        else:
            return idx


    def get_wave_powers(self):
        theta_peaks, _ = find_peaks(self.psd_autocorr[self.find_nearest(self.frequencies_autocorr, 4):self.find_nearest(self.frequencies_autocorr, 8)])
        alpha_peaks, _ = find_peaks(self.psd_autocorr[self.find_nearest(self.frequencies_autocorr, 8):self.find_nearest(self.frequencies_autocorr, 12)])
        low_beta_peaks, _ = find_peaks(self.psd_autocorr[self.find_nearest(self.frequencies_autocorr, 12):self.find_nearest(self.frequencies_autocorr, 21)])
        high_beta_peaks, _ = find_peaks(self.psd_autocorr[self.find_nearest(self.frequencies_autocorr, 21):self.find_nearest(self.frequencies_autocorr, 30)])

        if (len(theta_peaks) == 0) or (len(alpha_peaks) == 0) or (
                len(low_beta_peaks) == 0) or (len(high_beta_peaks) == 0):
            theta_freq = -100
            alpha_freq = -100
            low_beta_freq = -100
            high_beta_freq = -100

            theta_power = -100
            alpha_power = -100
            low_beta_power = -100
            high_beta_power = -100

        else:
            theta_peaks = theta_peaks + self.find_nearest(self.frequencies_autocorr, 4)
            alpha_peaks = alpha_peaks + self.find_nearest(self.frequencies_autocorr, 8)
            low_beta_peaks = low_beta_peaks + self.find_nearest(self.frequencies_autocorr, 12)
            high_beta_peaks = high_beta_peaks + self.find_nearest(self.frequencies_autocorr, 21)

            theta_peak_width = peak_widths(self.psd_autocorr, theta_peaks, rel_height=1)[0][
                self.psd_autocorr[theta_peaks].argmax()]
            alpha_peak_width = peak_widths(self.psd_autocorr, alpha_peaks, rel_height=1)[0][
                self.psd_autocorr[alpha_peaks].argmax()]
            low_beta_peak_width = peak_widths(self.psd_autocorr, low_beta_peaks, rel_height=1)[0][
                self.psd_autocorr[low_beta_peaks].argmax()]
            high_beta_peak_width = peak_widths(self.psd_autocorr, high_beta_peaks, rel_height=1)[0][
                self.psd_autocorr[high_beta_peaks].argmax()]

            # theta_peak = self.psd_autocorr[theta_peaks[self.psd_autocorr[theta_peaks].argmax()]]
            theta_freq = self.frequencies_autocorr[theta_peaks[self.psd_autocorr[theta_peaks].argmax()]]
            theta_power = sum(self.psd_autocorr[
                              self.find_nearest(self.frequencies_autocorr, theta_freq - theta_peak_width * 0.1 / 2):self.find_nearest(
                                  self.frequencies_autocorr, theta_freq + theta_peak_width * 0.1 / 2)])

            # alpha_peak = self.psd_autocorr[alpha_peaks[self.psd_autocorr[alpha_peaks].argmax()]]
            alpha_freq = self.frequencies_autocorr[alpha_peaks[self.psd_autocorr[alpha_peaks].argmax()]]
            alpha_power = sum(self.psd_autocorr[
                              self.find_nearest(self.frequencies_autocorr, alpha_freq - alpha_peak_width * 0.1 / 2):self.find_nearest(
                                  self.frequencies_autocorr, alpha_freq + alpha_peak_width * 0.1 / 2)])

            # low_beta_peak = self.psd_autocorr[low_beta_peaks[self.psd_autocorr[low_beta_peaks].argmax()]]
            low_beta_freq = self.frequencies_autocorr[low_beta_peaks[self.psd_autocorr[low_beta_peaks].argmax()]]
            low_beta_power = sum(self.psd_autocorr[self.find_nearest(self.frequencies_autocorr,
                                                           low_beta_freq - low_beta_peak_width * 0.1 / 2):self.find_nearest(
                self.frequencies_autocorr, low_beta_freq + low_beta_peak_width * 0.1 / 2)])

            # high_beta_peak = self.psd_autocorr[high_beta_peaks[self.psd_autocorr[high_beta_peaks].argmax()]]
            high_beta_freq = self.frequencies_autocorr[high_beta_peaks[self.psd_autocorr[high_beta_peaks].argmax()]]
            high_beta_power = sum(self.psd_autocorr[self.find_nearest(self.frequencies_autocorr,
                                                            high_beta_freq - high_beta_peak_width * 0.1 / 2):self.find_nearest(
                self.frequencies_autocorr, high_beta_freq + high_beta_peak_width * 0.1 / 2)])


        return theta_power, alpha_power, low_beta_power, high_beta_power


    def isi_plot(self):
        fig, ax = plt.subplots(1, 1)
        # Ensure that spikes are from the target neuron, probably by requiring spike sorting before plotting.
        spike_times = self.main_times[np.array(list(self.features_spiketrain_indices), dtype=int)]
        isi = np.log(np.diff(spike_times))

        X = np.ravel(isi)
        X = np.ravel(isi).reshape(-1, 1)
        M_best = GaussianMixture(n_components=2, covariance_type='spherical').fit(X)
        x = np.linspace(np.log(1.5/1000), 0, 10000)
        y = np.exp(M_best.score_samples(x.reshape(-1, 1)))
        y_individual = M_best.predict_proba(x.reshape(-1, 1)) *  y[:, np.newaxis]

        ax.hist(isi, bins=20, density=True)
        ax.plot(x, y)
        ax.plot(x, y_individual, 'k--')

        ax.set_ylabel('Density')
        ax.set_xlabel('log(ISI)')
        ax.set_title('Log-Interspike-Interval Histogram')
        ax.legend(['Gaussian Mixture'])

        return fig


    def oscillations_plot(self):
        fig, ax = plt.subplots(4, 1, sharex=True)
        ax[0].plot(self.spike_oscillations_theta_times, self.spike_oscillations_theta_wave)
        ax[1].plot(self.spike_oscillations_alpha_times, self.spike_oscillations_alpha_wave)
        ax[2].plot(self.spike_oscillations_low_beta_times, self.spike_oscillations_low_beta_wave)
        ax[3].plot(self.spike_oscillations_high_beta_times, self.spike_oscillations_high_beta_wave)
        ax[3].set_xlabel('Time (s)')
        ax[0].set_ylabel('Theta')
        ax[1].set_ylabel('Alpha')
        ax[2].set_ylabel('Low Beta')
        ax[3].set_ylabel('High Beta')
        ax[0].set_title('Spiketrain Oscillations Waveforms')

        fig.tight_layout()

        return fig


    def lfp_plot(self):
        fig = plt.figure()
        gs = fig.add_gridspec(4, 2)
        ax1 = fig.add_subplot(gs[:, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[1, 1])
        ax4 = fig.add_subplot(gs[2, 1])
        ax5 = fig.add_subplot(gs[3, 1])

        # Fix the time arrays for the LFP waveforms.
        ax1.plot(self.lfp_psd_freqs, self.lfp_psd_power)

        ax2.plot(np.arange(0, self.main_times[-1], self.main_times[-1]/len(self.lfp_theta_wave)), self.lfp_theta_wave)
        ax3.plot(np.arange(0, self.main_times[-1], self.main_times[-1]/len(self.lfp_alpha_wave)), self.lfp_alpha_wave)
        ax4.plot(np.arange(0, self.main_times[-1], self.main_times[-1]/len(self.lfp_low_beta_wave)), self.lfp_low_beta_wave)
        ax5.plot(np.arange(0, self.main_times[-1], self.main_times[-1]/len(self.lfp_high_beta_wave)), self.lfp_high_beta_wave)

        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('Spectral Power')
        ax1.set_title('Power Spectral Density')
        ax2.set_ylabel('Theta')
        ax3.set_ylabel('Alpha')
        ax4.set_ylabel('Low Beta')
        ax5.set_ylabel('High Beta')
        ax5.set_xlabel('Time (s)')
        ax2.set_title('LFP Waveforms')

        fig.tight_layout()

        return fig