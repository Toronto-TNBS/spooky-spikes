import numpy as np
from app.slots.backend import analysis

class ChannelData():

    def __init__(self, id, raw_signal, fs):

        self.id = id
        self.raw_signal = raw_signal
        self.fs = fs

        self.filtered_signal = None
        self.threshold_factor = None
        self.threshold = None
        self.spike_indices_all = None
        self.current_spike_indices = None
        self.signal_inverted = False
        self.spike_matrix = None

        self.spikesorting_clusters = {
            'r': None, 
            'b': None,
            'g': None,
            'c': None,
            'm': None,
            'y': None
        }
        self.selected_cluster_label = None

        # Features
        self.features_patterned = {'FR': None,
                                   'BI': None,
                                   'CV': None}
        self.features_qualitymetrics = {'SNR': None,
                                       'ISI': None,
                                       'SS': None}
        self.features_spiketrain = {'Power': {'Theta': None,
                                              'Alpha': None,
                                              'LowBeta': None,
                                              'HighBeta': None},
                                    'BurstDur': {'Theta': None,
                                              'Alpha': None,
                                              'LowBeta': None,
                                              'HighBeta': None}}
        self.features_lfp = {'Power': {'Theta': None,
                                       'Alpha': None,
                                       'LowBeta': None,
                                       'HighBeta': None},
                             'Wave': {'Theta': None,
                                      'Alpha': None,
                                      'LowBeta': None,
                                      'HighBeta': None},
                             'PSD': {'Freq': None,
                                     'Power': None}}
    

    def get_current_signal(self):
        # Returns active state of signal (filtered or raw) depending on whether filtering was executed.
        if self.filtered_signal is None:
            return self.raw_signal
        return self.filtered_signal
    

    def get_time_vector(self):
        return np.arange(len(self.get_current_signal())) / self.fs
    

    def get_spike_plot_colour(self):
        if self.selected_cluster_label is None:
            return 'r'
        return list(self.spikesorting_clusters.keys())[self.selected_cluster_label]
    

    def reset_spikesorting_clusters_data(self):
     for key in self.spikesorting_clusters.keys():
          self.spikesorting_clusters[key] = None
        

    def get_existing_clusters(self):
        return [(i, key) for i, key in enumerate(self.spikesorting_clusters.keys()) if self.spikesorting_clusters[key] is not None]
    

    def reset_channel(self):
        self.__init__(self.id, self.raw_signal, self.fs)
    

    def compute_features_patterned(self):
        self.features_patterned['FR'] = analysis.get_FR(self.current_spike_indices, self.fs)
        self.features_patterned['BI'] = analysis.get_BI(self.current_spike_indices)
        self.features_patterned['CV'] = analysis.get_CV(self.current_spike_indices)
    

    def compute_features_qualitymetrics(self, ss):
        # SS must be computed and assigned separately as it is required asynchronously with the rest of the metrics.
        if ss == None:
            ss = ''
        self.features_qualitymetrics['SNR'] = analysis.get_SNR(self.raw_signal, self.current_spike_indices, self.fs)
        self.features_qualitymetrics['ISI'] = analysis.get_pISIv(self.current_spike_indices, self.fs)
        self.features_qualitymetrics['SS'] = ss
    

    def compute_features_spiketrain(self):
        theta_power, alpha_power, lowbeta_power, highbeta_power = analysis.get_spiketrain_power(self.current_spike_indices, self.fs)
        theta_burst, alpha_burst, lowbeta_burst, highbeta_burst = analysis.get_spiketrain_burstduration(self.current_spike_indices, self.fs)
        
        self.features_spiketrain['Power']['Theta'] = theta_power
        self.features_spiketrain['Power']['Alpha'] = alpha_power
        self.features_spiketrain['Power']['LowBeta'] = lowbeta_power
        self.features_spiketrain['Power']['HighBeta'] = highbeta_power

        self.features_spiketrain['BurstDur']['Theta'] = theta_burst
        self.features_spiketrain['BurstDur']['Alpha'] = alpha_burst
        self.features_spiketrain['BurstDur']['LowBeta'] = lowbeta_burst
        self.features_spiketrain['BurstDur']['HighBeta'] = highbeta_burst
    

    def compute_features_lfp(self):
        powers, waves, psd = analysis.get_lfp_metrics(self.raw_signal, self.fs)

        self.features_lfp['Power']['Theta'] = powers[0]
        self.features_lfp['Power']['Alpha'] = powers[1]
        self.features_lfp['Power']['LowBeta'] = powers[2]
        self.features_lfp['Power']['HighBeta'] = powers[3]

        self.features_lfp['Wave']['Theta'] = waves[0]
        self.features_lfp['Wave']['Alpha'] = waves[1]
        self.features_lfp['Wave']['LowBeta'] = waves[2]
        self.features_lfp['Wave']['HighBeta'] = waves[3]

        self.features_lfp['PSD']['Freq'] = psd[0]
        self.features_lfp['PSD']['Power'] = psd[1]


    def reset_features(self):
        self.features_patterned = {'FR': None,
                                   'BI': None,
                                   'CV': None}
        self.features_qualitymetrics = {'SNR': None,
                                       'ISI': None,
                                       'SS': None}
        self.features_spiketrain = {'Power': {'Theta': None,
                                              'Alpha': None,
                                              'LowBeta': None,
                                              'HighBeta': None},
                                    'BurstDur': {'Theta': None,
                                              'Alpha': None,
                                              'LowBeta': None,
                                              'HighBeta': None}}
        self.features_lfp = {'Power': {'Theta': None,
                                       'Alpha': None,
                                       'LowBeta': None,
                                       'HighBeta': None},
                             'Wave': {'Theta': None,
                                      'Alpha': None,
                                      'LowBeta': None,
                                      'HighBeta': None},
                             'PSD': {'Freq': None,
                                     'Power': None}}


class FileData():
    
    def __init__(self, filepath, channels, timebase=None):
        
        self.filepath = filepath
        self.channels = channels    # Dictionary of ChannelData objects. Keys are channel ID's.

        self.timebase = timebase    # If Spike2 file loaded.