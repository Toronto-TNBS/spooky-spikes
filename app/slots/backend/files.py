import math
from app.slots.backend import data
import scipy.io as spio
import mat73
import pandas as pd
import numpy as np
import neo


def load_spike2(filepath):
    load = neo.io.Spike2IO(filepath).read_block().segments[0].analogsignals[0]
    fs = float(load.sampling_rate.magnitude)
    channels = {}
    for i, chid in enumerate(load.array_annotations["channel_ids"]):
        signal = load[:, i].magnitude.flatten()
        channel_data = data.ChannelData(int(chid), signal, fs)
        channels[int(chid)] = channel_data
    file_data = data.FileData(filepath, channels)

    return file_data


def load_mat(filepath):
    try:
        file = spio.loadmat(filepath)
    except:
        file = mat73.loadmat(filepath)

    channels = {}
    for var in file.keys():
        if var[:6] == 'values':
            try:
                int(var[6:])    # If this raises an error, invalid entry follows "values" (e.g. values1.5, valuesfijgn)
            except:
                continue
            
            ch_id = int(var[6:])
            signal = np.ndarray.flatten(file[var])
            fs = float(np.ndarray.flatten(file['fs' + str(ch_id)])[0])

            channels[ch_id] = data.ChannelData(ch_id, signal, fs)
    
    file_data = data.FileData(filepath, channels)
    return file_data


def load_csv(filepath):
    file = pd.read_csv(filepath)
    channels = {}
    for var in file.columns:
        if var[:6] == 'values':
            try:
                int(var[6:])    
            except:
                continue

            ch_id = int(var[6:])
            signal = file[var].values
            fs = float(file['fs' + str(ch_id)][0])

            channels[ch_id] = data.ChannelData(ch_id, signal, fs)

    file_data = data.FileData(filepath, channels)
    return file_data


def save_features(filepath, channel, file):
    filtered = False
    if channel.filtered_signal is not None:
        filtered = True

    num_clusters = len(channel.get_existing_clusters())

    properties = pd.DataFrame(
        {
            'Channel': list(file.channels.keys())[channel.id],
            'Threshold': channel.threshold,
            'Threshold_Factor': channel.threshold_factor,
            'Filtering': filtered,
            'Highpass': channel.hpcutoff,
            'Lowpass': channel.lpcutoff,
            'Inverted': channel.signal_inverted,
            'Number_Clusters': num_clusters,
            'SNR': channel.features_qualitymetrics['SNR'],
            'Percent_ISI_Violations': channel.features_qualitymetrics['ISI'],
            'Firing_Rate': channel.features_patterned['FR'],
            'Burst_Index': channel.features_patterned['BI'],
            'CoV': channel.features_patterned['CV'],
            'Silhouette': channel.features_qualitymetrics['SS'],
            'Spiketrain_Theta_Power': channel.features_spiketrain['Power']['Theta'],
            'Spiketrain_Alpha_Power': channel.features_spiketrain['Power']['Alpha'],
            'Spiketrain_Low_Beta_Power': channel.features_spiketrain['Power']['LowBeta'],
            'Spiketrain_High_Beta_Power': channel.features_spiketrain['Power']['HighBeta'],
            'Theta_Mean_Burst_Duration': channel.features_spiketrain['BurstDur']['Theta'],
            'Alpha_Mean_Burst_Duration': channel.features_spiketrain['BurstDur']['Alpha'],
            'Low_Beta_Mean_Burst_Duration': channel.features_spiketrain['BurstDur']['LowBeta'],
            'High_Beta_Mean_Burst_Duration': channel.features_spiketrain['BurstDur']['HighBeta'],
            'LFP_Theta_Power': channel.features_lfp['Power']['Theta'],
            'LFP_Alpha_Power': channel.features_lfp['Power']['Alpha'],
            'LFP_Low_Beta_Power': channel.features_lfp['Power']['LowBeta'],
            'LFP_High_Beta_Power': channel.features_lfp['Power']['HighBeta']
        },
        index=[0]
    )
    try:
        properties.to_csv(filepath, index=None)
        return 0
    except:
        return -1


def save_segment(filepath, signal, events, fs, timebase=None):
    filetype = filepath.split('.')[-1]
    if filetype == 'mat':
        status = save_mat(filepath, signal, events, fs)
    return status


def save_mat(filepath, signal, events, fs):
    try:
        variables = {'values': signal, 'spike_indices': events, 'fs': fs}
        spio.savemat(filepath, variables)
        return 0
    except:
        return -1