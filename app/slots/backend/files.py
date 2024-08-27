from sonpy import lib as son
import math
from app.slots.backend import data
import scipy.io as spio
import mat73
import pandas as pd
import numpy as np


def load_spike2(filepath):
    file = son.SonFile(filepath, True)

    wave_channels = [i for i in range(file.MaxChannels()) if file.ChannelType(i) == son.DataType.Adc or file.ChannelType == son.DataType.RealWave]
    channels = {}
    for channel_idx in wave_channels:
        timebase = file.GetTimeBase()    # s/tick
        fs = 1 / (file.ChannelDivide(channel_idx) * timebase)    # 1 / ([tick] [s/tick])
        num_points = math.floor(file.ChannelMaxTime(channel_idx) * timebase * fs)    # [tick] [s/tick] [sample/s]

        signal = file.ReadFloats(channel_idx, num_points, 0)

        channel_data = data.ChannelData(channel_idx, signal, fs)
        channels[channel_idx+1] = channel_data
    
    file_data = data.FileData(filepath, channels, timebase=timebase)

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
    if filetype == 'smr':
        status = save_smr(filepath, signal, events, fs, timebase)
    elif filetype == 'mat':
        status = save_mat(filepath, signal, events, fs)
    return status


def save_smr(filepath, signal, events, fs, timebase):
    events = events/fs
    wave_time_base = timebase
    wave_sample_ticks = int(1 / (fs * wave_time_base))

    filename = filepath
    new = son.SonFile(filename)

    wave_offset = 0
    wave_scale = 1
    wave_Y_low = -5
    wave_Y_high = 5
    wave_units = 'V'
    wave_title = 'Waveform'
    wave_Fs = fs
    time_start = 0

    new.SetTimeBase(wave_time_base)
    new_wave_channel = 0

    new.SetWaveChannel(new_wave_channel, wave_sample_ticks, son.DataType.Adc, wave_Fs)
    new.SetChannelTitle(new_wave_channel, wave_title)
    new.SetChannelUnits(new_wave_channel, wave_units)
    new.SetChannelScale(new_wave_channel, wave_scale)
    new.SetChannelOffset(new_wave_channel, wave_offset)
    new.SetChannelYRange(new_wave_channel, wave_Y_low, wave_Y_high)

    spike_channel = new.GetFreeChannel()
    max_event_rate = 1 / (wave_time_base * wave_sample_ticks)

    new.SetEventChannel(spike_channel, max_event_rate, son.DataType.EventFall)
    new.SetChannelTitle(spike_channel, 'ST-1')
    new.SetChannelUnits(spike_channel, wave_units)
    new.SetChannelScale(spike_channel, wave_scale)
    new.SetChannelOffset(spike_channel, wave_offset)
    new.SetChannelYRange(spike_channel, wave_Y_low, wave_Y_high)
    del new

    reload = son.SonFile(filename, False)
    scaled_wave = np.array(6553.6 * signal, dtype=int)
    write_wave = reload.WriteInts(new_wave_channel, scaled_wave, time_start)
    if write_wave < 0:
        return -1

    spike_ticks = [int(i / wave_time_base) for i in events]
    spike_times_formatted = np.array(spike_ticks, dtype=np.int64)
    write_train = reload.WriteEvents(spike_channel, spike_times_formatted)
    if write_train < 0:
        return -1

    del reload

    return 0


def save_mat(filepath, signal, events, fs):
    try:
        variables = {'values': signal, 'spike_indices': events, 'fs': fs}
        spio.savemat(filepath, variables)
        return 0
    except:
        return -1