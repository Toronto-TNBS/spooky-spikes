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