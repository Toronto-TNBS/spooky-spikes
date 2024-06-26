from sonpy import lib as son
import math
from slots.backend import data


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
