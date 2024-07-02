import numpy as np

class ChannelData():

    def __init__(self, id, raw_signal, fs):

        self.id = id
        self.raw_signal = raw_signal
        self.fs = fs

        self.filtered_signal = None
        self.threshold_factor = None
        self.threshold = None
        self.spike_indices = None
    

    def get_current_signal(self):
        # Returns active state of signal (filtered or raw) depending on whether filtering was executed.
        if self.filtered_signal is None:
            return self.raw_signal
        return self.filtered_signal
    

    def get_time_vector(self):
        return np.arange(len(self.get_current_signal())) / self.fs


class FileData():
    
    def __init__(self, filepath, channels, timebase=None):
        
        self.filepath = filepath
        self.channels = channels    # Dictionary of ChannelData objects. Keys are channel ID's.

        self.timebase = timebase    # If Spike2 file loaded.