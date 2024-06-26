
class ChannelData():

    def __init__(self, id, raw_signal, fs):

        self.id = id
        self.raw_signal = raw_signal
        self.fs = fs


class FileData():
    
    def __init__(self, filepath, channels, timebase=None):
        
        self.filepath = filepath
        self.channels = channels    # Dictionary of ChannelData objects. Keys are channel ID's.

        self.timebase = timebase    # If Spike2 file loaded.