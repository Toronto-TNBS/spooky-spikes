import scipy.signal as sps
import numpy as np


def butterworth_filter(signal, fs, hpcutoff=None, lpcutoff=None):
    if hpcutoff == None:
        b, a = sps.iirfilter(4, lpcutoff, btype='lowpass', fs=fs)
    elif lpcutoff == None:
        b, a = sps.iirfilter(4, hpcutoff, btype='highpass', fs=fs)
    else:
        b, a = sps.iirfilter(4, [hpcutoff, lpcutoff], fs=fs)
    
    return sps.filtfilt(b, a, signal)


def get_spike_threshold(signal, factor):
    # Account for inverted segments.
    return factor * np.median(np.abs(signal)) / 0.6745


def get_spike_indices(signal, threshold, fs):
    # Account later for inverted segments.
    return sps.find_peaks(-signal, height=threshold, distance=int(2/1000*fs))[0]