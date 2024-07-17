import scipy.signal as sps
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


def butterworth_filter(signal, fs, hpcutoff=None, lpcutoff=None):
    if hpcutoff == None:
        b, a = sps.iirfilter(4, lpcutoff, btype='lowpass', fs=fs)
    elif lpcutoff == None:
        b, a = sps.iirfilter(4, hpcutoff, btype='highpass', fs=fs)
    else:
        b, a = sps.iirfilter(4, [hpcutoff, lpcutoff], fs=fs)
    
    return sps.filtfilt(b, a, signal)


def get_spike_threshold(signal, factor):
    return factor * np.median(np.abs(signal)) / 0.6745


def get_spike_indices(signal, threshold, fs, inverted):
    if inverted:
        signal = -signal.copy()    # Signal inverted for spike detection, factor and threshold remain positive.
    return sps.find_peaks(signal, height=threshold, distance=int(2/1000*fs))[0]


def get_spike_matrix(signal, spike_indices, window_size):
    spike_indices = spike_indices[(spike_indices > window_size) & (spike_indices < spike_indices[-1]-window_size)]
    spike_matrix = np.empty((1, window_size*2))
    for i in spike_indices:
        spike_matrix = np.append(spike_matrix, [signal[i-window_size: i+window_size]], axis=0)
    spike_matrix = np.delete(spike_matrix, 0, axis=0)
    return spike_matrix, spike_indices


def spike_sorting(spike_matrix, num_clusters):

    def clustering(spike_matrix, num_clusters):
        pca_results = PCA(n_components=2).fit_transform(spike_matrix)
        kmeans = KMeans(n_clusters=num_clusters, random_state=0).fit(pca_results)
        cluster_labels = kmeans.labels_
        return pca_results, cluster_labels

    if num_clusters is None:
        scores = []
        for c in range(2, 7):
            features_temp, label_temp = clustering(spike_matrix, c)
            scores.append(silhouette_score(features_temp, label_temp))
        num_clusters = np.array(scores).argmax() + 2

    features, labels = clustering(spike_matrix, num_clusters)
    score = silhouette_score(features, labels)

    return features, labels, score