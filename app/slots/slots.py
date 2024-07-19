from PySide6.QtWidgets import *
from PySide6 import QtCore
from PySide6 import QtGui
from slots.backend import data
from slots.backend import files
from slots.backend import analysis
from slots.backend import miscellaneous as misc
import numpy as np
import pyqtgraph as pg


def button_file_clicked(app):
    filetypes = 'Spike2 Datafiles (*.smr *.smrx);;MAT-Files (*.mat);;Comma Separated Values (*.csv)'
    filepath = QFileDialog.getOpenFileName(caption='Open File', dir='.', filter=filetypes)[0]
    if filepath == '':
        return

    file_data = files.load_spike2(filepath)
    app.filedata = file_data

    app.label_file.setText(file_data.filepath.split('/')[-1])

    app.dropdown_channel.addItems([f'Ch {ch_id}' for ch_id in app.filedata.channels.keys()])    


def dropdown_channel_changed(app):
    ch_id = int(app.dropdown_channel.currentText().split(' ')[-1])

    if app.channeldata != None:
        app.channeldata.filtered_signal = None    # Reset filter settings of previously selected channel if exists.

    app.channeldata = app.filedata.channels[ch_id]
    
    app.plot2_tab_main.addItem(pg.PlotDataItem(np.arange(len(app.channeldata.raw_signal)) / app.channeldata.fs, app.channeldata.raw_signal, pen=pg.mkPen('cornflowerblue', width=1.5)))
    
    app.check_filtering.setEnabled(True)
    app.entry_tab_main_madfactor.setEnabled(True)
    app.button_tab_main_threshold.setEnabled(True)
    app.button_tab_features_lfpplot.setEnabled(True)
    app.button_tab_features_autocorrelationplot.setEnabled(True)

    # Need to address what happens to data objects when different channel is selected.


def check_filtering_changed(app):
    check_state = app.check_filtering.isChecked()

    app.entry_tab_main_hpcutoff.setEnabled(check_state)
    app.entry_tab_main_lpcutoff.setEnabled(check_state)
    app.button_tab_main_filtering.setEnabled(check_state)

    if not check_state:
        app.channeldata.filtered_signal = None
        current_signal_item = app.plot2_tab_main.listDataItems()[0]
        app.plot2_tab_main.removeItem(current_signal_item)
        app.plot2_tab_main.addItem(pg.PlotDataItem(np.arange(len(app.channeldata.raw_signal)) / app.channeldata.fs, app.channeldata.raw_signal, pen=pg.mkPen('cornflowerblue', width=1.5)))

        app.entry_tab_main_hpcutoff.clear()
        app.entry_tab_main_lpcutoff.clear()
        app.label_tab_main_filterstatus.setText('HP: None\nLP: None')

        # Re-run spike detection after filtering disabled.
        app.check_tab_main_spiketrain.setChecked(False)
        app.check_tab_main_eventtimes.setChecked(False)
        app.check_tab_main_thresholdbar.setChecked(False)
        button_tab_main_threshold_clicked(app)


def check_spikesorting_changed(app):
    check_state = app.check_spikesorting.isChecked()
    cluster_colours = list(app.channeldata.spikesorting_clusters.keys())
    
    misc.remove_plot_items(plotitem=app.plot_tab_spikesorting, itemtype=pg.ScatterPlotItem)    # Reset plot every time triggered.
    misc.reset_spikesorting_clusters_data(app)
    if not check_state:
        app.channeldata.silhouette_score = None
        app.label_tab_features_qualitymetrics_silhouette_status.setText('')
        app.dropdown_tab_main_cluster.clear()
        app.channeldata.spike_matrix = None    # Reset spike matrix for new spike sorting operation.

        # Reset to all spikes.
        app.channeldata.current_spike_indices = app.channeldata.spike_indices_all
        app.channeldata.selected_cluster_label = None
        app.check_tab_main_spiketrain.setChecked(False)
        app.check_tab_main_eventtimes.setChecked(False)
        app.check_tab_main_spiketrain.setChecked(True)
        app.check_tab_main_eventtimes.setChecked(True)
        return
    
    num_clusters = misc.check_convert_string(app.dropdown_tab_spikesorting_clusters.currentText(), int)
    
    current_signal = app.channeldata.get_current_signal()
    if app.channeldata.spike_matrix is None:    # Keep spike matrix in memory for efficiency.
        app.channeldata.spike_matrix, app.channeldata.spike_indices_all = analysis.get_spike_matrix(current_signal, app.channeldata.spike_indices_all, int(1/1000*app.channeldata.fs))
    pca_results, cluster_labels, silhouette_score = analysis.spike_sorting(app.channeldata.spike_matrix, num_clusters)

    for c in np.unique(cluster_labels):
        cluster_indices = np.where(cluster_labels == c)[0]
        cluster_plot_item = pg.ScatterPlotItem(x=pca_results[cluster_indices, 0], y=pca_results[cluster_indices, 1], pen=pg.mkPen(cluster_colours[c], width=1.25), brush=pg.mkBrush(cluster_colours[c]))
        app.plot_tab_spikesorting.addItem(cluster_plot_item)
        
        app.channeldata.spikesorting_clusters[cluster_colours[c]] = app.channeldata.spike_indices_all[cluster_indices]
    
    app.channeldata.silhouette_score = silhouette_score
    app.label_tab_features_qualitymetrics_silhouette_status.setText(str(round(app.channeldata.silhouette_score, 2)))

    existing_clusters = misc.get_existing_clusters(app)
    app.channeldata.selected_cluster_label = np.argmax([np.mean(app.channeldata.spikesorting_clusters[c]) for i, c in existing_clusters])

    cluster_colour_names = ['red', 'blue', 'green', 'cyan', 'magenta', 'yellow']
    app.dropdown_tab_main_cluster.clear()
    app.dropdown_tab_main_cluster.addItems(cluster_colour_names[:len(existing_clusters)] + ['none'])
    app.dropdown_tab_main_cluster.setCurrentIndex(app.channeldata.selected_cluster_label)


def check_invertthreshold_changed(app):
    app.entry_tab_main_madfactor.setText(str(-app.channeldata.threshold_factor))
    app.channeldata.signal_inverted = app.check_invertthreshold.isChecked()
    app.check_tab_main_spiketrain.setChecked(False)
    app.check_tab_main_eventtimes.setChecked(False)
    app.check_tab_main_thresholdbar.setChecked(False)

    app.check_spikesorting.setChecked(False)    # Turn off spike sorting.
    button_tab_main_threshold_clicked(app)    # Restarts the threshold analysis with inverted parameter.


def button_tab_main_filtering_clicked(app):

    hpcutoff = misc.check_convert_string(app.entry_tab_main_hpcutoff.text(), float)
    lpcutoff = misc.check_convert_string(app.entry_tab_main_lpcutoff.text(), float)
    
    app.label_tab_main_filterstatus.setText(f'HP: {hpcutoff}\nLP: {lpcutoff}')
    
    if hpcutoff == None and lpcutoff == None:
        app.channeldata.filtered_signal = None
        current_signal_item = app.plot2_tab_main.listDataItems()[0]
        app.plot2_tab_main.removeItem(current_signal_item)
        app.plot2_tab_main.addItem(pg.PlotDataItem(np.arange(len(app.channeldata.raw_signal)) / app.channeldata.fs, app.channeldata.raw_signal, pen=pg.mkPen('cornflowerblue', width=1.5)))
        return
    
    app.channeldata.filtered_signal = analysis.butterworth_filter(app.channeldata.raw_signal, app.channeldata.fs, hpcutoff, lpcutoff)
    
    current_signal_item = app.plot2_tab_main.listDataItems()[0]    # Will have to assess list ordering issues when I add threshold bar to this plot and removing objects.
    app.plot2_tab_main.removeItem(current_signal_item)
    app.plot2_tab_main.addItem(pg.PlotDataItem(np.arange(len(app.channeldata.filtered_signal)) / app.channeldata.fs, app.channeldata.filtered_signal, pen=pg.mkPen('cornflowerblue', width=1.5)))


def button_tab_main_threshold_clicked(app):

    factor = misc.check_convert_string(app.entry_tab_main_madfactor.text(), float)
    
    if factor == None:
        # Reset threshold analysis if nonsense value inserted.
        app.channeldata.threshold = None
        app.channeldata.threshold_factor = None
        app.channeldata.spike_indices_all = None
        app.check_tab_main_spiketrain.setChecked(False)
        app.check_tab_main_eventtimes.setChecked(False)
        app.check_tab_main_thresholdbar.setChecked(False)
        app.label_tab_main_thresholdstatus.setText(f'Threshold: None (Factor: None)')
        app.check_spikesorting.setEnabled(False)
        return
    
    current_signal = app.channeldata.get_current_signal()
    factor = abs(factor)    # Factor and threshold in backend always positive. Only signal is inverted in backend.

    threshold = analysis.get_spike_threshold(current_signal, factor)
    
    app.channeldata.threshold_factor = factor
    app.channeldata.threshold = threshold
    app.label_tab_main_thresholdstatus.setText(f'Threshold: {round(app.channeldata.threshold, 3)} (Factor: {app.channeldata.threshold_factor})')
    
    spike_indices = analysis.get_spike_indices(current_signal, app.channeldata.threshold, app.channeldata.fs, app.channeldata.signal_inverted)
    app.channeldata.spike_indices_all = spike_indices
    app.channeldata.current_spike_indices = app.channeldata.spike_indices_all

    misc.remove_plot_items(plotitem=app.plot2_tab_main, itemtype=pg.ScatterPlotItem)    # Remove previously plotted spikes if they exist.

    app.check_invertthreshold.setEnabled(True)
    app.check_tab_main_spiketrain.setEnabled(True)
    app.check_tab_main_eventtimes.setEnabled(True)
    app.check_tab_main_thresholdbar.setEnabled(True)
    app.check_spikesorting.setEnabled(True)
    # Following resets properly display spike-related contents on main plot at any point of UI event history.
    # In order for plot display items to be refreshed, checks must be set to False then True again.
    # Can now use this threshold-analysis-related function to perform analysis resets.
    app.check_tab_main_spiketrain.setChecked(False)
    app.check_tab_main_eventtimes.setChecked(False)
    app.check_tab_main_thresholdbar.setChecked(False)
    app.check_tab_main_spiketrain.setChecked(True)
    app.check_tab_main_eventtimes.setChecked(True)
    app.check_tab_main_thresholdbar.setChecked(True)


def check_tab_main_spiketrain_changed(app):
    check_state = app.check_tab_main_spiketrain.isChecked()
    if check_state:
        colour = app.channeldata.get_spike_plot_colour()
        peaks_scatter_item = pg.ScatterPlotItem(x=app.channeldata.get_time_vector()[app.channeldata.current_spike_indices], y=app.channeldata.get_current_signal()[app.channeldata.current_spike_indices], pen=pg.mkPen(colour), brush=pg.mkBrush(colour), size=5)
        app.plot2_tab_main.addItem(peaks_scatter_item)
        return
    misc.remove_plot_items(app.plot2_tab_main, pg.ScatterPlotItem)


def check_tab_main_eventtimes_changed(app):
    check_state = app.check_tab_main_eventtimes.isChecked()
    if check_state:
        colour = app.channeldata.get_spike_plot_colour()
        events_scatter_item = pg.ScatterPlotItem(x=app.channeldata.get_time_vector()[app.channeldata.current_spike_indices], y=[0]*len(app.channeldata.current_spike_indices), symbol=misc.vertical_line_symbol(), pen=pg.mkPen(colour), brush=pg.mkBrush(colour))
        app.plot1_tab_main.addItem(events_scatter_item)
        return
    misc.remove_plot_items(app.plot1_tab_main, pg.ScatterPlotItem)


def check_tab_main_thresholdbar_changed(app):
    check_state = app.check_tab_main_thresholdbar.isChecked()
    if check_state:
        threshold = app.channeldata.threshold
        if app.channeldata.signal_inverted:
            threshold *= -1
        app.threshold_bar_item = pg.InfiniteLine(pos=threshold, angle=0, pen=pg.mkPen('k', width=1, style=QtCore.Qt.DashLine))
        app.plot2_tab_main.addItem(app.threshold_bar_item)
        return
    app.plot2_tab_main.removeItem(app.threshold_bar_item)


def dropdown_tab_main_cluster_changed(app):
    if app.dropdown_tab_main_cluster.currentText() == 'none' or app.dropdown_tab_main_cluster.currentText() == '':
        # Empty string for if dropdown items are cleared.
        app.channeldata.current_spike_indices = app.channeldata.spike_indices_all
        app.channeldata.selected_cluster_label = None
    else:
        app.channeldata.selected_cluster_label = app.dropdown_tab_main_cluster.currentIndex()
        selected_key = list(app.channeldata.spikesorting_clusters.keys())[app.channeldata.selected_cluster_label]
        app.channeldata.current_spike_indices = app.channeldata.spikesorting_clusters[selected_key]

    app.check_tab_main_spiketrain.setChecked(False)
    app.check_tab_main_eventtimes.setChecked(False)
    app.check_tab_main_spiketrain.setChecked(True)
    app.check_tab_main_eventtimes.setChecked(True)


def dropdown_tab_spikesorting_clusters_changed(app):
    check_spikesorting_changed(app)