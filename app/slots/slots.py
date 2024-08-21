from PySide6.QtWidgets import *
from PySide6 import QtCore
from PySide6 import QtGui
from app.slots.backend import data
from app.slots.backend import files
from app.slots.backend import analysis
from app.slots.backend import miscellaneous as misc
import numpy as np
import pyqtgraph as pg
from sklearn.mixture import GaussianMixture
import matplotlib.pyplot as plt
from scipy.signal import hilbert


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

    if app.channeldata != None:    # If channel selected previously.
        # Reset analysis:
        app.check_spikesorting.setChecked(False)
        app.check_filtering.setChecked(False)
        app.check_invertthreshold.setChecked(False)
        
        app.entry_tab_main_madfactor.setText('')    # Set up for clearing operation by subsequent function.
        button_tab_main_threshold_clicked(app)    # After setup: resets plot disp checks, spikesorting, threshold, clears unit features.
        
        misc.remove_plot_items(app.plot2_tab_main, pg.PlotDataItem)
        app.channeldata.reset_channel()    # Reset full analysis for currently-selected channel. Run at end since above functions need existing values.


    app.channeldata = app.filedata.channels[ch_id]
    
    app.plot2_tab_main.addItem(pg.PlotDataItem(app.channeldata.get_time_vector(), app.channeldata.raw_signal, pen=pg.mkPen('cornflowerblue', width=1.5)))
    
    app.check_filtering.setEnabled(True)
    app.entry_tab_main_madfactor.setEnabled(True)
    app.button_tab_main_threshold.setEnabled(True)
    app.button_tab_features_lfpplot.setEnabled(True)

    # Compute and display features
    app.channeldata.compute_features_lfp()
    misc.update_lfp_features_display(app)


def check_filtering_changed(app):
    check_state = app.check_filtering.isChecked()

    app.entry_tab_main_hpcutoff.setEnabled(check_state)
    app.entry_tab_main_lpcutoff.setEnabled(check_state)
    app.button_tab_main_filtering.setEnabled(check_state)

    if not check_state:
        app.check_spikesorting.setChecked(False)

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
    app.channeldata.reset_spikesorting_clusters_data()
    if not check_state:
        app.channeldata.features_qualitymetrics['SS'] = None
        app.dropdown_tab_main_cluster.clear()    # This line triggers the dropdown function. Initiates computation with all spikes.
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
        app.channeldata.spike_matrix = analysis.get_spike_matrix(current_signal, app.channeldata.spike_indices_all, int(1/1000*app.channeldata.fs))
    pca_results, cluster_labels, silhouette_score = analysis.spike_sorting(app.channeldata.spike_matrix, num_clusters)

    for c in np.unique(cluster_labels):
        cluster_indices = np.where(cluster_labels == c)[0]
        cluster_plot_item = pg.ScatterPlotItem(x=pca_results[cluster_indices, 0], y=pca_results[cluster_indices, 1], pen=pg.mkPen(cluster_colours[c], width=1.25), brush=pg.mkBrush(cluster_colours[c]))
        app.plot_tab_spikesorting.addItem(cluster_plot_item)
        
        app.channeldata.spikesorting_clusters[cluster_colours[c]] = app.channeldata.spike_indices_all[cluster_indices]
    
    existing_clusters = app.channeldata.get_existing_clusters()
    app.channeldata.selected_cluster_label = np.argmax([np.mean(app.channeldata.spikesorting_clusters[c]) for i, c in existing_clusters])

    app.channeldata.features_qualitymetrics['SS'] = silhouette_score    # Need to set this before setting cluster index in dropdown, which performs computation.

    cluster_colour_names = ['red', 'blue', 'green', 'cyan', 'magenta', 'yellow']
    app.dropdown_tab_main_cluster.clear()
    app.dropdown_tab_main_cluster.addItems(cluster_colour_names[:len(existing_clusters)] + ['none'])
    app.dropdown_tab_main_cluster.setCurrentIndex(app.channeldata.selected_cluster_label)    # Setting index computes all features within corresponding slot.


def check_invertthreshold_changed(app):
    app.check_spikesorting.setChecked(False)
    app.channeldata.signal_inverted = app.check_invertthreshold.isChecked()
    button_tab_main_threshold_clicked(app)    # Restarts the threshold analysis with inverted parameter. Displays results on main plot.


def button_tab_main_filtering_clicked(app):
    app.check_spikesorting.setChecked(False)

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

    button_tab_main_threshold_clicked(app)    # Runs after (re)filtering according to user-defined parameters.


def button_tab_main_threshold_clicked(app):

    app.check_spikesorting.setChecked(False)    # Does not run function if already OFF.

    factor = misc.check_convert_string(app.entry_tab_main_madfactor.text(), float)
    
    if factor == None:
        # Reset threshold analysis if nonsense value inserted.
        app.check_spikesorting.setChecked(False)
        app.check_spikesorting.setEnabled(False)

        app.button_tab_features_isiplot.setEnabled(False)
        app.button_tab_features_spikeoscillationsplot.setEnabled(False)
        app.button_tab_features_autocorrelationplot.setEnabled(False)

        app.channeldata.threshold = None
        app.channeldata.threshold_factor = None
        app.channeldata.spike_indices_all = None
        app.check_tab_main_spiketrain.setChecked(False)
        app.check_tab_main_eventtimes.setChecked(False)
        app.check_tab_main_thresholdbar.setChecked(False)
        app.label_tab_main_thresholdstatus.setText(f'Threshold: None (Factor: None)')
        misc.update_unit_features_display(app, clear=True)    # Clear unit features.
        return
    
    current_signal = app.channeldata.get_current_signal()
    factor = abs(factor)    # Factor and threshold in backend always positive. Only signal is inverted in backend.

    threshold = analysis.get_spike_threshold(current_signal, factor)
    
    app.channeldata.threshold_factor = factor
    app.channeldata.threshold = threshold

    disp_threshold = app.channeldata.threshold
    # Need to display a negative threshold when inverted. Keep scaler positive. Checkbox indicates inversion.
    if app.channeldata.signal_inverted:
        disp_threshold *= -1
    app.label_tab_main_thresholdstatus.setText(f'Threshold: {round(disp_threshold, 3)} (Factor: {app.channeldata.threshold_factor})')
    
    spike_indices = analysis.get_spike_indices(current_signal, app.channeldata.threshold, app.channeldata.fs, app.channeldata.signal_inverted)
    app.channeldata.spike_indices_all = spike_indices
    app.channeldata.current_spike_indices = app.channeldata.spike_indices_all

    misc.remove_plot_items(plotitem=app.plot2_tab_main, itemtype=pg.ScatterPlotItem)    # Remove previously plotted spikes if they exist.

    app.check_invertthreshold.setEnabled(True)
    app.check_tab_main_spiketrain.setEnabled(True)
    app.check_tab_main_eventtimes.setEnabled(True)
    app.check_tab_main_thresholdbar.setEnabled(True)
    app.check_spikesorting.setEnabled(True)

    app.button_tab_features_isiplot.setEnabled(True)
    app.button_tab_features_spikeoscillationsplot.setEnabled(True)
    app.button_tab_features_autocorrelationplot.setEnabled(True)

    # Following resets properly display spike-related contents on main plot at any point of UI event history.
    # In order for plot display items to be refreshed, checks must be set to False then True again.
    # Can now use this threshold-analysis-related function to perform analysis resets.
    app.check_tab_main_spiketrain.setChecked(False)
    app.check_tab_main_eventtimes.setChecked(False)
    app.check_tab_main_thresholdbar.setChecked(False)
    app.check_tab_main_spiketrain.setChecked(True)
    app.check_tab_main_eventtimes.setChecked(True)
    app.check_tab_main_thresholdbar.setChecked(True)

    # Compute and display features
    app.channeldata.compute_features_patterned()
    app.channeldata.compute_features_qualitymetrics(ss=None)    # Threshold operation always run when spike sorting is off.
    app.channeldata.compute_features_spiketrain()
    misc.update_unit_features_display(app)


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

    # Compute and display features for newly selected cluster.
    # SS already computed from previous spike sorting operation. Re-use.
    app.channeldata.compute_features_patterned()
    app.channeldata.compute_features_qualitymetrics(ss=app.channeldata.features_qualitymetrics['SS'])
    app.channeldata.compute_features_spiketrain()
    misc.update_unit_features_display(app, ss=app.channeldata.features_qualitymetrics['SS'])


def dropdown_tab_spikesorting_clusters_changed(app):
    check_spikesorting_changed(app)


def button_tab_features_isiplot_clicked(app):    
    
    spike_times = app.channeldata.get_time_vector()[np.array(list(app.channeldata.current_spike_indices), dtype=int)]
    isi = np.log(np.diff(spike_times))

    X = np.ravel(isi)
    X = np.ravel(isi).reshape(-1, 1)
    M_best = GaussianMixture(n_components=2, covariance_type='spherical').fit(X)
    x = np.linspace(np.log(1.5/1000), 0, 10000)
    y = np.exp(M_best.score_samples(x.reshape(-1, 1)))
    y_individual = M_best.predict_proba(x.reshape(-1, 1)) *  y[:, np.newaxis]

    densities, edges, _ = plt.hist(isi, bins=20, density=True)

    layout = pg.GraphicsLayoutWidget()
    plot = layout.addPlot(row=0, col=0)
    plot.addLegend(offset=(-1, 1))

    plot.addItem(pg.BarGraphItem(x0=edges[:-1], width=np.diff(edges), height=densities, brush=pg.mkBrush('#e24a33')))
    plot.addItem(pg.PlotDataItem(x, y, name='Gaussian Mixture', pen=pg.mkPen('cornflowerblue', width=2)))
    plot.addItem(pg.PlotDataItem(x, y_individual[:, 0], pen=pg.mkPen('k', style=QtCore.Qt.DashLine, width=2)))
    plot.addItem(pg.PlotDataItem(x, y_individual[:, 1], pen=pg.mkPen('k', style=QtCore.Qt.DashLine, width=2)))

    plot.getAxis('left').setLabel('Density')
    plot.getAxis('bottom').setLabel('log(ISI)')
    plot.setTitle('Log-Interspike-Interval Histogram')
    layout.setBackground('white')

    app.generate_plot_window(layout, 'isi')


def button_tab_features_spikeoscillationsplot_clicked(app):

    waves, times, burst_threshold = analysis.get_spiketrain_burstduration(app.channeldata.current_spike_indices, app.channeldata.fs, waves=True)
    spike_oscillations_theta_wave = waves[0]
    spike_oscillations_alpha_wave = waves[1]
    spike_oscillations_low_beta_wave = waves[2]
    spike_oscillations_high_beta_wave = waves[3]
    spike_oscillations_theta_times = times[0]
    spike_oscillations_alpha_times = times[1]
    spike_oscillations_low_beta_times = times[2]
    spike_oscillations_high_beta_times = times[3]

    burst_theta = np.abs(hilbert(spike_oscillations_theta_wave))
    burst_alpha = np.abs(hilbert(spike_oscillations_alpha_wave))
    burst_low_beta = np.abs(hilbert(spike_oscillations_low_beta_wave))
    burst_high_beta = np.abs(hilbert(spike_oscillations_high_beta_wave))

    layout = pg.GraphicsLayoutWidget()
    plot0 = layout.addPlot(row=0, col=0)
    plot1 = layout.addPlot(row=1, col=0)
    plot2 = layout.addPlot(row=2, col=0)
    plot3 = layout.addPlot(row=3, col=0)
    # Need to make it such that they share common axes. Largest range gets linked.
    wave_colour = '#e24a33'

    plot0.addItem(pg.PlotDataItem(spike_oscillations_theta_times, spike_oscillations_theta_wave, pen=pg.mkPen(wave_colour, width=1.5)))
    plot1.addItem(pg.PlotDataItem(spike_oscillations_alpha_times, spike_oscillations_alpha_wave, pen=pg.mkPen(wave_colour, width=1.5)))
    plot2.addItem(pg.PlotDataItem(spike_oscillations_low_beta_times, spike_oscillations_low_beta_wave, pen=pg.mkPen(wave_colour, width=1.5)))
    plot3.addItem(pg.PlotDataItem(spike_oscillations_high_beta_times, spike_oscillations_high_beta_wave, pen=pg.mkPen(wave_colour, width=1.5)))

    plot0.addItem(pg.PlotDataItem(spike_oscillations_theta_times, burst_theta, pen=pg.mkPen('cornflowerblue', width=1.5)))
    plot1.addItem(pg.PlotDataItem(spike_oscillations_alpha_times, burst_alpha, pen=pg.mkPen('cornflowerblue', width=1.5)))
    plot2.addItem(pg.PlotDataItem(spike_oscillations_low_beta_times, burst_low_beta, pen=pg.mkPen('cornflowerblue', width=1.5)))
    plot3.addItem(pg.PlotDataItem(spike_oscillations_high_beta_times, burst_high_beta, pen=pg.mkPen('cornflowerblue', width=1.5)))

    plot0.addItem(pg.InfiniteLine(pos=burst_threshold, angle=0, pen=pg.mkPen('k', width=1.5)))
    plot1.addItem(pg.InfiniteLine(pos=burst_threshold, angle=0, pen=pg.mkPen('k', width=1.5)))
    plot2.addItem(pg.InfiniteLine(pos=burst_threshold, angle=0, pen=pg.mkPen('k', width=1.5)))
    plot3.addItem(pg.InfiniteLine(pos=burst_threshold, angle=0, pen=pg.mkPen('k', width=1.5)))

    def boundary_plotitems(times, wave, start, duration):
        y_range = [np.min(wave), np.max(wave)]
        x0 = times[start]
        x1 = times[start] + duration
        return pg.PlotDataItem([x0, x0], y_range), pg.PlotDataItem([x1, x1], y_range)

    above_threshold_theta = burst_theta > burst_threshold
    cross_indices_theta = np.where(np.diff(above_threshold_theta))[0]
    if len(cross_indices_theta) > 1:
        if burst_theta[cross_indices_theta[0] + 1] <= burst_threshold:
            # If so, remove the first index from cross_indices
            cross_indices_theta = cross_indices_theta[1:]
        for start, end in zip(cross_indices_theta[:-1:2], cross_indices_theta[1::2]):
            duration = spike_oscillations_theta_times[end] - spike_oscillations_theta_times[
                start]  # Calculate the duration of the burst
            if duration > 0.1:
                x0_item, x1_item = boundary_plotitems(spike_oscillations_theta_times, spike_oscillations_theta_wave, start, duration)
                plot0.addItem(pg.FillBetweenItem(x0_item, x1_item, brush=pg.mkBrush('lightsteelblue')))

    above_threshold_alpha = burst_alpha > burst_threshold
    cross_indices_alpha = np.where(np.diff(above_threshold_alpha))[0]
    if len(cross_indices_alpha) > 1:
        if burst_alpha[cross_indices_alpha[0] + 1] <= burst_threshold:
            # If so, remove the first index from cross_indices
            cross_indices_alpha = cross_indices_alpha[1:]
        for start, end in zip(cross_indices_alpha[:-1:2], cross_indices_alpha[1::2]):
            duration = spike_oscillations_alpha_times[end] - spike_oscillations_alpha_times[
                start]  # Calculate the duration of the burst
            if duration > 0.1:
                x0_item, x1_item = boundary_plotitems(spike_oscillations_alpha_times, spike_oscillations_alpha_wave, start, duration)
                plot1.addItem(pg.FillBetweenItem(x0_item, x1_item, brush=pg.mkBrush('lightsteelblue')))

    above_threshold_low_beta = burst_low_beta > burst_threshold
    cross_indices_low_beta = np.where(np.diff(above_threshold_low_beta))[0]
    if len(cross_indices_low_beta) > 1:
        if burst_low_beta[cross_indices_low_beta[0] + 1] <= burst_threshold:
            # If so, remove the first index from cross_indices
            cross_indices_low_beta = cross_indices_low_beta[1:]
        for start, end in zip(cross_indices_low_beta[:-1:2], cross_indices_low_beta[1::2]):
            duration = spike_oscillations_low_beta_times[end] - spike_oscillations_low_beta_times[
                start]  # Calculate the duration of the burst
            if duration > 0.1:
                x0_item, x1_item = boundary_plotitems(spike_oscillations_low_beta_times, spike_oscillations_low_beta_wave, start, duration)
                plot2.addItem(pg.FillBetweenItem(x0_item, x1_item, brush=pg.mkBrush('lightsteelblue')))

    above_threshold_high_beta = burst_high_beta > burst_threshold
    cross_indices_high_beta = np.where(np.diff(above_threshold_high_beta))[0]
    if len(cross_indices_high_beta) > 1:
        if burst_high_beta[cross_indices_high_beta[0] + 1] <= burst_threshold:
            # If so, remove the first index from cross_indices
            cross_indices_high_beta = cross_indices_high_beta[1:]
        for start, end in zip(cross_indices_high_beta[:-1:2], cross_indices_high_beta[1::2]):
            duration = spike_oscillations_high_beta_times[end] - spike_oscillations_high_beta_times[
                start]  # Calculate the duration of the burst
            if duration > 0.1:
                x0_item, x1_item = boundary_plotitems(spike_oscillations_high_beta_times, spike_oscillations_high_beta_wave, start, duration)
                plot3.addItem(pg.FillBetweenItem(x0_item, x1_item, brush=pg.mkBrush('lightsteelblue')))

    plot3.getAxis('bottom').setLabel('Time (s)')
    plot0.getAxis('left').setLabel('Theta')
    plot1.getAxis('left').setLabel('Alpha')
    plot2.getAxis('left').setLabel('Low Beta')
    plot3.getAxis('left').setLabel('High Beta')
    plot0.setTitle('Spiketrain Oscillations Waveforms')
    layout.setBackground('white')

    # ymax = np.max([spike_oscillations_theta_wave, spike_oscillations_alpha_wave, spike_oscillations_low_beta_wave, spike_oscillations_high_beta_wave])
    # ymin = np.min([spike_oscillations_theta_wave, spike_oscillations_alpha_wave, spike_oscillations_low_beta_wave, spike_oscillations_high_beta_wave])
    # xmin = spike_oscillations_theta_times[0]
    # xmax = spike_oscillations_theta_times[-1]
    # # plot0.setXRange(xmin, xmax)
    # plot1.setXRange(xmin, xmax)
    # plot2.setXRange(xmin, xmax)
    # plot3.setXRange(xmin, xmax)
    # plot0.setYRange(ymin, ymax)
    # plot1.setYRange(ymin, ymax)
    # plot2.setYRange(ymin, ymax)
    # plot3.setYRange(ymin, ymax)

    plot1.setXLink(plot0)
    plot2.setXLink(plot0)
    plot3.setXLink(plot0)
    plot1.setYLink(plot0)
    plot2.setYLink(plot0)
    plot3.setYLink(plot0)

    app.generate_plot_window(layout, 'spikeoscillations')


def button_tab_features_lfpplot_clicked(app):
    pass


def button_tab_features_autocorrelationplot_clicked(app):
    pass