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
    
    app.check_spikesorting.setEnabled(True)
    app.check_filtering.setEnabled(True)
    app.check_invertsegment.setEnabled(True)
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
        return
    
    current_signal = app.channeldata.get_current_signal()
    time_vector = app.channeldata.get_time_vector()
    
    threshold = analysis.get_spike_threshold(current_signal, factor)
    
    app.channeldata.threshold_factor = factor
    app.channeldata.threshold = threshold
    app.label_tab_main_thresholdstatus.setText(f'Threshold: {round(app.channeldata.threshold, 3)} (Factor: {app.channeldata.threshold_factor})')
    
    spike_indices = analysis.get_spike_indices(current_signal, app.channeldata.threshold, app.channeldata.fs)
    app.channeldata.spike_indices = spike_indices

    # Update main plot
    def vertical_line_symbol():
        symbol = QtGui.QPainterPath()
        symbol.moveTo(0, -10)
        symbol.lineTo(0, 10)
        return symbol
        

    misc.remove_plot_item(plotitem=app.plot2_tab_main, itemtype=pg.ScatterPlotItem)    # Remove previously plotted spikes if they exist.
    peaks_scatter_item = pg.ScatterPlotItem(x=time_vector[spike_indices], y=current_signal[spike_indices], pen=pg.mkPen('r'), brush=pg.mkBrush('r'), size=5)
    app.plot2_tab_main.addItem(peaks_scatter_item)
    threshold_bar_item = pg.InfiniteLine(pos=-app.channeldata.threshold, angle=0, pen=pg.mkPen('k', width=1, style=QtCore.Qt.DashLine))
    app.plot2_tab_main.addItem(threshold_bar_item)
    # misc.plot_spike_events(plotitem=app.plot1_tab_main, spike_times=time_vector[spike_indices])
    events_scatter_item = pg.ScatterPlotItem(x=time_vector[spike_indices], y=[0]*len(spike_indices), pen=pg.mkPen('r'), brush=pg.mkBrush('r'))
    app.plot1_tab_main.addItem(events_scatter_item)

    # Will this execute the slots connected to the checkboxes?
    app.check_tab_main_spiketrain.setChecked(True)
    app.check_tab_main_eventtimes.setChecked(True)
    app.check_tab_main_thresholdbar.setChecked(True)
    app.check_tab_main_spiketrain.setEnabled(True)
    app.check_tab_main_eventtimes.setEnabled(True)
    app.check_tab_main_thresholdbar.setEnabled(True)
