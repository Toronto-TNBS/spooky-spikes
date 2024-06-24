from PySide6.QtWidgets import *
from PySide6 import QtCore
from PySide6.QtCore import Slot
import pyqtgraph as pg
import numpy as np

with open(f'{__file__[:-7]}/styles/styles.css', 'r') as file:
    STYLESHEET = file.read()


class App(QApplication):

    def __init__(self):
        super().__init__()

        self.setStyleSheet(STYLESHEET)

        self.window_main, self.grid_main = self.init_main_window()
        self.window_main.setWindowTitle('Spooky Spikes')

        self.frame_control, self.grid_control = self.init_control_frame()
        self.grid_main.addWidget(self.frame_control, 0, 0)

        self.frame_analysis, self.grid_analysis = self.init_analysis_frame()
        self.grid_main.addWidget(self.frame_analysis, 0, 1, 1, 2)


        ## Control Frame Widgets
        self.header_file = QLabel('Source Data')
        self.header_file.setProperty('class', 'header')
        self.grid_control.addWidget(self.header_file, 0, 0, 1, 3)

        self.label_file = QLabel('')
        self.label_file.setStyleSheet('background-color: white;')
        self.grid_control.addWidget(self.label_file, 1, 0, 1, 2)

        self.button_file = QPushButton('Choose file')
        self.grid_control.addWidget(self.button_file, 1, 2, 1, 1)

        self.header_channel = QLabel('Source Channel')
        self.header_channel.setProperty('class', 'header')
        self.grid_control.addWidget(self.header_channel, 2, 0, 1, 3)

        self.dropdown_channel = QComboBox()
        self.grid_control.addWidget(self.dropdown_channel, 3, 0, 1, 3)

        self.header_processing = QLabel('Processing')
        self.header_processing.setProperty('class', 'header')
        self.grid_control.addWidget(self.header_processing, 4, 0, 1, 3)

        self.checkbox_filtering = QCheckBox('Filtering')
        self.grid_control.addWidget(self.checkbox_filtering, 5, 0, 1, 3)

        self.checkbox_spikesorting = QCheckBox('Spike Sorting')
        self.grid_control.addWidget(self.checkbox_spikesorting, 6, 0, 1, 3)

        self.check_invertsegment = QCheckBox('Invert Segment')
        self.grid_control.addWidget(self.check_invertsegment, 7, 0, 1, 3)

        self.header_export = QLabel('Export')
        self.header_export.setProperty('class', 'header')
        self.grid_control.addWidget(self.header_export, 8, 0, 1, 3)

        self.dropdown_export = QComboBox()
        self.grid_control.addWidget(self.dropdown_export, 9, 0, 1, 3)

        self.grid_control.setRowStretch(self.grid_control.rowCount(), 1)    # Set empty row as stretching, so it takes up remainder of space at bottom.


        ## Analysis Frame Widgets
        self.grid_tab_main = QGridLayout()
        self.grid_tab_spikesorting = QGridLayout()
        self.grid_tab_features = QGridLayout()

        self.tabs = QTabWidget()
        self.grid_analysis.addWidget(self.tabs, 0, 0)

        widget_tab_main = QWidget()
        widget_tab_main.setLayout(self.grid_tab_main)
        self.tabs.addTab(widget_tab_main, 'Processing')

        widget_tab_spikesorting = QWidget()
        widget_tab_spikesorting.setLayout(self.grid_tab_spikesorting)
        self.tabs.addTab(widget_tab_spikesorting, 'Spike Sorting')

        widget_tab_features = QWidget()
        widget_tab_features.setLayout(self.grid_tab_features)
        self.tabs.addTab(widget_tab_features, 'Features')

        # Main Tab
        self.subheader_tab_main_filtering = QLabel('Digital Filtering')
        self.subheader_tab_main_filtering.setProperty('class', 'subheader')
        self.grid_tab_main.addWidget(self.subheader_tab_main_filtering, 0, 0, 1, 4)

        self.label_tab_main_hpcutoff = QLabel('Highpass cut-off (Hz)')
        self.label_tab_main_hpcutoff.setProperty('class', 'label')
        self.grid_tab_main.addWidget(self.label_tab_main_hpcutoff, 1, 0)

        self.label_tab_main_lpcutoff = QLabel('Lowpass cut-off (Hz)')
        self.label_tab_main_lpcutoff.setProperty('class', 'label')
        self.grid_tab_main.addWidget(self.label_tab_main_lpcutoff, 1, 1)
        
        self.entry_tab_main_hpcutoff = QLineEdit()
        self.grid_tab_main.addWidget(self.entry_tab_main_hpcutoff, 2, 0)
        
        self.entry_tab_main_lpcutoff = QLineEdit()
        self.grid_tab_main.addWidget(self.entry_tab_main_lpcutoff, 2, 1)
        
        self.button_tab_main_filtering = QPushButton('Set')
        self.grid_tab_main.addWidget(self.button_tab_main_filtering, 2, 2)
        
        self.label_tab_main_filterstatus = QLabel('HP: None\nLP: None')
        self.grid_tab_main.addWidget(self.label_tab_main_filterstatus, 2, 3)

        self.subheader_tab_main_threshold = QLabel('Threshold')
        self.subheader_tab_main_threshold.setProperty('class', 'subheader')
        self.grid_tab_main.addWidget(self.subheader_tab_main_threshold, 3, 0, 1, 4)

        self.label_tab_main_madfactor = QLabel('Scaling Factor')
        self.label_tab_main_madfactor.setProperty('class', 'label')
        self.grid_tab_main.addWidget(self.label_tab_main_madfactor, 4, 0)

        self.entry_tab_main_madfactor = QLineEdit()
        self.grid_tab_main.addWidget(self.entry_tab_main_madfactor, 5, 0)

        self.button_tab_main_threshold = QPushButton('Set')
        self.grid_tab_main.addWidget(self.button_tab_main_threshold, 5, 1)

        self.label_tab_main_thresholdstatus = QLabel('Threshold: None (Factor: None)')
        self.grid_tab_main.addWidget(self.label_tab_main_thresholdstatus, 5, 2, 1, 2)

        self.subheader_tab_main_plot = QLabel('Plot Display')
        self.subheader_tab_main_plot.setProperty('class', 'subheader')
        self.grid_tab_main.addWidget(self.subheader_tab_main_plot, 6, 0, 1, 3)

        self.subheader_tab_main_cluster = QLabel('Select Cluster')
        self.subheader_tab_main_cluster.setProperty('class', 'subheader')
        self.grid_tab_main.addWidget(self.subheader_tab_main_cluster, 6, 3)

        self.check_tab_main_spiketrain = QCheckBox('Spike Train')
        self.grid_tab_main.addWidget(self.check_tab_main_spiketrain, 7, 0)

        self.check_tab_main_eventtimes = QCheckBox('Event Times')
        self.grid_tab_main.addWidget(self.check_tab_main_eventtimes, 7, 1)

        self.check_tab_main_thresholdbar = QCheckBox('Threshold Bar')
        self.grid_tab_main.addWidget(self.check_tab_main_thresholdbar, 7, 2)

        self.dropdown_tab_main_cluster = QComboBox()
        self.grid_tab_main.addWidget(self.dropdown_tab_main_cluster, 7, 3)

        self.plot1_tab_main, self.plot2_tab_main, self.plot_layout_tab_main = self.init_tab_main_plot()
        self.frame_plot_tab_main = QFrame()
        self.frame_plot_tab_main.setProperty('class', 'frame-plot')
        self.frame_plot_tab_main.setFrameShape(QFrame.Box)
        self.grid_frame_plot_tab_main = QGridLayout()
        self.frame_plot_tab_main.setLayout(self.grid_frame_plot_tab_main)
        self.grid_frame_plot_tab_main.addWidget(self.plot_layout_tab_main, 0, 0)
        self.grid_tab_main.addWidget(self.frame_plot_tab_main, 8, 0, 1, 4)
    
        
        # Spike Sorting Tab
        self.subheader_tab_spikesorting_clusters = QLabel('Number of Clusters')
        self.subheader_tab_spikesorting_clusters.setProperty('class', 'subheader')
        self.grid_tab_spikesorting.addWidget(self.subheader_tab_spikesorting_clusters, 0, 0)

        self.dropdown_tab_spikesorting_clusters = QComboBox()
        self.grid_tab_spikesorting.addWidget(self.dropdown_tab_spikesorting_clusters, 1, 0)

        self.plot_tab_spikesorting, self.plot_layout_tab_spikesorting = self.init_tab_spikesorting_plot()
        self.frame_plot_tab_spikesorting = QFrame()
        self.frame_plot_tab_spikesorting.setFrameShape(QFrame.Box)
        self.frame_plot_tab_spikesorting.setProperty('class', 'frame-plot')
        self.grid_frame_plot_tab_spikesorting = QGridLayout()
        self.frame_plot_tab_spikesorting.setLayout(self.grid_frame_plot_tab_spikesorting)
        self.grid_frame_plot_tab_spikesorting.addWidget(self.plot_layout_tab_spikesorting)
        self.grid_tab_spikesorting.addWidget(self.frame_plot_tab_spikesorting, 2, 0, 1, 3)

        # Features Tab

        # Features Patterned Frame
        self.subheader_tab_features_patterned = QLabel('Signal Patterns')
        self.subheader_tab_features_patterned.setProperty('class', 'subheader')
        self.grid_tab_features.addWidget(self.subheader_tab_features_patterned, 0, 0)

        self.frame_tab_features_patterned, self.grid_tab_features_patterned = self.init_features_frame_grid()
        self.grid_tab_features.addWidget(self.frame_tab_features_patterned, 1, 0)
        
        self.label_tab_features_patterned_firingrate = QLabel('Firing rate (Hz): ')
        self.grid_tab_features_patterned.addWidget(self.label_tab_features_patterned_firingrate, 0, 0)
        self.label_tab_features_patterned_firingrate_status = QLabel('')
        self.label_tab_features_patterned_firingrate_status.setProperty('class', 'label-features-status')
        self.grid_tab_features_patterned.addWidget(self.label_tab_features_patterned_firingrate_status, 0, 1)

        self.label_tab_features_patterned_burstindex = QLabel('Burst index:')
        self.grid_tab_features_patterned.addWidget(self.label_tab_features_patterned_burstindex, 1, 0)
        self.label_tab_features_patterned_burstindex_status = QLabel('')
        self.label_tab_features_patterned_burstindex_status.setProperty('class', 'label-features-status')
        self.grid_tab_features_patterned.addWidget(self.label_tab_features_patterned_burstindex_status, 1, 1)

        self.label_tab_features_patterned_cov = QLabel('Coefficient of variation:')
        self.grid_tab_features_patterned.addWidget(self.label_tab_features_patterned_cov, 2, 0)
        self.label_tab_features_patterned_cov_status = QLabel('')
        self.label_tab_features_patterned_cov_status.setProperty('class', 'label-features-status')
        self.grid_tab_features_patterned.addWidget(self.label_tab_features_patterned_cov_status, 2, 1)

        # Features Quality Metrics Frame
        self.subheader_tab_features_qualitymetrics = QLabel('Quality Metrics')
        self.subheader_tab_features_qualitymetrics.setProperty('class', 'subheader')
        self.grid_tab_features.addWidget(self.subheader_tab_features_qualitymetrics, 0, 1)

        self.frame_tab_features_qualitymetrics, self.grid_tab_features_qualitymetrics = self.init_features_frame_grid()
        self.grid_tab_features.addWidget(self.frame_tab_features_qualitymetrics, 1, 1)

        self.label_tab_features_qualitymetrics_snr = QLabel('Signal to noise ratio:')
        self.grid_tab_features_qualitymetrics.addWidget(self.label_tab_features_qualitymetrics_snr, 0, 0)
        self.label_tab_features_qualitymetrics_snr_status = QLabel('')
        self.label_tab_features_qualitymetrics_snr_status.setProperty('class', 'label-features-status')
        self.grid_tab_features_qualitymetrics.addWidget(self.label_tab_features_qualitymetrics_snr_status, 0, 1)

        self.label_tab_features_qualitymetrics_isi = QLabel('Percent ISI violations:')
        self.grid_tab_features_qualitymetrics.addWidget(self.label_tab_features_qualitymetrics_isi, 1, 0)
        self.label_tab_features_qualitymetrics_isi_status = QLabel('')
        self.label_tab_features_qualitymetrics_isi_status.setProperty('class', 'label-features-status')
        self.grid_tab_features_qualitymetrics.addWidget(self.label_tab_features_qualitymetrics_isi_status, 1, 1)

        self.label_tab_features_qualitymetrics_silhouette = QLabel('Silhouette score:')
        self.grid_tab_features_qualitymetrics.addWidget(self.label_tab_features_qualitymetrics_silhouette, 2, 0)
        self.label_tab_features_qualitymetrics_silhouette_status = QLabel('')
        self.label_tab_features_qualitymetrics_silhouette_status.setProperty('class', 'label-features-status')
        self.grid_tab_features_qualitymetrics.addWidget(self.label_tab_features_qualitymetrics_silhouette_status, 2, 1)

        # Features ISI Plot Frame
        self.frame_tab_features_isiplot, self.grid_tab_features_isiplot = self.init_features_frame_grid()
        self.grid_tab_features.addWidget(self.frame_tab_features_isiplot, 1, 2)

        self.label_tab_features_isiplot = QLabel('ISI Plot')
        self.grid_tab_features_isiplot.addWidget(self.label_tab_features_isiplot, 0, 0)

        self.button_tab_features_isiplot = QPushButton('Generate')
        self.grid_tab_features_isiplot.addWidget(self.button_tab_features_isiplot, 1, 0)

        # Features Spiketrain Oscillations Frame
        self.subheader_tab_features_spiketrain = QLabel('Spiketrain Oscillations')
        self.subheader_tab_features_spiketrain.setProperty('class', 'subheader')
        self.grid_tab_features.addWidget(self.subheader_tab_features_spiketrain, 2, 0)

        self.frame_tab_features_spiketrain, self.grid_tab_features_spiketrain = self.init_features_frame_grid()
        self.grid_tab_features.addWidget(self.frame_tab_features_spiketrain, 3, 0, 1, 2)

        self.label_tab_features_spiketrain_power = QLabel('Power:')
        self.grid_tab_features_spiketrain.addWidget(self.label_tab_features_spiketrain_power, 1, 0)

        self.label_tab_features_spiketrain_burstduration = QLabel('Burst duration:')
        self.grid_tab_features_spiketrain.addWidget(self.label_tab_features_spiketrain_burstduration, 2, 0)

        self.label_tab_features_spiketrain_theta = QLabel('Theta')
        self.grid_tab_features_spiketrain.addWidget(self.label_tab_features_spiketrain_theta, 0, 1)
        self.label_tab_features_spiketrain_thetapower_status = QLabel('')
        self.label_tab_features_spiketrain_thetapower_status.setProperty('class', 'label-features-status')
        self.grid_tab_features_spiketrain.addWidget(self.label_tab_features_spiketrain_thetapower_status, 1, 1)
        self.label_tab_features_spiketrain_thetaburstduration_status = QLabel('')
        self.label_tab_features_spiketrain_thetaburstduration_status.setProperty('class', 'label-features-status')
        self.grid_tab_features_spiketrain.addWidget(self.label_tab_features_spiketrain_thetaburstduration_status, 2, 1)

        self.label_tab_features_spiketrain_alpha = QLabel('Alpha')
        self.grid_tab_features_spiketrain.addWidget(self.label_tab_features_spiketrain_alpha, 0, 2)
        self.label_tab_features_spiketrain_alphapower_status = QLabel('')
        self.label_tab_features_spiketrain_alphapower_status.setProperty('class', 'label-features-status')
        self.grid_tab_features_spiketrain.addWidget(self.label_tab_features_spiketrain_alphapower_status, 1, 2)
        self.label_tab_features_spiketrain_alphaburstduration_status = QLabel('')
        self.label_tab_features_spiketrain_alphaburstduration_status.setProperty('class', 'label-features-status')
        self.grid_tab_features_spiketrain.addWidget(self.label_tab_features_spiketrain_alphaburstduration_status, 2, 2)

        self.label_tab_features_spiketrain_lowbeta = QLabel('Low-beta')
        self.grid_tab_features_spiketrain.addWidget(self.label_tab_features_spiketrain_lowbeta, 0, 3)
        self.label_tab_features_spiketrain_lowbetapower_status = QLabel('')
        self.label_tab_features_spiketrain_lowbetapower_status.setProperty('class', 'label-features-status')
        self.grid_tab_features_spiketrain.addWidget(self.label_tab_features_spiketrain_lowbetapower_status, 1, 3)
        self.label_tab_features_spiketrain_lowbetaburstduration_status = QLabel('')
        self.label_tab_features_spiketrain_lowbetaburstduration_status.setProperty('class', 'label-features-status')
        self.grid_tab_features_spiketrain.addWidget(self.label_tab_features_spiketrain_lowbetaburstduration_status, 2, 3)

        self.label_tab_features_spiketrain_highbeta = QLabel('High-beta')
        self.grid_tab_features_spiketrain.addWidget(self.label_tab_features_spiketrain_highbeta, 0, 4)
        self.label_tab_features_spiketrain_highbetapower_status = QLabel('')
        self.label_tab_features_spiketrain_highbetapower_status.setProperty('class', 'label-features-status')
        self.grid_tab_features_spiketrain.addWidget(self.label_tab_features_spiketrain_highbetapower_status, 1, 4)
        self.label_tab_features_spiketrain_highbetaburstduration_status = QLabel('')
        self.label_tab_features_spiketrain_highbetaburstduration_status.setProperty('class', 'label-features-status')
        self.grid_tab_features_spiketrain.addWidget(self.label_tab_features_spiketrain_highbetaburstduration_status, 2, 4)

        # Features Spiketrain Oscillations Plot
        self.frame_tab_features_spikeoscillationsplot, self.grid_tab_features_spikeoscillationsplot = self.init_features_frame_grid()
        self.grid_tab_features.addWidget(self.frame_tab_features_spikeoscillationsplot, 3, 2)

        self.label_tab_features_spikeoscillationsplot = QLabel('Oscillations Plot')
        self.grid_tab_features_spikeoscillationsplot.addWidget(self.label_tab_features_spikeoscillationsplot, 0, 0)
        
        self.button_tab_features_spikeoscillationsplot = QPushButton('Generate')
        self.grid_tab_features_spikeoscillationsplot.addWidget(self.button_tab_features_spikeoscillationsplot, 1, 0)

        # Features LFP Power Frame
        self.subheader_tab_features_lfp = QLabel('LFP Power')
        self.subheader_tab_features_lfp.setProperty('class', 'subheader')
        self.grid_tab_features.addWidget(self.subheader_tab_features_lfp, 4, 0)

        self.frame_tab_features_lfp, self.grid_tab_features_lfp = self.init_features_frame_grid()
        self.grid_tab_features.addWidget(self.frame_tab_features_lfp, 5, 0, 1, 2)

        self.label_tab_features_lfp_power = QLabel('Power:')
        self.grid_tab_features_lfp.addWidget(self.label_tab_features_lfp_power, 1, 0)

        self.label_tab_features_lfp_theta = QLabel('Theta')
        self.grid_tab_features_lfp.addWidget(self.label_tab_features_lfp_theta, 0, 1)
        self.label_tab_features_lfp_thetapower_status = QLabel('')
        self.label_tab_features_lfp_thetapower_status.setProperty('class', 'label-features-status')
        self.grid_tab_features_lfp.addWidget(self.label_tab_features_lfp_thetapower_status, 1, 1)

        self.label_tab_features_lfp_alpha = QLabel('Alpha')
        self.grid_tab_features_lfp.addWidget(self.label_tab_features_lfp_alpha, 0, 2)
        self.label_tab_features_lfp_alphapower_status = QLabel('')
        self.label_tab_features_lfp_alphapower_status.setProperty('class', 'label-features-status')
        self.grid_tab_features_lfp.addWidget(self.label_tab_features_lfp_alphapower_status, 1, 2)

        self.label_tab_features_lfp_lowbeta = QLabel('Low-beta')
        self.grid_tab_features_lfp.addWidget(self.label_tab_features_lfp_lowbeta, 0, 3)
        self.label_tab_features_lfp_lowbetapower_status = QLabel('')
        self.label_tab_features_lfp_lowbetapower_status.setProperty('class', 'label-features-status')
        self.grid_tab_features_lfp.addWidget(self.label_tab_features_lfp_lowbetapower_status, 1, 3)

        self.label_tab_features_lfp_highbeta = QLabel('High-beta')
        self.grid_tab_features_lfp.addWidget(self.label_tab_features_lfp_highbeta, 0, 4)
        self.label_tab_features_lfp_highbetapower_status = QLabel('')
        self.label_tab_features_lfp_highbetapower_status.setProperty('class', 'label-features-status')
        self.grid_tab_features_lfp.addWidget(self.label_tab_features_lfp_highbetapower_status, 1, 4)

        # Features LFP Plot
        self.frame_tab_features_lfpplot, self.grid_tab_features_lfpplot = self.init_features_frame_grid()
        self.grid_tab_features.addWidget(self.frame_tab_features_lfpplot, 5, 2)

        self.label_tab_features_lfpplot = QLabel('LFP Plot')
        self.grid_tab_features_lfpplot.addWidget(self.label_tab_features_lfpplot, 0, 0)
        
        self.button_tab_features_lfpplot = QPushButton('Generate')
        self.grid_tab_features_lfpplot.addWidget(self.button_tab_features_lfpplot, 1, 0)

        # Features Autocorrelation Frame
        self.subheader_tab_features_autocorrelation = QLabel('Autocorrelation')
        self.subheader_tab_features_autocorrelation.setProperty('class', 'subheader')
        self.grid_tab_features.addWidget(self.subheader_tab_features_autocorrelation, 6, 0)

        self.frame_tab_features_autocorrelationplot, self.grid_tab_features_autocorrelationplot = self.init_features_frame_grid()
        self.grid_tab_features.addWidget(self.frame_tab_features_autocorrelationplot, 7, 0)

        self.label_tab_features_autocorrelationplot = QLabel('Autocorrelation Plot')
        self.grid_tab_features_autocorrelationplot.addWidget(self.label_tab_features_autocorrelationplot, 0, 0)

        self.button_tab_features_autocorrelationplot = QPushButton('Generate')
        self.grid_tab_features_autocorrelationplot.addWidget(self.button_tab_features_autocorrelationplot, 1, 0)

        self.grid_tab_features.setRowStretch(self.grid_tab_features.rowCount(), 1)    # Stretch empty row at end to fill bottom space.


        
        self.exec()
    

    def init_main_window(self):
        window = QWidget()
        window.resize(900, 600)

        grid = QGridLayout()
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setHorizontalSpacing(0)

        window.setLayout(grid)
        window.show()
        return window, grid
    

    def init_control_frame(self):
        frame = QFrame()
        frame.setObjectName('control')
        frame.setFrameShape(QFrame.Box)    # Edit this. Make it look nicer.
        grid = QGridLayout()
        frame.setLayout(grid)
        return frame, grid
    

    def init_analysis_frame(self):
        frame = QFrame()
        frame.setObjectName('analysis')
        grid = QGridLayout()
        grid.setContentsMargins(0, 0, 0, 0)
        frame.setLayout(grid)
        return frame, grid


    def init_features_frame_grid(self):
        frame = QFrame()
        frame.setProperty('class', 'frame-features')
        frame.setFrameShape(QFrame.Box)
        grid = QGridLayout()
        frame.setLayout(grid)
        return frame, grid
    

    def init_tab_main_plot(self):
        data = np.random.normal(0, 1, int(1e3))

        layout = pg.GraphicsLayoutWidget()
        plot1 = layout.addPlot(row=0, col=0)
        plot2 = layout.addPlot(y=data, row=1, col=0, rowspan=2, pen=pg.mkPen('cornflowerblue', width=1.5))

        for i in np.linspace(0, len(data), 100):
            event = pg.InfiniteLine(pos=i, pen=pg.mkPen('r'))
            plot1.addItem(event)
        
        plot2.setClipToView(True)
        plot2.setDownsampling(True)
        
        layout.setBackground('#EDEDED')
        plot1.setProperty('class', 'plot-main')
        plot2.setProperty('class', 'plot-main')

        layout.ci.layout.setRowStretchFactor(0, 1)
        layout.ci.layout.setRowStretchFactor(1, 10)

        plot1.showAxes(selection=['bottom', 'left'], showValues=False)
        plot1.getAxis('left').setLabel('Spike Train')
        plot1.getAxis('left').setTicks([])
        plot1.getAxis('right').setTicks([])
        plot1.getAxis('top').setTicks([])
        plot1.getAxis('bottom').setTicks([])
        plot1.hideAxis('right')
        plot1.hideAxis('top')

        plot2.getAxis('bottom').setLabel('Time (s)')
        plot2.getAxis('left').setLabel('Magnitude')

        plot1.setXLink(plot2)
        
        return plot1, plot2, layout
    

    def init_tab_spikesorting_plot(self):
        data = np.random.normal(0, 1, int(1e2))
        layout = pg.GraphicsLayoutWidget()
        plot = layout.addPlot(row=0, col=0)
        item = pg.ScatterPlotItem(x=data, y=data, pen=pg.mkPen('cornflowerblue', width=1.25), brush=pg.mkBrush('cornflowerblue'))    # Pen for border, Brush for fill.
        plot.addItem(item)

        layout.setBackground('#EDEDED')
        plot.setClipToView(True)
        plot.setDownsampling(True)

        plot.getAxis('bottom').setLabel('Principal Component 1')
        plot.getAxis('left').setLabel('Principal Component 2')

        return plot, layout



App()