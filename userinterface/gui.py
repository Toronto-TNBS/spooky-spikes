from tkinter import *
import tkinter.ttk as ttk
import tkinter.filedialog as tkf
import tkinter.messagebox as tkm
from userinterface.gui_styles import GUIStyles
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from spikedata.readsavedata import ReadSpike
from spikedata.spikeanalysis import SpikeAnalysis
import numpy as np
import sys
import base64
from images import base64_image

read = ReadSpike()
spike = SpikeAnalysis()


class App(GUIStyles):

    def __init__(self):
        super().__init__()


        self.title('Spooky Spikes Dashboard')

        img = base64_image.logo
        img = base64.b64decode(img)
        img = PhotoImage(data=img)
        self.wm_iconphoto(True, img)
        # self.iconbitmap('images/tnbs_logo.ico')

        self.left_panel = Canvas(self, bg=self.left_panel_bg, highlightthickness=0)
        self.left_panel.grid(column=0, row=0, columnspan=2, rowspan=25, sticky='nsew')

        self.right_panel = Canvas(self, bg=self.right_panel_bg, highlightthickness=0)
        self.right_panel.grid(column=2, row=0, rowspan=25, sticky='nsew')


        self.header_source_data = ttk.Label(
            text='Source Data:',
            style='headers.TLabel',
            background=self.left_panel_bg,
        )
        self.header_source_data.grid(column=0, row=0, sticky='w', padx=self.header_align_pad)

        self.label_file = ttk.Label(
            text='',
            background='white',
            width=32,
            style='file.TLabel'
        )
        self.label_file.grid(column=0, row=1, columnspan=2, sticky='ew', padx=self.widget_align_pad)

        self.button_file = ttk.Button(
            text='Choose file',
            style='TButton',
            command=self.button_file_press
        )
        self.button_file.grid(column=1, row=1, sticky='e', padx=self.widget_align_pad)

        self.header_channel = ttk.Label(
            text='Source Channel:',
            style='headers.TLabel',
            background=self.left_panel_bg,
        )
        self.header_channel.grid(column=0, row=2, sticky='w', padx=self.header_align_pad)

        # read open channels, or all of them, regardless of whether off or on
        self.channels_available = []
        self.dropdown_selection = StringVar()
        self.dropdown_menu = None
        self.dropdown_channel = ttk.Menubutton(
            direction='below',
            text='Channels',
            style='TMenubutton',
            width=40
        )
        self.dropdown_channel.grid(column=0, row=3, columnspan=2, sticky='we', padx=self.widget_align_pad)

        self.header_parameters = ttk.Label(
            text='Parameters:',
            style='headers.TLabel',
            background=self.left_panel_bg,
        )
        self.header_parameters.grid(column=0, row=4, sticky='w', padx=self.header_align_pad)

        self.check_filtering_status = IntVar()
        self.check_filtering = ttk.Checkbutton(
            command=self.check_filtering_choice,
            style='left_check.TCheckbutton',
            text='Filtering',
            variable=self.check_filtering_status,
        )
        self.check_filtering.grid(column=0, row=5, sticky='w', padx=self.widget_align_pad)
        # use padding here, much better than style

        self.check_spikesorting_status = IntVar()
        self.check_spikesorting = ttk.Checkbutton(
            command=self.check_spikesorting_choice,
            style='left_check.TCheckbutton',
            text='Spike Sorting',
            variable=self.check_spikesorting_status
        )
        self.check_spikesorting.grid(column=0, row=6, sticky='w', padx=self.widget_align_pad)

        self.check_invertsegment_status = IntVar()
        self.check_invertsegment = ttk.Checkbutton(
            command=self.check_invertsegment_choice,
            style='left_check.TCheckbutton',
            text='Invert Threshold',
            variable=self.check_invertsegment_status
        )
        self.check_invertsegment.grid(column=0, row=7, sticky='w', padx=self.widget_align_pad)

        self.header_quality_metrics = ttk.Label(
            text='Quality Metrics:',
            style='headers.TLabel',
            background=self.left_panel_bg
        )
        self.header_quality_metrics.grid(column=0, row=8, sticky='w', padx=self.header_align_pad)

        self.label_snr = ttk.Label(
            text='Signal to noise ratio:',
            style='outputs.TLabel',
            background=self.left_panel_bg,
            width=self.label_output_name_width
        )
        self.label_snr.grid(column=0, row=9, sticky='w', padx=self.widget_align_pad, pady=self.label_output_pady)

        self.label_snr_output = ttk.Label(
            text='',
            style='outputs.TLabel',
            background='white',
            width=self.label_output_width
        )
        self.label_snr_output.grid(column=1, row=9, sticky='we', padx=self.widget_align_pad, pady=self.label_output_pady)

        self.label_interspike = ttk.Label(
            text='Percent ISI violations:',
            style='outputs.TLabel',
            background=self.left_panel_bg,
            width=self.label_output_name_width
        )
        self.label_interspike.grid(column=0, row=10, sticky='w', padx=self.widget_align_pad, pady=self.label_output_pady)

        self.label_interspike_output = ttk.Label(
            text='',
            style='outputs.TLabel',
            background='white',
            width=self.label_output_width
        )
        self.label_interspike_output.grid(column=1, row=10, sticky='we', padx=self.widget_align_pad, pady=self.label_output_pady)

        self.label_silhouette = ttk.Label(
            text='Silhouette Score:',
            style='outputs.TLabel',
            background=self.left_panel_bg,
            width=self.label_output_name_width
        )
        self.label_silhouette.grid(column=0, row=11, sticky='w', padx=self.widget_align_pad, pady=self.label_output_pady)

        self.label_silhouette_output = ttk.Label(
            text='',
            style='outputs.TLabel',
            background='white',
            width=self.label_output_width
        )
        self.label_silhouette_output.grid(column=1, row=11, sticky='we', padx=self.widget_align_pad, pady=self.label_output_pady)

        self.header_features = ttk.Label(
            text='Features:',
            style='headers.TLabel',
            background=self.left_panel_bg
        )
        self.header_features.grid(column=0, row=12, sticky='w', padx=self.header_align_pad)

        self.label_firing_rate = ttk.Label(
            text='Firing rate (Hz):',
            style='outputs.TLabel',
            background=self.left_panel_bg,
            width=self.label_output_name_width
        )
        self.label_firing_rate.grid(column=0, row=13, sticky='w', padx=self.widget_align_pad, pady=self.label_output_pady)

        self.label_firing_rate_output = ttk.Label(
            text='',
            style='outputs.TLabel',
            background='white',
            width=self.label_output_width
        )
        self.label_firing_rate_output.grid(column=1, row=13, sticky='we', padx=self.widget_align_pad, pady=self.label_output_pady)

        self.label_burst_index = ttk.Label(
            text='Burst index:',
            style='outputs.TLabel',
            background=self.left_panel_bg,
            width=self.label_output_name_width
        )
        self.label_burst_index.grid(column=0, row=14, sticky='w', padx=self.widget_align_pad, pady=self.label_output_pady)

        self.label_burst_index_output = ttk.Label(
            text='',
            style='outputs.TLabel',
            background='white',
            width=self.label_output_width
        )
        self.label_burst_index_output.grid(column=1, row=14, sticky='we', padx=self.widget_align_pad, pady=self.label_output_pady)

        self.label_variation_coefficient = ttk.Label(
            text='Coefficient of variation:',
            style='outputs.TLabel',
            background=self.left_panel_bg,
            width=self.label_output_name_width
        )
        self.label_variation_coefficient.grid(column=0, row=15, sticky='w', padx=self.widget_align_pad, pady=self.label_output_pady)

        self.label_variation_coefficient_output = ttk.Label(
            text='',
            style='outputs.TLabel',
            background='white',
            width=self.label_output_width
        )
        self.label_variation_coefficient_output.grid(column=1, row=15, sticky='we', padx=self.widget_align_pad, pady=self.label_output_pady)

        self.label_theta = ttk.Label(
            text='Theta oscillatory power:',
            style='outputs.TLabel',
            background=self.left_panel_bg,
            width=self.label_output_name_width
        )
        self.label_theta.grid(column=0, row=17, sticky='w', padx=self.widget_align_pad, pady=self.label_output_pady)

        self.label_theta_output = ttk.Label(
            text='',
            style='outputs.TLabel',
            background='white',
            width=self.label_output_width
        )
        self.label_theta_output.grid(column=1, row=17, sticky='we', padx=self.widget_align_pad, pady=self.label_output_pady)

        self.label_alpha = ttk.Label(
            text='Alpha oscillatory power:',
            style='outputs.TLabel',
            background=self.left_panel_bg,
            width=self.label_output_name_width
        )
        self.label_alpha.grid(column=0, row=18, sticky='w', padx=self.widget_align_pad, pady=self.label_output_pady)

        self.label_alpha_output = ttk.Label(
            text='',
            style='outputs.TLabel',
            background='white',
            width=self.label_output_width
        )
        self.label_alpha_output.grid(column=1, row=18, sticky='we', padx=self.widget_align_pad, pady=self.label_output_pady)

        self.label_low_beta = ttk.Label(
            text='Low beta oscillatory power:',
            style='outputs.TLabel',
            background=self.left_panel_bg,
            width=self.label_output_name_width
        )
        self.label_low_beta.grid(column=0, row=19, sticky='w', padx=self.widget_align_pad, pady=self.label_output_pady)

        self.label_low_beta_output = ttk.Label(
            text='',
            style='outputs.TLabel',
            background='white',
            width=self.label_output_width
        )
        self.label_low_beta_output.grid(column=1, row=19, sticky='we', padx=self.widget_align_pad, pady=self.label_output_pady)

        self.label_high_beta = ttk.Label(
            text='High beta oscillatory power:',
            style='outputs.TLabel',
            background=self.left_panel_bg,
            width=self.label_output_name_width
        )
        self.label_high_beta.grid(column=0, row=20, sticky='w', padx=self.widget_align_pad, pady=self.label_output_pady)

        self.label_high_beta_output = ttk.Label(
            text='',
            style='outputs.TLabel',
            background='white',
            width=self.label_output_width
        )
        self.label_high_beta_output.grid(column=1, row=20, sticky='we', padx=self.widget_align_pad, pady=self.label_output_pady)

        self.header_export = ttk.Label(
            text='Export:',
            style='headers.TLabel',
            background=self.left_panel_bg
        )
        self.header_export.grid(column=0, row=22, sticky='w', padx=self.header_align_pad)

        self.dropdown_save_menu = Menu(tearoff=False)
        self.dropdown_save_menu.add_command(
            label='Save Properties',
            command=self.dropdown_save_properties,
            font=self.menu_item_font
        )
        self.dropdown_save_menu.add_command(
            label='Save Segment (All Spikes)',
            command=self.dropdown_save_segment_all_spikes,
            font=self.menu_item_font
        )
        self.dropdown_save_menu.add_command(
            label=f'Save Segment (Selected Cluster)',
            command=self.dropdown_save_segment,
            font=self.menu_item_font
        )
        self.dropdown_save = ttk.Menubutton(
            direction='below',
            text='Save Options',
            menu=self.dropdown_save_menu,
            style='TMenubutton',
            width=40
        )
        self.dropdown_save.grid(column=0, row=23, columnspan=2, sticky='we', padx=self.widget_align_pad, pady=5)

        self.num_rows = self.grid_size()[1]
        self.num_columns = self.grid_size()[0]
        self.dynamic_resize()

        self.notebook = ttk.Notebook(self)
        self.notebook.grid(column=2, row=0, rowspan=self.num_rows, sticky='nesw')
        self.tab_main = Frame(self.notebook, bg='white')
        self.tab_spikesorting = Frame(self.notebook, bg='white')
        self.tab_psd = Frame(self.notebook, bg='white')
        self.tab_oscillations = Frame(self.notebook, bg='white')
        self.tab_main.grid(column=2, row=0, sticky='nesw')
        self.tab_spikesorting.grid(column=2, row=0, sticky='nesw')
        self.tab_psd.grid(column=2, row=0, sticky='nesw')
        self.tab_oscillations.grid(column=2, row=0, sticky='nesw')
        self.notebook.add(self.tab_main, text='Main')
        self.notebook.add(self.tab_spikesorting, text='Spike Sorting')
        self.notebook.add(self.tab_psd, text='PSD')
        self.notebook.add(self.tab_oscillations, text='Oscillations')

        # Background tab section dividers
        self.filtering_bg = Canvas(self.tab_main, bg='white', height=120)
        self.filtering_bg.grid(column=1, row=0, columnspan=4, rowspan=3, sticky='nesw')

        self.threshold_bg = Canvas(self.tab_main, bg='white', height=120)
        self.threshold_bg.grid(column=1, row=3, columnspan=4, rowspan=3, sticky='nesw')

        self.plot_display_bg = Canvas(self.tab_main, bg='white', height=90)
        self.plot_display_bg.grid(column=1, row=6, columnspan=3, rowspan=2, sticky='nesw')

        self.cluster_select_bg = Canvas(self.tab_main, bg='white', height=90, width=1)
        self.cluster_select_bg.grid(column=4, row=6, rowspan=2, sticky='nesw')

        self.header_main_filtering = ttk.Label(
            self.tab_main,
            text='Signal Filtering:',
            style='right_headers.TLabel',
            background=self.right_panel_bg
        )
        self.header_main_filtering.grid(column=1, row=0, padx=self.header_align_pad, sticky='ws')

        self.label_set_cutoffs = ttk.Label(
            self.tab_main,
            text='Highpass cut-off (Hz):',
            background=self.right_panel_bg,
            style='outputs.TLabel'
        )
        self.label_set_cutoffs.grid(column=1, row=1, padx=self.widget_align_pad, sticky='w')

        self.entry_highpass_cutoff = ttk.Entry(
            self.tab_main,
            justify='left',
            state='disabled',
            width=20,
            style='unchecked.TEntry'
        )
        self.entry_highpass_cutoff.grid(column=1, row=2, sticky='w', padx=self.widget_align_pad)

        self.label_set_lowpass_cutoff = ttk.Label(
            self.tab_main,
            text='Lowpass cut-off (Hz):',
            background=self.right_panel_bg,
            style='outputs.TLabel'
        )
        self.label_set_lowpass_cutoff.grid(column=2, row=1, sticky='w')

        self.entry_lowpass_cutoff = ttk.Entry(
            self.tab_main,
            justify='left',
            state='disabled',
            width=20,
            style='unchecked.TEntry'
        )
        self.entry_lowpass_cutoff.grid(column=2, row=2, sticky='w')

        self.button_save_filter_cutoffs = ttk.Button(
            self.tab_main,
            text='Set',
            command=self.button_set_cutoffs_press,
            style='TButton'
        )
        self.button_save_filter_cutoffs.grid(column=3, row=2, sticky='w', padx=self.widget_align_pad)

        self.label_set_cutoffs = ttk.Label(
            self.tab_main,
            text='HP: 0 Hz\nLP: 0 Hz',
            style='right_outputs.TLabel',
            background=self.right_panel_bg,
        )
        self.label_set_cutoffs.grid(column=4, row=2, sticky='w')

        self.header_threshold = ttk.Label(
            self.tab_main,
            text='Threshold:',
            style='right_headers.TLabel',
            background=self.right_panel_bg
        )
        self.header_threshold.grid(column=1, row=3, sticky='ws', padx=self.header_align_pad)

        self.label_threshold_factor = ttk.Label(
            self.tab_main,
            text='MAD Factor:',
            style='outputs.TLabel',
            background=self.right_panel_bg
        )
        self.label_threshold_factor.grid(column=1, row=4, sticky='w', padx=self.widget_align_pad)

        self.entry_threshold_factor = ttk.Entry(
            self.tab_main,
            justify='left',
            state='disabled',
            width=20,
            style='unchecked.TEntry',
        )
        self.entry_threshold_factor.grid(column=1, row=5, sticky='w', padx=self.widget_align_pad)

        self.button_save_threshold_factor = ttk.Button(
            self.tab_main,
            text='Set',
            command=self.button_set_threshold_press,
            style='TButton'
        )
        self.button_save_threshold_factor.grid(column=2, row=5, sticky='w')

        self.label_threshold = ttk.Label(
            self.tab_main,
            text='Threshold: 0 V [Factor: 0]',
            style='right_outputs.TLabel',
            background=self.right_panel_bg
        )
        self.label_threshold.grid(column=3, row=5, columnspan=2, sticky='w', padx=self.widget_align_pad)

        self.header_plot_display = ttk.Label(
            self.tab_main,
            text='Plot Display:',
            style='right_headers.TLabel',
            background=self.right_panel_bg
        )
        self.header_plot_display.grid(column=1, row=6, sticky='ws', padx=self.header_align_pad)

        self.check_spiketrain_status = IntVar()
        self.check_spiketrain = ttk.Checkbutton(
            self.tab_main,
            text='Spike Train',
            command=self.check_spiketrain_choice,
            style='right_check.TCheckbutton',
            variable=self.check_spiketrain_status
        )
        self.check_spiketrain.grid(column=1, row=7, sticky='w', padx=self.widget_align_pad)

        self.check_thresholdbar_status = IntVar()
        self.check_thresholdbar = ttk.Checkbutton(
            self.tab_main,
            text='Threshold Bar',
            command=self.check_thresholdbar_choice,
            style='right_check.TCheckbutton',
            variable=self.check_thresholdbar_status
        )
        self.check_thresholdbar.grid(column=3, row=7, sticky='w', padx=self.widget_align_pad)

        self.check_eventpeaks_status = IntVar()
        self.check_eventpeaks = ttk.Checkbutton(
            self.tab_main,
            text='Event Peaks',
            command=self.check_eventpeaks_choice,
            style='right_check.TCheckbutton',
            variable=self.check_eventpeaks_status
        )
        self.check_eventpeaks.grid(column=2, row=7, sticky='w')

        self.header_select_cluster = ttk.Label(
            self.tab_main,
            text='Select Cluster:',
            style='right_headers.TLabel',
            background=self.right_panel_bg
        )
        self.header_select_cluster.grid(column=4, row=6, sticky='ws', padx=self.header_align_pad)

        self.dropdown_selected_colour_cluster = StringVar()
        self.dropdown_colour_cluster = ttk.Menubutton(
            self.tab_main,
            direction='below',
            text='Cluster',
            menu=None,
            style='clusters.TMenubutton',
            width=18
        )
        self.dropdown_colour_cluster.grid(column=4, row=7, columnspan=1, sticky='w', padx=self.widget_align_pad)

        # fillers to widen plot and shrink space between widgets above
        ttk.Label(self.tab_main, background=self.right_panel_bg, width=15).grid(column=0, row=0)
        ttk.Label(self.tab_main, background=self.right_panel_bg, width=50).grid(column=4, row=0)
        # This plot will be for testing display only. More functions and classes will be made to support plotting.
        self.tab_main_plot()
        self.plot_toolbar_main()


        self.header_number_clusters = ttk.Label(
            self.tab_spikesorting,
            text='Number of Clusters:',
            style='right_headers.TLabel',
            background=self.right_panel_bg
        )
        self.header_number_clusters.grid(column=1, row=0, sticky='w', padx=self.header_align_pad, pady=self.right_header_pady)

        # self.entry_number_clusters = ttk.Entry(
        #     self.tab_spikesorting,
        #     justify='left',
        #     width=20,
        #     state='disabled',
        #     style='unchecked.TEntry'
        # )
        # self.entry_number_clusters.grid(column=1, row=1, sticky='w', padx=self.widget_align_pad)
        #
        # self.button_save_number_clusters = ttk.Button(
        #     self.tab_spikesorting,
        #     text='Set',
        #     command=self.button_set_number_clusters_press,
        #     style='TButton'
        # )
        # self.button_save_number_clusters.grid(column=2, row=1, sticky='w')
        #
        # self.label_number_clusters = ttk.Label(
        #     self.tab_spikesorting,
        #     text='Clusters: 0',
        #     style='outputs.TLabel',
        #     background=self.right_panel_bg
        # )
        # self.label_number_clusters.grid(column=3, row=1, sticky='w', padx=self.widget_align_pad)

        self.dropdown_selected_cluster = StringVar()
        self.desired_clusters_menu = None
        self.dropdown_desired_cluster = ttk.Menubutton(
            master=self.tab_spikesorting,
            direction='below',
            text='Number of Clusters',
            menu=self.desired_clusters_menu,
            style='clusters.TMenubutton',
            width=20
        )
        self.dropdown_desired_cluster.grid(column=1, row=2, columnspan=1, sticky='we', padx=self.widget_align_pad,
                                           pady=self.right_header_pady)

        # self.entry_desired_clusters = ttk.Entry(
        #     self.tab_spikesorting,
        #     justify='left',
        #     width=20,
        #     state='disabled',
        #     style='unchecked.TEntry'
        # )
        # self.entry_desired_clusters.grid(column=1, row=3, sticky='w', padx=self.widget_align_pad)
        #
        # self.button_save_desired_clusters = ttk.Button(
        #     self.tab_spikesorting,
        #     text='Set',
        #     style='TButton',
        #     command=self.button_set_desired_clusters_press
        # )
        # self.button_save_desired_clusters.grid(column=2, row=3, sticky='w')
        #
        # self.label_desired_clusters = ttk.Label(
        #     self.tab_spikesorting,
        #     text='Desired: 0',
        #     style='outputs.TLabel',
        #     background=self.right_panel_bg
        # )
        # self.label_desired_clusters.grid(column=3, row=3, sticky='w', padx=self.widget_align_pad)

        self.tab_spikesorting_plot()

        ttk.Label(self.tab_spikesorting, background=self.right_panel_bg, width=15).grid(column=0, row=0)
        ttk.Label(self.tab_spikesorting, background=self.right_panel_bg, width=100).grid(column=4, row=0)

        self.plot_toolbar_spikesorting()

        # self.header_desired_clusters = ttk.Label(
        #     self.tab_spikesorting,
        #     text='Desired Clusters:',
        #     style='right_headers.TLabel',
        #     background=self.right_panel_bg
        # )
        # self.header_desired_clusters.grid(column=1, row=2, sticky='w', padx=self.header_align_pad,
        #                                   pady=self.right_header_pady)


        self.header_psd_fFrom = ttk.Label(
            self.tab_psd,
            text='fFrom:',
            style='right_headers.TLabel',
            background=self.right_panel_bg
        )
        self.header_psd_fFrom.grid(column=1, row=0, sticky='w', padx=self.header_align_pad, pady=self.right_header_pady)

        self.entry_psd_fFrom = ttk.Entry(
            self.tab_psd,
            justify='left',
            width=20,
            style='unchecked.TEntry',
            state='disabled'
        )
        self.entry_psd_fFrom.grid(column=1, row=1, sticky='w', padx=self.widget_align_pad)

        self.header_psd_fTo = ttk.Label(
            self.tab_psd,
            text='fTo:',
            style='right_headers.TLabel',
            background=self.right_panel_bg
        )
        self.header_psd_fTo.grid(column=2, row=0, sticky='w', pady=self.right_header_pady)

        self.entry_psd_fTo = ttk.Entry(
            self.tab_psd,
            justify='left',
            width=20,
            style='unchecked.TEntry',
            state='disabled'
        )
        self.entry_psd_fTo.grid(column=2, row=1, sticky='w')

        self.button_save_psd_fFrom_fTo = ttk.Button(
            self.tab_psd,
            text='Set',
            command=self.button_set_psd_fFromfTo_press,
            style='TButton'
        )
        self.button_save_psd_fFrom_fTo.grid(column=3, row=1, sticky='w', padx=self.widget_align_pad)

        self.label_psd_fFrom_fTo = ttk.Label(
            self.tab_psd,
            text='fFrom: 0 s\nfTo: 0 s',
            style='outputs.TLabel',
            background=self.right_panel_bg
        )
        self.label_psd_fFrom_fTo.grid(column=4, row=1, sticky='w')

        ttk.Label(self.tab_psd, background=self.right_panel_bg, width=15).grid(column=0, row=0)
        ttk.Label(self.tab_psd, background=self.right_panel_bg, width=50).grid(column=4, row=0)

        self.tab_psd_plot()
        self.plot_toolbar_psd()


        # self.header_lag_time = ttk.Label(
        #     self.tab_oscillations,
        #     text='Lag Time:',
        #     style='right_headers.TLabel',
        #     background=self.right_panel_bg
        # )
        # self.header_lag_time.grid(column=1, row=0, sticky='w', padx=self.header_align_pad, pady=self.right_header_pady)
        #
        # self.entry_lag_time = ttk.Entry(
        #     self.tab_oscillations,
        #     justify='left',
        #     width=20,
        #     style='unchecked.TEntry',
        #     state='disabled'
        # )
        # self.entry_lag_time.grid(column=1, row=1, sticky='w', padx=self.widget_align_pad)
        #
        # self.header_time_interval = ttk.Label(
        #     self.tab_oscillations,
        #     text='Time Interval:',
        #     style='right_headers.TLabel',
        #     background=self.right_panel_bg
        # )
        # self.header_time_interval.grid(column=2, row=0, sticky='w',
        #                                pady=self.right_header_pady)
        #
        # self.entry_time_interval = ttk.Entry(
        #     self.tab_oscillations,
        #     justify='left',
        #     width=20,
        #     style='unchecked.TEntry',
        #     state='disabled'
        # )
        # self.entry_time_interval.grid(column=2, row=1, sticky='w')
        #
        # self.button_save_lag_time_interval = ttk.Button(
        #     self.tab_oscillations,
        #     text='Set',
        #     command=self.button_set_lag_time_interval_press,
        #     style='TButton'
        # )
        # self.button_save_lag_time_interval.grid(column=3, row=1, sticky='w', padx=self.widget_align_pad)
        #
        # self.label_lag_time_interval = ttk.Label(
        #     self.tab_oscillations,
        #     text='Lag: 0.5 s\nInterval: 0.01 s',
        #     style='outputs.TLabel',
        #     background=self.right_panel_bg
        # )
        # self.label_lag_time_interval.grid(column=4, row=1, sticky='w')

        # ttk.Label(self.tab_oscillations, background=self.right_panel_bg, width=15).grid(column=0, row=0)
        # ttk.Label(self.tab_oscillations, background=self.right_panel_bg, width=50).grid(column=4, row=0)

        self.tab_oscillations_plot()
        self.plot_toolbar_oscillations()





        self.event_loop()

    def event_loop(self):
        while True:
            self.update()
            try:
                self.state()
            except:
                sys.exit()


    def dynamic_resize(self):
        for i in range(self.num_columns):
            Grid.columnconfigure(self, index=i, weight=1)
        for i in range(self.num_rows):
            Grid.rowconfigure(self, index=i, weight=1)


    def button_file_press(self):
        supported_filetypes = [
            ('Spike2 Datafile 32-bit (*.smr)', '*.smr'),
            ('Spike2 Datafile 64-bit (*.smrx)', '*.smrx')
        ]
        filepath = tkf.askopenfilename(
            title='Open File',
            filetypes=supported_filetypes
        )
        filename = filepath.split('/')[-1]
        file_ext = filepath.split('.')[-1]
        self.label_file['text'] = filename

        if filepath == '' or file_ext != 'smr':
            self.reset_all_parameters()
            read.channel = 0
            read.read_filepath = ''

            if len(self.channels_available) != 0:
                for i in range(len(read.get_channel_list())):
                    self.dropdown_menu.delete(0)
                self.dropdown_channel['text'] = 'Channels'

            self.entry_threshold_factor['style'] = 'unchecked.TEntry'
            self.entry_threshold_factor['state'] = 'disabled'

            if filepath == '':
                tkm.showerror(
                    title='File Error',
                    message='No file selected.'
                )
            elif file_ext != 'smr':
                tkm.showerror(
                    title='File Error',
                    message=f'Invalid file extension (*.{file_ext}).'
                )
        else:
            self.reset_all_parameters()

            read.read_filepath = filepath

            self.channels_available = read.get_channel_list()
            self.dropdown_menu = Menu(tearoff=False)
            for i in self.channels_available:
                self.dropdown_menu.add_radiobutton(
                    label=f'{i}',
                    variable=self.dropdown_selection,
                    command=self.dropdown_channel_choice,
                    font=self.menu_item_font
                )

            read.get_lowest_wave_channel()

            self.dropdown_channel['menu'] = self.dropdown_menu
            self.dropdown_channel['text'] = f'{self.channels_available[read.channel]}'
            raw_wave = read.read_smr_waves()
            spike.main_times = raw_wave[0]
            spike.main_magnitudes = raw_wave[1]
            spike.main_fs = raw_wave[2]

            self.delete_plot(canvas=self.plot_canvas_main)
            self.tab_main_plot()
            self.plot_toolbar_main()

            self.entry_threshold_factor['state'] = 'normal'
            self.entry_threshold_factor['style'] = 'TEntry'
            self.entry_threshold_factor.delete(0, END)
            self.entry_threshold_factor.insert(0, '2')

            spike.psd_good_file = True
            frequencies, psd = spike.get_psd(magnitudes=spike.main_magnitudes, fs=spike.main_fs)
            spike.psd_frequencies = frequencies
            spike.psd_power = psd
            self.entry_psd_fFrom['state'] = 'normal'
            self.entry_psd_fFrom['style'] = 'TEntry'
            self.entry_psd_fFrom.delete(0, END)
            self.entry_psd_fFrom.insert(0, '0')
            self.entry_psd_fTo['state'] = 'normal'
            self.entry_psd_fTo['style'] = 'TEntry'
            self.entry_psd_fTo.delete(0, END)
            self.entry_psd_fTo.insert(0, '100')
            self.label_psd_fFrom_fTo['text'] = 'fFrom: 0 Hz\nfTo: 100 Hz'

            spike.psd_plot_xlim = [0, 100]
            self.delete_plot(self.plot_canvas_psd)
            self.tab_psd_plot()
            self.plot_toolbar_psd()


    def dropdown_channel_choice(self):
        print(f'Channel: {self.dropdown_selection.get()}')
        selected_channel = self.dropdown_selection.get()
        self.dropdown_channel['text'] = selected_channel

        read.channel = int(selected_channel.split(' ')[0]) - 1

        self.reset_all_parameters()

        raw_wave = read.read_smr_waves()
        if raw_wave == -1:
            print('Error: the selected channel does not contain any wave data.')
            tkm.showerror(
                title='Type Error',
                message='The channel you selected contains no wave data.'
            )
        else:
            spike.main_times = raw_wave[0]
            spike.main_magnitudes = raw_wave[1]
            spike.main_fs = raw_wave[2]

            self.delete_plot(self.plot_canvas_main)
            self.tab_main_plot()
            self.plot_toolbar_main()

            spike.psd_good_file = True
            frequencies, psd = spike.get_psd(magnitudes=spike.main_magnitudes, fs=spike.main_fs)
            spike.psd_frequencies = frequencies
            spike.psd_power = psd

            spike.psd_plot_xlim = [0, 100]
            self.delete_plot(self.plot_canvas_psd)
            self.tab_psd_plot()
            self.plot_toolbar_psd()


    def check_filtering_choice(self):
        print(f'Filtering: {self.check_filtering_status.get()}')    # IntVar object has .get() method
        if self.check_filtering_status.get() == 1:
            self.entry_highpass_cutoff['state'] = 'normal'
            self.entry_lowpass_cutoff['state'] = 'normal'
            self.entry_highpass_cutoff['style'] = 'TEntry'
            self.entry_lowpass_cutoff['style'] = 'TEntry'

            read.filtering = True
        else:
            spike.main_magnitudes = read.read_smr_waves()[1]
            threshold_unfilt = spike.get_main_threshold(magnitudes=spike.main_magnitudes,
                                                        factor=spike.main_threshold_factor)
            spike.main_threshold = threshold_unfilt
            self.label_threshold['text'] = f'Threshold: {round(threshold_unfilt, 2)} V ' \
                                           f'[Factor: {self.entry_threshold_factor.get()}]'

            self.entry_lowpass_cutoff.delete(0, END)
            self.entry_lowpass_cutoff.insert(0, '')
            self.entry_highpass_cutoff.delete(0, END)
            self.entry_highpass_cutoff.insert(0, '')
            self.label_set_cutoffs['text'] = f'HP: 0 Hz\nLP: 0 Hz'

            self.entry_highpass_cutoff['state'] = 'disabled'
            self.entry_lowpass_cutoff['state'] = 'disabled'
            self.entry_highpass_cutoff['style'] = 'unchecked.TEntry'
            self.entry_lowpass_cutoff['style'] = 'unchecked.TEntry'

            read.filtering = False

            if spike.spikesorting_plot_disp:
                self.update_spikesorting_parameters()
                self.delete_plot(self.plot_canvas_spikesorting)
                self.tab_spikesorting_plot()
                self.plot_toolbar_spikesorting()

            self.delete_plot(canvas=self.plot_canvas_main)
            self.tab_main_plot()
            self.plot_toolbar_main()

            # Need conditional because if filtering unchecked and no threshold set, uses all points > 0 V
            # in the loops and reduces speed.
            if spike.main_threshold_set and spike.lag_time != 0 and spike.time_interval != 0:
                self.update_oscillations_parameters()
                self.delete_plot(self.plot_canvas_oscillations)
                self.tab_oscillations_plot()
                self.plot_toolbar_oscillations()


    def check_spikesorting_choice(self):
        print(f'Spike sorting: {self.check_spikesorting_status.get()}')
        if self.check_spikesorting_status.get() == 1 and spike.main_threshold_set:
            spike.spikesorting_plot_disp = True

            self.update_spikesorting_parameters()
            self.delete_plot(canvas=self.plot_canvas_spikesorting)
            self.tab_spikesorting_plot()
            self.plot_toolbar_spikesorting()

            self.dropdown_colour_cluster_menu = Menu(tearoff=False)
            colours = ['Red', 'Blue', 'Orange', 'Cyan', 'Purple', 'Brown']
            for i in colours:
                self.dropdown_colour_cluster_menu.add_radiobutton(
                    label=f'{i}',
                    variable=self.dropdown_selected_colour_cluster,
                    command=self.dropdown_colour_cluster_choice,
                    font=self.menu_item_font
                )
            self.dropdown_colour_cluster['menu'] = self.dropdown_colour_cluster_menu

            self.desired_clusters_menu = Menu(tearoff=False)
            for i in range(2, 7):
                self.desired_clusters_menu.add_radiobutton(
                    label=f'{i}',
                    variable=self.dropdown_selected_cluster,
                    command=self.dropdown_desired_clusters_choice,
                    font=self.menu_item_font
                )
            self.dropdown_desired_cluster['menu'] = self.desired_clusters_menu

            read.number_clusters = spike.number_desired_clusters

            # Needs to be here because self.tab_spikesorting_plot() runs spike.spike_sorting() which is
            # essential for the main plot to update, since it provides the colour labels.
            self.delete_plot(self.plot_canvas_main)
            self.tab_main_plot()
            self.plot_toolbar_main()

            self.dropdown_colour_cluster['text'] = spike.selected_cluster.title()

            # Since coloured spikes are selected automatically from checking spike sorting option,
            # must update parameters according to selected cluster.
            self.update_coloured_spikesorting_parameters()
            # Oscillations tab no longer passing un-sorted peaks, but automatically-selected cluster peaks.
            self.update_coloured_oscillations_parameters()
            self.delete_plot(self.plot_canvas_oscillations)
            self.tab_oscillations_plot()
            self.plot_toolbar_oscillations()
            # Update left panel
            self.update_properties()

            if self.check_spiketrain_status.get() == 0:
                self.check_spiketrain.invoke()
            if self.check_eventpeaks_status.get() == 0:
                self.check_eventpeaks.invoke()
            if self.check_thresholdbar_status.get() == 0:
                self.check_thresholdbar.invoke()

        elif self.check_spikesorting_status.get() == 1:
            print('Error: no threshold has been set. Spike-sorting is unavailable.')
        else:
            for i in range(2, 7):
                self.desired_clusters_menu.delete(0)
            self.dropdown_desired_cluster['text'] = 'Number of Clusters'

            for i in range(6):
                self.dropdown_colour_cluster_menu.delete(0)
            self.dropdown_colour_cluster['text'] = 'Cluster'

            spike.spikesorting_plot_disp = False
            spike.selected_cluster = None
            spike.selected_cluster_peaks = []

            self.delete_plot(canvas=self.plot_canvas_spikesorting)
            self.tab_spikesorting_plot()
            self.plot_toolbar_spikesorting()

            self.delete_plot(self.plot_canvas_main)
            self.tab_main_plot()
            self.plot_toolbar_main()

            self.update_properties()

            self.label_silhouette_output['text'] = ''


    def check_invertsegment_choice(self):
        if self.check_invertsegment_status.get() == 1:
            spike.segment_inverted = True
            read.segment_inverted = True
            # spike.main_magnitudes *= -1
            spike.main_threshold *= - 1

            self.update_spikesorting_parameters()
            self.delete_plot(self.plot_canvas_spikesorting)
            self.tab_spikesorting_plot()
            self.plot_toolbar_spikesorting()

            self.delete_plot(self.plot_canvas_main)
            self.tab_main_plot()
            self.plot_toolbar_main()

            self.update_oscillations_parameters()
            self.delete_plot(self.plot_canvas_oscillations)
            self.tab_oscillations_plot()
            self.plot_toolbar_oscillations()
        else:
            spike.segment_inverted = False
            read.segment_inverted = False
            # spike.main_magnitudes *= -1
            spike.main_threshold *= -1

            self.update_spikesorting_parameters()
            self.delete_plot(self.plot_canvas_spikesorting)
            self.tab_spikesorting_plot()
            self.plot_toolbar_spikesorting()

            self.delete_plot(self.plot_canvas_main)
            self.tab_main_plot()
            self.plot_toolbar_main()

            self.update_oscillations_parameters()
            self.delete_plot(self.plot_canvas_oscillations)
            self.tab_oscillations_plot()
            self.plot_toolbar_oscillations()


    def ask_save_filepath(self, colour=None, properties=False):
        if colour == None and properties == True:
            ask_save = tkm.askyesno(
                title='Save Properties',
                message=f'Edit this message. Threshold: '
                        f'{round(spike.main_threshold, 2)} V.\n\n'
                        f'Would you like to save?'
            )
        elif colour == None:
            ask_save = tkm.askyesno(
                title='Save Segment',
                message=f'This segment contains wave data and a spike-train obtained using the set threshold of '
                        f'{round(spike.main_threshold, 2)} V.\n\n'
                        f'Would you like to save?'
            )
        else:
            ask_save = tkm.askyesno(
                title=f'Save Segment ({colour.title()})',
                message=f'This segment contains wave data and a spike-train obtained using the set threshold of '
                        f'{round(spike.main_threshold, 2)} V.\n\n'
                        f'Would you like to save?'
            )

        if ask_save:
            if properties == True:
                supported_filetypes = [
                    ('Comma Separated Values (*.csv)', '*.smr')
                ]
            else:
                supported_filetypes = [
                    ('Spike2 Datafile 32-bit (*.smr)', '*.smr')
                ]
            filetype_choice = StringVar()
            save_filepath = tkf.asksaveasfilename(
                title='Save As...',
                filetypes=supported_filetypes,
                typevariable=filetype_choice
            )
            if save_filepath.endswith('.smr') and filetype_choice.get() == 'Spike2 Datafile 32-bit (*.smr)':
                pass
            elif not save_filepath.endswith('.smr') and filetype_choice.get() == 'Spike2 Datafile 32-bit (*.smr)':
                save_filepath += '.smr'
            # elif save_filepath.endswith('.smrx') and filetype_choice.get() == 'Spike2 Datafile 64-bit (*.smrx)':
            #     pass
            # elif not save_filepath.endswith('.smrx') and filetype_choice.get() == 'Spike2 Datafile 64-bit (*.smrx)':
            #     save_filepath += '.smrx'

            if save_filepath.endswith('.csv') and filetype_choice.get() == 'Comma Separated Values (*.csv)':
                pass
            elif not save_filepath.endswith('.csv') and filetype_choice.get() == 'Comma Separated Values (*.csv)':
                save_filepath += '.csv'

            return save_filepath


    def dropdown_save_properties(self):
        if read.read_filepath == '':
            tkm.showerror(
                title='Data Error',
                message='There is nothing to save since no file was imported.'
            )
        else:
            read.save_filepath = self.ask_save_filepath(colour=None, properties=True)
            if read.save_filepath == '':
                tkm.showerror(
                    title='File Error',
                    message='No path selected.'
                )
            else:
                read.save_properties()


    def dropdown_save_segment_all_spikes(self):
        if len(spike.main_event_peak_times) == 0:
            tkm.showerror(
                title='Data Error',
                message=f'There are no spikes to save. Please check to see if a threshold was set.'
            )
            save = -1
        else:
            read.save_filepath = self.ask_save_filepath()
            if read.save_filepath == '':
                tkm.showerror(
                    title='File Error',
                    message='No path selected.'
                )
                save = -1
            else:
                save = read.save_segment(events=spike.main_event_peak_times)
        print(f'Save status code: {save}')


    def dropdown_save_segment(self):
        if self.dropdown_selected_colour_cluster.get() == 'Red' or spike.selected_cluster == 'red':
            spike.selected_cluster = 'red'
            spike.selected_cluster_peaks = spike.reds
        elif self.dropdown_selected_colour_cluster.get() == 'Blue' or spike.selected_cluster == 'blue':
            spike.selected_cluster = 'blue'
            spike.selected_cluster_peaks = spike.blues
        elif self.dropdown_selected_colour_cluster.get() == 'Orange' or spike.selected_cluster == 'orange':
            spike.selected_cluster = 'orange'
            spike.selected_cluster_peaks = spike.oranges
        elif self.dropdown_selected_colour_cluster.get() == 'Cyan' or spike.selected_cluster == 'cyan':
            spike.selected_cluster = 'cyan'
            spike.selected_cluster_peaks = spike.cyans
        elif self.dropdown_selected_colour_cluster.get() == 'Purple' or spike.selected_cluster == 'purple':
            spike.selected_cluster = 'purple'
            spike.selected_cluster_peaks = spike.purples
        elif self.dropdown_selected_colour_cluster.get() == 'Brown' or spike.selected_cluster == 'brown':
            spike.selected_cluster = 'brown'
            spike.selected_cluster_peaks = spike.browns
        elif spike.selected_cluster == None:
            spike.selected_cluster_peaks = []

        if len(spike.selected_cluster_peaks) == 0:
            if spike.selected_cluster == None:
                tkm.showerror(
                    title='Error',
                    message=f'There are no clusters to save. Spike sorting was not activated.'
                )
                save = -1
            else:
                tkm.showerror(
                    title='Data Error',
                    message=f'There are no {spike.selected_cluster.lower()} spikes to save.'
                )
                save = -1
        else:
            read.save_filepath = self.ask_save_filepath(colour=f'{spike.selected_cluster.lower()}')
            if read.save_filepath == '':
                tkm.showerror(
                    title='File Error',
                    message='No path selected.'
                )
                save = -1
            else:
                save = read.save_segment(events=spike.selected_cluster_peaks)
        print(f'Save status code: {save}')


    def button_set_cutoffs_press(self):
        print(f'Highpass cut-off: {self.entry_highpass_cutoff.get()}')
        print(f'Lowpass cut-off: {self.entry_lowpass_cutoff.get()}')
        mincutoff = self.entry_highpass_cutoff.get()
        maxcutoff = self.entry_lowpass_cutoff.get()
        self.label_set_cutoffs['text'] = f'HP: {round(float(mincutoff), 2)} Hz\n' \
                                         f'LP: {round(float(maxcutoff), 2)} Hz'
        filt_wave = spike.filter_wave(wave=spike.main_magnitudes,
                                      mincutoff=mincutoff, maxcutoff=maxcutoff, fs=spike.main_fs)
        spike.main_magnitudes = filt_wave

        threshold = spike.get_main_threshold(magnitudes=spike.main_magnitudes,
                                             factor=spike.main_threshold_factor)
        spike.main_threshold = threshold
        read.threshold = threshold
        read.highpass = float(mincutoff)
        read.lowpass = float(maxcutoff)

        if spike.num_threshold_sets == 0:
            self.entry_threshold_factor.delete(0, END)
            self.label_threshold['text'] = f'Threshold: {round(threshold, 2)} V ' \
                                       f'[Factor: 0]'
            self.entry_threshold_factor.insert(0, '2')
        else:
            self.label_threshold['text'] = f'Threshold: {round(threshold, 2)} V ' \
                                           f'[Factor: {self.entry_threshold_factor.get()}]'

        if spike.spikesorting_plot_disp:
            self.update_spikesorting_parameters()
            self.delete_plot(self.plot_canvas_spikesorting)
            self.tab_spikesorting_plot()
            self.plot_toolbar_spikesorting()

        self.delete_plot(self.plot_canvas_main)
        self.tab_main_plot()
        self.plot_toolbar_main()

        if spike.main_threshold != 0:
            self.update_spikesorting_parameters()
            self.delete_plot(self.plot_canvas_spikesorting)
            self.tab_spikesorting_plot()
            self.plot_toolbar_spikesorting()

        if spike.main_threshold != 0 and spike.lag_time != 0 and spike.time_interval != 0:
            self.update_oscillations_parameters()
            self.delete_plot(self.plot_canvas_oscillations)
            self.tab_oscillations_plot()
            self.plot_toolbar_oscillations()


    def button_set_threshold_press(self):
        if self.entry_threshold_factor.get() == '' or self.entry_threshold_factor.get() == 0:
            spike.main_threshold = 0
            spike.main_threshold_factor = 0
            spike.main_disp_threshold = False
            spike.main_threshold_set = False
            spike.spikesorting_plot_disp = False    # Requires valid threshold.
            self.label_threshold['text'] = f'Threshold: {round(spike.main_threshold, 2)} V ' \
                                           f'[Factor: {spike.main_threshold_factor}]'

            if self.check_spikesorting_status.get() == 1:
                # Spike sorting clear plot bc no threshold
                self.check_spikesorting.invoke()
            else:
                # Spike sorting clear plot bc no threshold
                # self.delete_plot(canvas=self.plot_canvas_spikesorting)
                # self.tab_spikesorting_plot()
                # self.plot_toolbar_spikesorting()

                # Main plot reset if spike sorting was enabled.
                self.delete_plot(self.plot_canvas_main)
                self.tab_main_plot()
                self.plot_toolbar_main()

            # Oscillations clear plot bc no threshold
            self.delete_plot(canvas=self.plot_canvas_oscillations)
            self.tab_oscillations_plot()
            self.plot_toolbar_oscillations()

            # self.entry_lag_time.delete(0, END)
            # self.entry_lag_time.insert(0, '')
            # self.entry_lag_time['state'] = 'disabled'
            # self.entry_lag_time['style'] = 'unchecked.TEntry'
            # self.entry_time_interval.delete(0, END)
            # self.entry_time_interval.insert(0, '')
            # self.entry_time_interval['state'] = 'disabled'
            # self.entry_time_interval['style'] = 'unchecked.TEntry'
            # self.label_lag_time_interval['text'] = 'Lag: 0.5 s\nInterval: 0.01 s'

            self.label_snr_output['text'] = ''
            self.label_interspike_output['text'] = ''
            self.label_firing_rate_output['text'] = ''
            self.label_burst_index_output['text'] = ''
            self.label_variation_coefficient_output['text'] = ''
            self.label_theta_output['text'] = ''
            self.label_alpha_output['text'] = ''
            self.label_low_beta_output['text'] = ''
            self.label_high_beta_output['text'] = ''

            read.threshold = None
            read.threshold_factor = None

        else:
            spike.num_threshold_sets += 1
            spike.main_threshold_set = True
            spike.main_threshold_factor = float(self.entry_threshold_factor.get())
            spike.main_threshold = spike.get_main_threshold(magnitudes=spike.main_magnitudes,
                                                 factor=spike.main_threshold_factor)
            self.label_threshold['text'] = f'Threshold: {round(spike.main_threshold, 2)} V ' \
                                           f'[Factor: {spike.main_threshold_factor}]'

            events = spike.get_event_peaks(
                times=spike.main_times,
                magnitudes=spike.main_magnitudes,
                threshold=spike.main_threshold
            )
            spike.main_event_peak_times = events[0]
            spike.main_event_peak_mags = events[1]

            # self.entry_lag_time['state'] = 'normal'
            # self.entry_lag_time['style'] = 'TEntry'
            # self.entry_time_interval['state'] = 'normal'
            # self.entry_time_interval['style'] = 'TEntry'

            # Spike sorting update with newly set threshold.
            self.update_spikesorting_parameters()
            self.delete_plot(canvas=self.plot_canvas_spikesorting)
            self.tab_spikesorting_plot()
            self.plot_toolbar_spikesorting()

            # Oscillation plot update with changing threshold.
            self.update_oscillations_parameters()
            self.delete_plot(self.plot_canvas_oscillations)
            self.tab_oscillations_plot()
            self.plot_toolbar_oscillations()

            self.update_properties()

            read.threshold = spike.main_threshold
            read.threshold_factor = spike.main_threshold_factor


        self.delete_plot(self.plot_canvas_main)
        self.tab_main_plot()
        self.plot_toolbar_main()


    def check_spiketrain_choice(self):
        print(f'Spike Train: {self.check_spiketrain_status.get()}')
        if self.check_spiketrain_status.get() == 1:
            spike.main_spiketrain_disp = True
        else:
            spike.main_spiketrain_disp = False
        self.delete_plot(canvas=self.plot_canvas_main)
        self.tab_main_plot()
        self.plot_toolbar_main()


    def check_thresholdbar_choice(self):
        print(f'Threshold Bar: {self.check_thresholdbar_status.get()}')
        if self.check_thresholdbar_status.get() == 1:
            spike.main_disp_threshold = True
        else:
            spike.main_disp_threshold = False
        self.delete_plot(canvas=self.plot_canvas_main)
        self.tab_main_plot()
        self.plot_toolbar_main()


    def check_eventpeaks_choice(self):
        if spike.main_threshold == 0:
            spike.main_event_peaks_disp = False
        elif self.check_eventpeaks_status.get() == 1:
            spike.main_event_peaks_disp = True

        self.delete_plot(canvas=self.plot_canvas_main)
        self.tab_main_plot()
        self.plot_toolbar_main()


    def dropdown_colour_cluster_choice(self):
        self.dropdown_colour_cluster['text'] = self.dropdown_selected_colour_cluster.get()
        spike.selected_cluster = self.dropdown_selected_colour_cluster.get().lower()

        # self.update_spikesorting_parameters()
        # self.delete_plot(self.plot_canvas_spikesorting)
        # self.tab_spikesorting_plot()
        # self.plot_toolbar_spikesorting()

        self.delete_plot(self.plot_canvas_main)
        self.tab_main_plot()
        self.plot_toolbar_main()

        self.update_coloured_spikesorting_parameters()

        self.update_coloured_oscillations_parameters()
        self.delete_plot(self.plot_canvas_oscillations)
        self.tab_oscillations_plot()
        self.plot_toolbar_oscillations()

        self.update_properties()


    def tab_main_plot(self):
        fig = spike.main_plot()
        self.plot_canvas_main = FigureCanvasTkAgg(fig, master=self.tab_main)
        self.plot_canvas_main.get_tk_widget().grid(column=0, row=8, columnspan=5, sticky='nsew')
        self.plot_canvas_main.draw()


    def plot_toolbar_main(self):
        self.toolbar_frame = Frame(self.tab_main)
        self.toolbar_frame.grid(column=0, row=9, columnspan=5)
        self.toolbar = NavigationToolbar2Tk(self.plot_canvas_main, self.toolbar_frame)
        self.toolbar.configure(background='white')
        self.toolbar._message_label.configure(background='white')
        for i in self.toolbar.winfo_children():
            i.configure(background='white', bd=0)
        # self.toolbar.deletecommand(f'.!notebook.!frame{frame_num}.!frame.!navigationtoolbar2tk.!button4')
        # self.toolbar.deletecommand(f'.!notebook.!frame{frame_num}.!frame.!navigationtoolbar2tk.!label')
        # self.toolbar.deletecommand(f'.!notebook.!frame{frame_num}.!frame.!navigationtoolbar2tk.!label2')


    def delete_plot(self, canvas):
        for i in canvas.get_tk_widget().find_all():
            canvas.get_tk_widget().delete(i)


    # def button_set_number_clusters_press(self):
    #     self.label_number_clusters['text'] = f'Clusters: {self.entry_number_clusters.get()}'


    def dropdown_desired_clusters_choice(self):
        self.dropdown_desired_cluster['text'] = f'{self.dropdown_selected_cluster.get()}'
        spike.number_desired_clusters = int(self.dropdown_selected_cluster.get())
        read.number_clusters = spike.number_desired_clusters

        self.update_spikesorting_parameters()
        self.delete_plot(self.plot_canvas_spikesorting)
        self.tab_spikesorting_plot()
        self.plot_toolbar_spikesorting()

        self.delete_plot(self.plot_canvas_main)
        self.tab_main_plot()
        self.plot_toolbar_main()

        self.update_coloured_spikesorting_parameters()
        self.update_coloured_oscillations_parameters()
        self.delete_plot(self.plot_canvas_oscillations)
        self.tab_oscillations_plot()
        self.plot_toolbar_oscillations()

        self.update_properties()


    def tab_spikesorting_plot(self):
        fig = spike.spikesorting_plot()
        self.plot_canvas_spikesorting = FigureCanvasTkAgg(fig, master=self.tab_spikesorting)
        self.plot_canvas_spikesorting.get_tk_widget().grid(column=0, row=3, columnspan=5, sticky='nsew')
        self.plot_canvas_spikesorting.draw()


    def plot_toolbar_spikesorting(self):
        self.toolbar_frame = Frame(self.tab_spikesorting)
        self.toolbar_frame.grid(column=0, row=4, columnspan=5)
        self.toolbar = NavigationToolbar2Tk(self.plot_canvas_spikesorting, self.toolbar_frame)
        self.toolbar.configure(background='white')
        self.toolbar._message_label.configure(background='white')
        for i in self.toolbar.winfo_children():
            i.configure(background='white', bd=0)


    def button_set_psd_fFromfTo_press(self):
        self.label_psd_fFrom_fTo['text'] = f'fFrom: {self.entry_psd_fFrom.get()} Hz\n' \
                                           f'fTo: {self.entry_psd_fTo.get()} Hz'
        spike.psd_plot_xlim = [int(self.entry_psd_fFrom.get()), int(self.entry_psd_fTo.get())]
        self.delete_plot(self.plot_canvas_psd)
        self.tab_psd_plot()
        self.plot_toolbar_psd()


    def tab_psd_plot(self):
        fig = spike.psd_plot()
        self.plot_canvas_psd = FigureCanvasTkAgg(fig, master=self.tab_psd)
        self.plot_canvas_psd.get_tk_widget().grid(column=0, row=2, columnspan=5, sticky='nesw')
        self.plot_canvas_psd.draw()


    def plot_toolbar_psd(self):
        self.toolbar_frame = Frame(self.tab_psd)
        self.toolbar_frame.grid(column=0, row=3, columnspan=5)
        self.toolbar = NavigationToolbar2Tk(self.plot_canvas_psd, self.toolbar_frame)
        self.toolbar.configure(background='white')
        self.toolbar._message_label.configure(background='white')
        for i in self.toolbar.winfo_children():
            i.configure(background='white', bd=0)


    # def button_set_lag_time_interval_press(self):
    #     self.label_lag_time_interval['text'] = f'Lag: {self.entry_lag_time.get()} s\n' \
    #                                            f'Interval: {self.entry_time_interval.get()} s'
    #     if spike.main_threshold_set:
    #         spike.lag_time = float(self.entry_lag_time.get())
    #         read.lag_time = spike.lag_time
    #         spike.time_interval = float(self.entry_time_interval.get())
    #         read.time_interval = spike.time_interval
    #
    #         self.update_oscillations_parameters()
    #         self.delete_plot(self.plot_canvas_oscillations)
    #         self.tab_oscillations_plot()
    #         self.plot_toolbar_oscillations()
    #     else:
    #         print('Error: no threshold has been set. Oscillation analysis unavailable.')


    def tab_oscillations_plot(self):
        fig = spike.oscillations_plot()
        self.plot_canvas_oscillations = FigureCanvasTkAgg(fig, master=self.tab_oscillations)
        self.plot_canvas_oscillations.get_tk_widget().grid(column=0, row=0, sticky='nesw')
        self.plot_canvas_oscillations.draw()


    def plot_toolbar_oscillations(self):
        self.toolbar_frame = Frame(self.tab_oscillations)
        self.toolbar_frame.grid(column=0, row=5, columnspan=5)
        self.toolbar = NavigationToolbar2Tk(self.plot_canvas_oscillations, self.toolbar_frame)
        self.toolbar.configure(background='white')
        self.toolbar._message_label.configure(background='white')
        for i in self.toolbar.winfo_children():
            i.configure(background='white', bd=0)


    def update_spikesorting_parameters(self):
        events = spike.get_event_peaks(magnitudes=spike.main_magnitudes,
                                       times=spike.main_times,
                                       threshold=spike.main_threshold)
        events, spikes = spike.get_spike_windows(
            magnitudes=spike.main_magnitudes,
            event_indices=events[2],
            spike_window=int(10 ** -3 * spike.main_fs))

        spike.peak_indices = events
        spike.spike_matrix = spikes


    def update_coloured_spikesorting_parameters(self):
        if spike.selected_cluster == 'red':
            coloured_events = spike.reds_indices
        elif spike.selected_cluster == 'blue':
            coloured_events = spike.blues_indices
        elif spike.selected_cluster == 'orange':
            coloured_events = spike.oranges_indices
        elif spike.selected_cluster == 'cyan':
            coloured_events = spike.cyans_indices
        elif spike.selected_cluster == 'purple':
            coloured_events = spike.purples_indices
        elif spike.selected_cluster == 'brown':
            coloured_events = spike.browns_indices

        events, spikes = spike.get_spike_windows(
            magnitudes=spike.main_magnitudes,
            event_indices=np.array(coloured_events),
            spike_window=int(10 ** -3 * spike.main_fs))

        spike.coloured_peak_indices = events
        spike.coloured_spike_matrix = spikes


    def update_oscillations_parameters(self):
        events = spike.get_event_peaks(magnitudes=spike.main_magnitudes,
                                       times=spike.main_times,
                                       threshold=spike.main_threshold)
        spike.frequencies_autocorr, spike.psd_autocorr, spike.binned_time, spike.autocorr = spike.get_spike_oscillations(
            magnitudes=spike.main_magnitudes,
            event_indices=events[2],
            fs=spike.main_fs,
            lag_time=spike.lag_time,
            time_interval=spike.time_interval
        )


    def update_coloured_oscillations_parameters(self):
        if spike.selected_cluster == 'red':
            coloured_events = spike.reds_indices
        elif spike.selected_cluster == 'blue':
            coloured_events = spike.blues_indices
        elif spike.selected_cluster == 'orange':
            coloured_events = spike.oranges_indices
        elif spike.selected_cluster == 'cyan':
            coloured_events = spike.cyans_indices
        elif spike.selected_cluster == 'purple':
            coloured_events = spike.purples_indices
        elif spike.selected_cluster == 'brown':
            coloured_events = spike.browns_indices

        spike.frequencies_autocorr, spike.psd_autocorr, spike.binned_time, spike.autocorr = spike.get_spike_oscillations(
            magnitudes=spike.main_magnitudes,
            event_indices=np.array(coloured_events),
            fs=spike.main_fs,
            lag_time=spike.lag_time,
            time_interval=spike.time_interval
        )


    def update_properties(self):
        snr = spike.get_snr()
        read.snr = snr
        self.label_snr_output['text'] = round(snr, 2)
        percent_isi_violations = spike.get_percent_isi_violations()
        read.percent_isi_violations = percent_isi_violations
        self.label_interspike_output['text'] = round(percent_isi_violations, 2)
        firing_rate = spike.get_firing_rate()
        read.firing_rate = firing_rate
        self.label_firing_rate_output['text'] = round(firing_rate, 2)
        burst_index = spike.get_burst_index()
        read.burst_index = burst_index
        self.label_burst_index_output['text'] = round(burst_index, 2)
        cov = spike.get_variation_coefficient()
        read.cov = cov
        self.label_variation_coefficient_output['text'] = round(cov, 2)
        if spike.spikesorting_plot_disp:
            silhouette = spike.silhouette
            read.silhouette = silhouette
            self.label_silhouette_output['text'] = round(silhouette, 2)

        theta_power, alpha_power, low_beta_power, high_beta_power = spike.get_wave_powers()
        self.label_theta_output['text'] = round(theta_power, 2)
        read.theta_power = theta_power
        self.label_alpha_output['text'] = round(alpha_power, 2)
        read.alpha_power = alpha_power
        self.label_low_beta_output['text'] = round(low_beta_power, 2)
        read.low_beta_power = low_beta_power
        self.label_high_beta_output['text'] = round(high_beta_power, 2)
        read.high_beta_power = high_beta_power


    def reset_all_parameters(self):
        if self.check_filtering_status.get() == 1:
            self.check_filtering.invoke()
        if self.check_spikesorting_status.get() == 1:
            self.check_spikesorting.invoke()
        if self.check_invertsegment_status.get() == 1:
            self.check_invertsegment.invoke()

        spike.main_magnitudes = []
        spike.main_times = []
        spike.main_events = []
        spike.main_event_peak_times = []
        spike.main_event_peak_mags = []
        spike.main_fs = 0
        spike.main_threshold_factor = 0
        spike.main_threshold = 0
        spike.main_threshold_set = False
        spike.main_disp_threshold = False
        spike.main_event_peaks_disp = False
        spike.main_spiketrain_disp = False
        spike.num_threshold_sets = 0

        spike.segment_inverted = False

        self.entry_highpass_cutoff.delete(0, END)
        self.entry_highpass_cutoff.insert(0, '')
        self.entry_lowpass_cutoff.delete(0, END)
        self.entry_lowpass_cutoff.insert(0, '')
        self.label_set_cutoffs['text'] = 'HP: 0 Hz\nLP: 0 Hz'
        self.entry_threshold_factor.delete(0, END)
        self.entry_threshold_factor.insert(0, '')
        self.label_threshold['text'] = 'Threshold: 0.0 V [Factor: 0]'

        self.delete_plot(self.plot_canvas_main)
        self.tab_main_plot()
        self.plot_toolbar_main()

        if self.check_spiketrain_status.get() == 1:
            self.check_spiketrain.invoke()
        if self.check_eventpeaks_status.get() == 1:
            self.check_eventpeaks.invoke()
        if self.check_thresholdbar_status.get() == 1:
            self.check_thresholdbar.invoke()

        spike.spikesorting_plot_disp = False
        spike.spike_matrix = [[]]
        spike.peak_indices = []
        spike.number_desired_clusters = 0
        spike.selected_cluster_peaks = []
        spike.selected_cluster = None

        self.delete_plot(self.plot_canvas_spikesorting)
        self.tab_spikesorting_plot()
        self.plot_toolbar_spikesorting()

        spike.psd_good_file = False
        self.psd_frequencies = []
        self.psd_power = []
        self.psd_plot_xlim = [0, 100]
        self.entry_psd_fFrom['style'] = 'unchecked.TEntry'
        self.entry_psd_fFrom.delete(0, END)
        self.entry_psd_fFrom.insert(0, '0')
        self.entry_psd_fFrom['state'] = 'disabled'
        self.entry_psd_fTo['style'] = 'unchecked.TEntry'
        self.entry_psd_fTo.delete(0, END)
        self.entry_psd_fTo.insert(0, '100')
        self.entry_psd_fTo['state'] = 'disabled'
        self.label_psd_fFrom_fTo['text'] = 'fFrom: 0 Hz\nfTo: 100 Hz'

        self.delete_plot(self.plot_canvas_psd)
        self.tab_psd_plot()
        self.plot_toolbar_psd()

        # spike.lag_time = 0.5
        # spike.time_interval = 0.01
        self.frequencies_autocorr = []
        self.psd_autocorr = []
        self.binned_time = []
        self.autocorr = []
        # self.entry_lag_time.delete(0, END)
        # self.entry_lag_time.insert(0, '')
        # self.entry_lag_time['state'] = 'disabled'
        # self.entry_lag_time['style'] = 'unchecked.TEntry'
        # self.entry_time_interval.delete(0, END)
        # self.entry_time_interval.insert(0, '')
        # self.entry_time_interval['state'] = 'disabled'
        # self.entry_time_interval['style'] = 'unchecked.TEntry'
        # self.label_lag_time_interval['text'] = f'Lag: 0.5 s\nInverval: 0.01 s'

        self.delete_plot(self.plot_canvas_oscillations)
        self.tab_oscillations_plot()
        self.plot_toolbar_oscillations()

        read.snr = 0
        read.percent_isi_violations = 0
        read.firing_rate = 0
        read.burst_index = 0
        read.cov = 0
        read.silhouette = 0
        read.delta_power = 0
        read.theta_power = 0
        read.alpha_power = 0
        read.low_beta_power = 0
        read.high_beta_power = 0
        read.gamma_power = 0

        read.threshold = None
        read.threshold_factor = None
        read.filtering = False
        read.highpass = None
        read.lowpass = None
        read.segment_inverted = False
        read.number_clusters = None
        # read.lag_time = 0.5
        # read.time_interval = 0.01

        self.label_snr_output['text'] = ''
        self.label_interspike_output['text'] = ''
        self.label_silhouette_output['text'] = ''
        self.label_firing_rate_output['text'] = ''
        self.label_burst_index_output['text'] = ''
        self.label_variation_coefficient_output['text'] = ''
        self.label_theta_output['text'] = ''
        self.label_alpha_output['text'] = ''
        self.label_low_beta_output['text'] = ''
        self.label_high_beta_output['text'] = ''
