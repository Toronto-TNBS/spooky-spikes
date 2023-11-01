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
from spikedata.AnalyzeMER import AnalyzeMER
from scipy.signal import welch

read = ReadSpike()
spike = SpikeAnalysis()
mer = AnalyzeMER()


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

        self.header_export = ttk.Label(
            text='Export:',
            style='headers.TLabel',
            background=self.left_panel_bg
        )
        self.header_export.grid(column=0, row=8, sticky='w', padx=self.header_align_pad)

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
        self.dropdown_save.grid(column=0, row=9, columnspan=2, sticky='we', padx=self.widget_align_pad, pady=5)

        self.num_rows = self.grid_size()[1]
        self.num_cols = self.grid_size()[0]
        self.dynamic_resize(self, self.num_rows, self.num_cols)

        self.notebook = ttk.Notebook(self)
        self.notebook.grid(column=2, row=0, rowspan=self.num_rows, sticky='nesw')
        self.tab_main = Frame(self.notebook, bg='white')
        self.tab_spikesorting = Frame(self.notebook, bg='white')
        self.tab_features = Frame(self.notebook, bg='white')
        self.tab_main.grid(column=2, row=0, sticky='nesw')
        self.tab_spikesorting.grid(column=2, row=0, sticky='nesw')
        self.tab_features.grid(column=2, row=0, sticky='nesw')
        self.notebook.add(self.tab_main, text='Main')
        self.notebook.add(self.tab_spikesorting, text='Spike Sorting')
        self.notebook.add(self.tab_features, text='Features')

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

        self.tab_spikesorting_plot()

        ttk.Label(self.tab_spikesorting, background=self.right_panel_bg, width=15).grid(column=0, row=0)
        ttk.Label(self.tab_spikesorting, background=self.right_panel_bg, width=100).grid(column=4, row=0)

        self.plot_toolbar_spikesorting()

        # --- Features Tab --- #

        self.label_patterned = ttk.Label(
            master=self.tab_features,
            text='Patterned',
            style='right_headers.TLabel',
            background='white'
        )
        self.label_patterned.grid(column=0, row=0, sticky='w', padx=self.header_align_pad, pady=self.features_headers_pady)

        self.frame_patterned = Frame(self.tab_features, bg=self.left_panel_bg)
        self.frame_patterned.grid(column=0, row=1, padx=self.frame_padx)

        self.label_firing_rate = ttk.Label(
            master=self.frame_patterned,
            text='Firing rate (Hz):',
            style='outputs.TLabel',
            background=self.left_panel_bg,
        )
        self.label_firing_rate.grid(column=0, row=0, sticky='we', padx=self.header_align_pad, pady=self.label_output_pady)

        self.label_firing_rate_output = ttk.Label(
            master=self.frame_patterned,
            text='',
            style='outputs.TLabel',
            background='white',
            width=self.label_output_width
        )
        self.label_firing_rate_output.grid(column=1, row=0, sticky='w', padx=self.header_align_pad, pady=self.label_output_pady)

        self.label_burst_index = ttk.Label(
            master=self.frame_patterned,
            text='Burst index:',
            style='outputs.TLabel',
            background=self.left_panel_bg,
        )
        self.label_burst_index.grid(column=0, row=1, sticky='we', padx=self.header_align_pad, pady=self.label_output_pady)

        self.label_burst_index_output = ttk.Label(
            master=self.frame_patterned,
            text='',
            style='outputs.TLabel',
            background='white',
            width=self.label_output_width
        )
        self.label_burst_index_output.grid(column=1, row=1, sticky='w', padx=self.header_align_pad, pady=self.label_output_pady)

        self.label_variation_coefficient = ttk.Label(
            master=self.frame_patterned,
            text='Coefficient of variation:',
            style='outputs.TLabel',
            background=self.left_panel_bg,
        )
        self.label_variation_coefficient.grid(column=0, row=2, sticky='we', padx=self.header_align_pad, pady=self.label_output_pady)

        self.label_variation_coefficient_output = ttk.Label(
            master=self.frame_patterned,
            text='',
            style='outputs.TLabel',
            background='white',
            width=self.label_output_width
        )
        self.label_variation_coefficient_output.grid(column=1, row=2, sticky='w', padx=self.header_align_pad, pady=self.label_output_pady)


        self.label_quality = ttk.Label(
            master=self.tab_features,
            text='Quality Metrics',
            style='right_headers.TLabel',
            background='white'
        )
        self.label_quality.grid(column=1, row=0, sticky='w', pady=self.features_headers_pady)

        self.frame_quality = Frame(self.tab_features, bg=self.left_panel_bg)
        self.frame_quality.grid(column=1, row=1)

        self.label_snr = ttk.Label(
            master=self.frame_quality,
            text='Signal to noise ratio:',
            style='outputs.TLabel',
            background=self.left_panel_bg,
        )
        self.label_snr.grid(column=0, row=0, sticky='we', padx=self.header_align_pad, pady=self.label_output_pady)

        self.label_snr_output = ttk.Label(
            master=self.frame_quality,
            text='',
            style='outputs.TLabel',
            background='white',
            width=self.label_output_width
        )
        self.label_snr_output.grid(column=1, row=0, sticky='w', padx=self.header_align_pad, pady=self.label_output_pady)

        self.label_interspike = ttk.Label(
            master=self.frame_quality,
            text='Percent ISI violations:',
            style='outputs.TLabel',
            background=self.left_panel_bg,
        )
        self.label_interspike.grid(column=0, row=1, sticky='we', padx=self.header_align_pad, pady=self.label_output_pady)

        self.label_interspike_output = ttk.Label(
            master=self.frame_quality,
            text='',
            style='outputs.TLabel',
            background='white',
            width=self.label_output_width
        )
        self.label_interspike_output.grid(column=1, row=1, sticky='w', padx=self.header_align_pad, pady=self.label_output_pady)

        self.label_silhouette = ttk.Label(
            master=self.frame_quality,
            text='Silhouette Score:',
            style='outputs.TLabel',
            background=self.left_panel_bg,
        )
        self.label_silhouette.grid(column=0, row=2, sticky='we', padx=self.header_align_pad, pady=self.label_output_pady)

        self.label_silhouette_output = ttk.Label(
            master=self.frame_quality,
            text='',
            style='outputs.TLabel',
            background='white',
            width=self.label_output_width
        )
        self.label_silhouette_output.grid(column=1, row=2, sticky='w', padx=self.header_align_pad, pady=self.label_output_pady)


        self.frame_patterned_aux = Frame(self.tab_features, bg=self.left_panel_bg)
        self.frame_patterned_aux.grid(column=2, row=1, sticky='wens', padx=self.frame_padx)

        self.label_plot_isi = ttk.Label(
            master=self.frame_patterned_aux,
            text='ISI Plot',
            style='outputs.TLabel',
            background=self.left_panel_bg,
            anchor='center'
        )
        self.label_plot_isi.grid(column=0, row=0, sticky='we', pady=self.right_header_pady)

        self.button_plot_isi = ttk.Button(
            master=self.frame_patterned_aux,
            text='Generate',
            style='TButton',
            command=self.button_plot_isi_press
        )
        self.button_plot_isi.grid(column=0, row=1, sticky='ew', padx=self.widget_align_pad, pady=self.right_header_pady)


        self.label_oscillations = ttk.Label(
            master=self.tab_features,
            text='Spiketrain Oscillations',
            style='right_headers.TLabel',
            background='white'
        )
        self.label_oscillations.grid(column=0, row=2, sticky='w', padx=self.header_align_pad, pady=self.features_headers_pady)

        self.frame_oscillations = Frame(self.tab_features, bg=self.left_panel_bg)
        self.frame_oscillations.grid(column=0, row=3, columnspan=2, sticky='nswe', padx=self.frame_padx)

        self.label_oscillations_power = ttk.Label(
            master=self.frame_oscillations,
            text='Power:',
            style='outputs.TLabel',
            background=self.left_panel_bg
        )
        self.label_oscillations_power.grid(column=0, row=1, sticky='we', padx=self.header_align_pad, pady=self.label_output_pady)

        self.label_oscillations_burst = ttk.Label(
            master=self.frame_oscillations,
            text='Burst Duration:',
            style='outputs.TLabel',
            background=self.left_panel_bg
        )
        self.label_oscillations_burst.grid(column=0, row=2, sticky='we', padx=self.header_align_pad, pady=self.label_output_pady)

        self.label_theta = ttk.Label(
            master=self.frame_oscillations,
            text='Theta',
            style='outputs.TLabel',
            background=self.left_panel_bg,
            anchor='center'
        )
        self.label_theta.grid(column=1, row=0, sticky='we', pady=self.label_output_pady)

        self.label_theta_power_output = ttk.Label(
            master=self.frame_oscillations,
            text='',
            style='outputs.TLabel',
            background='white',
            width=self.label_output_width
        )
        self.label_theta_power_output.grid(column=1, row=1, sticky='we')

        self.label_theta_burst_output = ttk.Label(
            master=self.frame_oscillations,
            text='',
            style='outputs.TLabel',
            background='white',
            width=self.label_output_width
        )
        self.label_theta_burst_output.grid(column=1, row=2, sticky='we', pady=self.label_output_pady)

        self.label_alpha = ttk.Label(
            master=self.frame_oscillations,
            text='Alpha',
            style='outputs.TLabel',
            background=self.left_panel_bg,
            anchor='center'
        )
        self.label_alpha.grid(column=2, row=0, sticky='we', padx=self.header_align_pad,
                              pady=self.label_output_pady)

        self.label_alpha_power_output = ttk.Label(
            master=self.frame_oscillations,
            text='',
            style='outputs.TLabel',
            background='white',
            width=self.label_output_width
        )
        self.label_alpha_power_output.grid(column=2, row=1, sticky='we', padx=self.header_align_pad)

        self.label_alpha_burst_output = ttk.Label(
            master=self.frame_oscillations,
            text='',
            style='outputs.TLabel',
            background='white',
            width=self.label_output_width
        )
        self.label_alpha_burst_output.grid(column=2, row=2, sticky='we', padx=self.header_align_pad,
                                           pady=self.label_output_pady)

        self.label_lowbeta = ttk.Label(
            master=self.frame_oscillations,
            text='Low-beta',
            style='outputs.TLabel',
            background=self.left_panel_bg,
            anchor='center'
        )
        self.label_lowbeta.grid(column=3, row=0, sticky='we', pady=self.label_output_pady)

        self.label_lowbeta_power_output = ttk.Label(
            master=self.frame_oscillations,
            text='',
            style='outputs.TLabel',
            background='white',
            width=self.label_output_width
        )
        self.label_lowbeta_power_output.grid(column=3, row=1, sticky='we')

        self.label_lowbeta_burst_output = ttk.Label(
            master=self.frame_oscillations,
            text='',
            style='outputs.TLabel',
            background='white',
            width=self.label_output_width
        )
        self.label_lowbeta_burst_output.grid(column=3, row=2, sticky='we', pady=self.label_output_pady)

        self.label_highbeta = ttk.Label(
            master=self.frame_oscillations,
            text='High-beta',
            style='outputs.TLabel',
            background=self.left_panel_bg,
            anchor='center'
        )
        self.label_highbeta.grid(column=4, row=0, sticky='we', padx=self.header_align_pad,
                                pady=self.label_output_pady)

        self.label_highbeta_power_output = ttk.Label(
            master=self.frame_oscillations,
            text='',
            style='outputs.TLabel',
            background='white',
            width=self.label_output_width
        )
        self.label_highbeta_power_output.grid(column=4, row=1, sticky='we', padx=self.header_align_pad)

        self.label_highbeta_burst_output = ttk.Label(
            master=self.frame_oscillations,
            text='',
            style='outputs.TLabel',
            background='white',
            width=self.label_output_width
        )
        self.label_highbeta_burst_output.grid(column=4, row=2, sticky='we', padx=self.header_align_pad,
                                             pady=self.label_output_pady)


        self.frame_oscillations_aux = Frame(self.tab_features, bg=self.left_panel_bg)
        self.frame_oscillations_aux.grid(column=2, row=3, sticky='nswe')

        self.label_oscillations_method = ttk.Label(
            master=self.frame_oscillations_aux,
            text='Oscillations Plot',
            style='outputs.TLabel',
            background=self.left_panel_bg,
            anchor='center'
        )
        self.label_oscillations_method.grid(column=0, row=0, sticky='we', padx=self.header_align_pad, pady=self.label_output_pady)

        # self.oscillations_method_menu = Menu(tearoff=False)
        # self.oscillations_method_menu.add_command(
        #     label='Method 1',
        #     # command=self.dropdown_save_properties,
        #     font=self.menu_item_font
        # )
        # self.oscillations_method_menu.add_command(
        #     label='Method 2',
        #     # command=self.some_function,
        #     font=self.menu_item_font
        # )
        # self.dropdown_oscillations_method = ttk.Menubutton(
        #     master=self.frame_oscillations_aux,
        #     direction='below',
        #     text='Select Method',
        #     menu=self.oscillations_method_menu,
        #     style='TMenubutton',
        # )
        # self.dropdown_oscillations_method.grid(column=0, row=1, sticky='we', padx=self.header_align_pad)

        self.button_plot_oscillations = ttk.Button(
            master=self.frame_oscillations_aux,
            text='Generate',
            style='TButton',
            command=self.button_plot_oscillations_press
        )
        self.button_plot_oscillations.grid(column=0, row=1, sticky='ew', padx=self.widget_align_pad, pady=self.label_output_pady)


        self.label_lfp = ttk.Label(
            master=self.tab_features,
            text='LFP Power',
            style='right_headers.TLabel',
            background='white'
        )
        self.label_lfp.grid(column=0, row=4, sticky='w', padx=self.header_align_pad,
                                     pady=self.features_headers_pady)

        self.frame_lfp = Frame(self.tab_features, bg=self.left_panel_bg)
        self.frame_lfp.grid(column=0, row=5, columnspan=2, sticky='nesw', padx=self.frame_padx)

        self.label_lfp_power = ttk.Label(
            master=self.frame_lfp,
            text='Power:',
            style='outputs.TLabel',
            background=self.left_panel_bg
        )
        self.label_lfp_power.grid(column=0, row=1, sticky='we', padx=self.header_align_pad,
                                           pady=self.label_output_pady)

        self.label_theta_lfp = ttk.Label(
            master=self.frame_lfp,
            text='Theta',
            style='outputs.TLabel',
            background=self.left_panel_bg,
            anchor='center'
        )
        self.label_theta_lfp.grid(column=1, row=0, sticky='we', pady=self.label_output_pady)

        self.label_theta_lfp_output = ttk.Label(
            master=self.frame_lfp,
            text='',
            style='outputs.TLabel',
            background='white',
            width=self.label_output_width
        )
        self.label_theta_lfp_output.grid(column=1, row=1, sticky='we')


        self.label_alpha_lfp = ttk.Label(
            master=self.frame_lfp,
            text='Alpha',
            style='outputs.TLabel',
            background=self.left_panel_bg,
            anchor='center'
        )
        self.label_alpha_lfp.grid(column=2, row=0, sticky='we', padx=self.header_align_pad,
                              pady=self.label_output_pady)

        self.label_alpha_lfp_output = ttk.Label(
            master=self.frame_lfp,
            text='',
            style='outputs.TLabel',
            background='white',
            width=self.label_output_width
        )
        self.label_alpha_lfp_output.grid(column=2, row=1, sticky='we', padx=self.header_align_pad)

        self.label_lowbeta_lfp = ttk.Label(
            master=self.frame_lfp,
            text='Low-beta',
            style='outputs.TLabel',
            background=self.left_panel_bg,
            anchor='center'
        )
        self.label_lowbeta_lfp.grid(column=3, row=0, sticky='we', pady=self.label_output_pady)

        self.label_lowbeta_lfp_output = ttk.Label(
            master=self.frame_lfp,
            text='',
            style='outputs.TLabel',
            background='white',
            width=self.label_output_width
        )
        self.label_lowbeta_lfp_output.grid(column=3, row=1, sticky='we')

        self.label_highbeta_lfp = ttk.Label(
            master=self.frame_lfp,
            text='High-beta',
            style='outputs.TLabel',
            background=self.left_panel_bg,
            anchor='center'
        )
        self.label_highbeta_lfp.grid(column=4, row=0, sticky='we', padx=self.header_align_pad,
                                 pady=self.label_output_pady)

        self.label_highbeta_lfp_output = ttk.Label(
            master=self.frame_lfp,
            text='',
            style='outputs.TLabel',
            background='white',
            width=self.label_output_width
        )
        self.label_highbeta_lfp_output.grid(column=4, row=1, sticky='we', padx=self.header_align_pad)


        self.frame_lfp_aux = Frame(self.tab_features, bg=self.left_panel_bg)
        self.frame_lfp_aux.grid(column=2, row=5, sticky='nesw')

        self.label_plot_lfp = ttk.Label(
            master=self.frame_lfp_aux,
            text='LFP Plot',
            style='outputs.TLabel',
            background=self.left_panel_bg,
            anchor='center'
        )
        self.label_plot_lfp.grid(column=0, row=0, sticky='we', pady=self.right_header_pady)

        self.button_plot_lfp = ttk.Button(
            master=self.frame_lfp_aux,
            text='Generate',
            style='TButton',
            command=self.button_plot_lfp_press
        )
        self.button_plot_lfp.grid(column=0, row=1, sticky='ew', padx=self.widget_align_pad, pady=self.right_header_pady)


        # self.tab_psd_plot()
        # self.plot_toolbar_psd()
        # self.tab_oscillations_plot()
        # self.plot_toolbar_oscillations()




        self.event_loop()

    def event_loop(self):
        while True:
            self.update()
            try:
                self.state()
            except:
                sys.exit()


    def dynamic_resize(self, root, rows, cols):
        for i in range(cols):
            Grid.columnconfigure(root, index=i, weight=1)
        for i in range(rows):
            Grid.rowconfigure(root, index=i, weight=1)


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
            self.update_lfp_properties(spike.main_magnitudes, spike.main_fs)


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
            self.update_lfp_properties(spike.main_magnitudes, spike.main_fs)


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
            self.update_basic_properties()
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

            self.update_basic_properties()
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


    def ask_save_filepath(self, colour=None, properties=False):
        if colour == None and properties == True:
            ask_save = tkm.askyesno(
                title='Save Properties',
                message=f'Features will be saved in a CSV file.\n\nWould you like to continue?'
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


    def button_set_threshold_press(self):
        if self.entry_threshold_factor.get() == '' or self.entry_threshold_factor.get() == 0:
            spike.main_threshold = 0
            spike.main_threshold_factor = 0
            spike.main_disp_threshold = False
            spike.main_threshold_set = False
            spike.spikesorting_plot_disp = False    # Requires valid threshold.
            self.label_threshold['text'] = f'Threshold: {round(spike.main_threshold, 2)} V ' \
                                           f'[Factor: {spike.main_threshold_factor}]'

            spike.main_event_peak_times = []
            spike.main_event_peak_mags = []

            if self.check_spikesorting_status.get() == 1:
                self.check_spikesorting.invoke()
            if self.check_spiketrain_status.get() == 1:
                self.check_spiketrain.invoke()
            if self.check_eventpeaks_status.get() == 1:
                self.check_eventpeaks.invoke()
            if self.check_thresholdbar_status.get() == 1:
                self.check_thresholdbar.invoke()

            # self.entry_lag_time.delete(0, END)
            # self.entry_lag_time.insert(0, '')
            # self.entry_lag_time['state'] = 'disabled'
            # self.entry_lag_time['style'] = 'unchecked.TEntry'
            # self.entry_time_interval.delete(0, END)
            # self.entry_time_interval.insert(0, '')
            # self.entry_time_interval['state'] = 'disabled'
            # self.entry_time_interval['style'] = 'unchecked.TEntry'
            # self.label_lag_time_interval['text'] = 'Lag: 0.5 s\nInterval: 0.01 s'

            read.threshold = None
            read.threshold_factor = None

        else:
            print('Enter threshold set process')
            spike.num_threshold_sets += 1
            spike.main_threshold_set = True
            spike.main_threshold_factor = float(self.entry_threshold_factor.get())
            print(f'Threshold factor retrieved: {spike.main_threshold_factor}')
            spike.main_threshold = spike.get_main_threshold(magnitudes=spike.main_magnitudes,
                                                 factor=spike.main_threshold_factor)
            print(f'Threshold calculated: {spike.main_threshold}')
            self.label_threshold['text'] = f'Threshold: {round(spike.main_threshold, 2)} V ' \
                                           f'[Factor: {spike.main_threshold_factor}]'

            events = spike.get_event_peaks(
                times=spike.main_times,
                magnitudes=spike.main_magnitudes,
                threshold=spike.main_threshold
            )
            print('Spike events computed.')
            spike.main_event_peak_times = events[0]
            spike.main_event_peak_mags = events[1]
            spike.features_spiketrain_indices = events[2]    # Peak indices.

            # self.entry_lag_time['state'] = 'normal'
            # self.entry_lag_time['style'] = 'TEntry'
            # self.entry_time_interval['state'] = 'normal'
            # self.entry_time_interval['style'] = 'TEntry'

            # Spike sorting update with newly set threshold.
            self.update_spikesorting_parameters()
            print('Spike sorting parameters updated')
            self.delete_plot(canvas=self.plot_canvas_spikesorting)
            self.tab_spikesorting_plot()
            self.plot_toolbar_spikesorting()
            print('SS plot updated.')

            self.update_basic_properties()
            self.update_properties()
            print('Features updated.')

            read.threshold = spike.main_threshold
            read.threshold_factor = spike.main_threshold_factor

        print('Updating main plot.')
        self.delete_plot(self.plot_canvas_main)
        self.tab_main_plot()
        self.plot_toolbar_main()
        print('Main plot updated.')


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

        self.update_basic_properties()
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

        self.update_basic_properties()
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


    def update_spikesorting_parameters(self):
        events = spike.get_event_peaks(magnitudes=spike.main_magnitudes,
                                       times=spike.main_times,
                                       threshold=spike.main_threshold)
        spike.features_spiketrain_indices = events[2]

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


    def update_lfp_properties(self, signal, fs):
        fs_lfp = 250
        raw_data_lfp = mer.get_LFP_data(signal, fs, fs_lfp)
        raw_data_lfp = (raw_data_lfp - np.mean(raw_data_lfp)) / np.std(raw_data_lfp)
        spike.lfp_theta_wave, read.lfp_theta_power = mer.get_LFP_power(raw_data_lfp, fs_lfp, 4, 8)
        spike.lfp_alpha_wave, read.lfp_alpha_power = mer.get_LFP_power(raw_data_lfp, fs_lfp, 8, 12)
        spike.lfp_low_beta_wave, read.lfp_low_beta_power = mer.get_LFP_power(raw_data_lfp, fs_lfp, 12, 21)
        spike.lfp_high_beta_wave, read.lfp_high_beta_power = mer.get_LFP_power(raw_data_lfp, fs_lfp, 21, 30)
        spike.lfp_psd_freqs, spike.lfp_psd_power = spike.get_psd(magnitudes=raw_data_lfp, fs=fs_lfp)

        self.label_theta_lfp_output['text'] = round(read.lfp_theta_power, 2)
        self.label_alpha_lfp_output['text'] = round(read.lfp_alpha_power, 2)
        self.label_lowbeta_lfp_output['text'] = round(read.lfp_low_beta_power, 2)
        self.label_highbeta_lfp_output['text'] = round(read.lfp_high_beta_power, 2)


    def update_basic_properties(self):
        # The main plot update method defines the colour-specific index arrays. Must run this method after updating main plot.
        if spike.spikesorting_plot_disp:
            if spike.selected_cluster == 'red':
                events = spike.reds_indices
            elif spike.selected_cluster == 'blue':
                events = spike.blues_indices
            elif spike.selected_cluster == 'orange':
                events = spike.oranges_indices
            elif spike.selected_cluster == 'cyan':
                events = spike.cyans_indices
            elif spike.selected_cluster == 'purple':
                events = spike.purples_indices
            elif spike.selected_cluster == 'brown':
                events = spike.browns_indices

            # This method must always be run before update_properties() in order to target correct spikes.
            spike.features_spiketrain_indices = events    # Orients analysis towards spikes obtained through sorting.

            silhouette = spike.silhouette
            read.silhouette = silhouette
            self.label_silhouette_output['text'] = round(silhouette, 2)

        elif spike.main_threshold_set:
            events = spike.features_spiketrain_indices
        else:
            self.label_snr_output['text'] = ''
            self.label_interspike_output['text'] = ''
            self.label_firing_rate_output['text'] = ''
            self.label_burst_index_output['text'] = ''
            self.label_variation_coefficient_output['text'] = ''
            return

        events = np.array(events) / spike.main_fs

        firing_rate = mer.get_FR(events)
        burst_index = mer.get_BI(events)
        cov = mer.get_CV(events)
        snr = mer.get_snr(spike.main_magnitudes, spike.main_fs, events, spike.main_times)
        isi_violations = spike.get_percent_isi_violations()

        read.snr = snr
        read.percent_isi_violations = isi_violations
        read.firing_rate = firing_rate
        read.burst_index = burst_index
        read.cov = cov

        self.label_snr_output['text'] = round(snr, 2)
        self.label_interspike_output['text'] = round(isi_violations, 2)
        self.label_firing_rate_output['text'] = round(firing_rate, 2)
        self.label_burst_index_output['text'] = round(burst_index, 2)
        self.label_variation_coefficient_output['text'] = round(cov, 2)


    def update_properties(self):
        if not spike.main_threshold_set:
            self.label_theta_power_output['text'] = ''
            self.label_alpha_power_output['text'] = ''
            self.label_lowbeta_power_output['text'] = ''
            self.label_highbeta_power_output['text'] = ''
            self.label_theta_burst_output['text'] = ''
            self.label_alpha_burst_output['text'] = ''
            self.label_lowbeta_burst_output['text'] = ''
            self.label_highbeta_burst_output['text'] = ''
            return

        events = np.array(spike.features_spiketrain_indices) / spike.main_fs
        # --- Spike train power (UNCOMMENT THIS WHEN PROBLEM FIXED)
        print('Kaneoke functions.')
        theta_spike_power = mer.kaneoke_oscillation_power(events, 4, 8, 10e-3, 500e-3)
        alpha_spike_power = mer.kaneoke_oscillation_power(events, 8, 12, 10e-3, 500e-3)
        low_beta_spike_power = mer.kaneoke_oscillation_power(events, 12, 21, 10e-3, 500e-3)
        high_beta_spike_power = mer.kaneoke_oscillation_power(events, 21, 30, 10e-3, 500e-3)

        # --- Burst duration (ALSO COMMENTED OUT)
        print('Waveform spiketrain oscillation functions.')
        spike.spike_oscillations_theta_wave, spike.spike_oscillations_theta_times, theta_wave_power = mer.waveform_spiketrain_oscillation(events, 4, 8)
        spike.spike_oscillations_alpha_wave, spike.spike_oscillations_alpha_times, alpha_wave_power = mer.waveform_spiketrain_oscillation(events, 8, 12)
        spike.spike_oscillations_low_beta_wave, spike.spike_oscillations_low_beta_times, low_beta_wave_power = mer.waveform_spiketrain_oscillation(events,
                                                                                                     12, 21)
        spike.spike_oscillations_high_beta_wave, spike.spike_oscillations_high_beta_times, high_beta_wave_power = mer.waveform_spiketrain_oscillation(
            events, 21, 30)
        print('Burst threshold function.')
        spiketrain_burst_threshold = mer.burst_threshold(events, data_type='spiketrain')

        theta_wave_envelope = mer.get_waveform_envelope(spike.spike_oscillations_theta_wave)
        alpha_wave_envelope = mer.get_waveform_envelope(spike.spike_oscillations_alpha_wave)
        low_beta_wave_envelope = mer.get_waveform_envelope(spike.spike_oscillations_low_beta_wave)
        high_beta_wave_envelope = mer.get_waveform_envelope(spike.spike_oscillations_high_beta_wave)
        print('Final burst calculation.')
        theta_spiketrain_mean_burst_duration = \
        mer.burst_features(spike.spike_oscillations_theta_times, theta_wave_envelope, spiketrain_burst_threshold)[1]
        alpha_spiketrain_mean_burst_duration = \
        mer.burst_features(spike.spike_oscillations_alpha_times, alpha_wave_envelope, spiketrain_burst_threshold)[1]
        low_beta_spiketrain_mean_burst_duration = \
        mer.burst_features(spike.spike_oscillations_low_beta_times, low_beta_wave_envelope, spiketrain_burst_threshold)[1]
        high_beta_spiketrain_mean_burst_duration = \
        mer.burst_features(spike.spike_oscillations_high_beta_times, high_beta_wave_envelope, spiketrain_burst_threshold)[1]
        # --- Burst duration (END)

        # theta_power, alpha_power, low_beta_power, high_beta_power = spike.get_wave_powers()
        self.label_theta_power_output['text'] = round(theta_spike_power, 2)
        self.label_theta_burst_output['text'] = round(theta_spiketrain_mean_burst_duration, 2)
        read.spiketrain_theta_power = theta_spike_power
        read.spiketrain_theta_burst = theta_spiketrain_mean_burst_duration
        self.label_alpha_power_output['text'] = round(alpha_spike_power, 2)
        self.label_alpha_burst_output['text'] = round(alpha_spiketrain_mean_burst_duration, 2)
        read.spiketrain_alpha_power = alpha_spike_power
        read.spiketrain_alpha_burst = alpha_spiketrain_mean_burst_duration
        self.label_lowbeta_power_output['text'] = round(low_beta_spike_power, 2)
        self.label_lowbeta_burst_output['text'] = round(low_beta_spiketrain_mean_burst_duration, 2)
        read.spiketrain_low_beta_power = low_beta_spike_power
        read.spiketrain_low_beta_burst = low_beta_spiketrain_mean_burst_duration
        self.label_highbeta_power_output['text'] = round(high_beta_spike_power, 2)
        self.label_highbeta_burst_output['text'] = round(high_beta_spiketrain_mean_burst_duration, 2)
        read.spiketrain_high_beta_power = high_beta_spike_power
        read.spiketrain_high_beta_burst = high_beta_spiketrain_mean_burst_duration


    def button_plot_isi_press(self):
        # Limit when this is called based on state of data (filtered, spike sorted, etc.)
        # Calculations of histogram will depend on cluster selected (if exists), filtering, etc.
        # Make it such that it can be called regardless of spike sorting being present.
        fig = spike.isi_plot()
        isi_win = Tk()
        isi_win.title('Inter-spike Intervals')
        isi_canvas = FigureCanvasTkAgg(fig, master=isi_win)
        isi_canvas.get_tk_widget().grid(column=0, row=0, sticky='nesw')
        isi_canvas.draw()

        toolbar_frame = Frame(isi_win)
        toolbar_frame.grid(column=0, row=1)
        toolbar = NavigationToolbar2Tk(isi_canvas, toolbar_frame)
        toolbar.configure(background='white')
        toolbar._message_label.configure(background='white')
        for i in toolbar.winfo_children():
            i.configure(background='white', bd=0)

        self.dynamic_resize(isi_win, 2, 1)
        isi_win.mainloop()


    def button_plot_oscillations_press(self):
        fig = spike.oscillations_plot()
        isi_win = Tk()
        isi_win.title('Spiketrain Oscillations')
        isi_canvas = FigureCanvasTkAgg(fig, master=isi_win)
        isi_canvas.get_tk_widget().grid(column=0, row=0, sticky='nesw')
        isi_canvas.draw()

        toolbar_frame = Frame(isi_win)
        toolbar_frame.grid(column=0, row=1)
        toolbar = NavigationToolbar2Tk(isi_canvas, toolbar_frame)
        toolbar.configure(background='white')
        toolbar._message_label.configure(background='white')
        for i in toolbar.winfo_children():
            i.configure(background='white', bd=0)

        self.dynamic_resize(isi_win, 2, 1)
        isi_win.mainloop()


    def button_plot_lfp_press(self):
        fig = spike.lfp_plot()
        lfp_win = Tk()
        lfp_win.title('LFP Oscillations')
        lfp_canvas = FigureCanvasTkAgg(fig, master=lfp_win)
        lfp_canvas.get_tk_widget().grid(column=0, row=0, sticky='nesw')
        lfp_canvas.draw()

        toolbar_frame = Frame(lfp_win)
        toolbar_frame.grid(column=0, row=1)
        toolbar = NavigationToolbar2Tk(lfp_canvas, toolbar_frame)
        toolbar.configure(background='white')
        toolbar._message_label.configure(background='white')
        for i in toolbar.winfo_children():
            i.configure(background='white', bd=0)

        self.dynamic_resize(lfp_win, 2, 1)
        lfp_win.mainloop()


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
        spike.lfp_psd_freqs = []
        spike.lfp_psd_power = []
        spike.lfp_theta_wave = []
        spike.lfp_alpha_wave = []
        spike.lfp_low_beta_wave = []
        spike.lfp_high_beta_wave = []

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

        read.snr = 0
        read.percent_isi_violations = 0
        read.firing_rate = 0
        read.burst_index = 0
        read.cov = 0
        read.silhouette = 0
        read.spiketrain_theta_power = 0
        read.spiketrain_alpha_power = 0
        read.spiketrain_low_beta_power = 0
        read.spiketrain_high_beta_power = 0
        read.spiketrain_theta_burst = 0
        read.spiketrain_alpha_burst = 0
        read.spiketrain_low_beta_burst = 0
        read.spiketrain_high_beta_burst = 0
        read.lfp_theta_power = 0
        read.lfp_alpha_power = 0
        read.lfp_low_beta_power = 0
        read.lfp_high_beta_power = 0

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
        self.label_theta_power_output['text'] = ''
        self.label_alpha_power_output['text'] = ''
        self.label_lowbeta_power_output['text'] = ''
        self.label_highbeta_power_output['text'] = ''
        self.label_theta_burst_output['text'] = ''
        self.label_alpha_burst_output['text'] = ''
        self.label_lowbeta_burst_output['text'] = ''
        self.label_highbeta_burst_output['text'] = ''
        self.label_theta_lfp_output['text'] = ''
        self.label_alpha_lfp_output['text'] = ''
        self.label_lowbeta_lfp_output['text'] = ''
        self.label_highbeta_lfp_output['text'] = ''
