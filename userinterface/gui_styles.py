import tkinter.ttk as ttk
from tkinter import *

class GUIStyles(Tk):

    def __init__(self):
        super().__init__()

        self.style = ttk.Style()
        self.style.theme_use('clam')    # if you want different styles for different widget types, use more objects

        self.left_panel_bg = '#EDEDED' #'#E5E5E5'
        self.right_panel_bg = 'white'
        self.accent_colour = '#205375'
        self.header_align_pad = 20
        self.right_header_pady = 15
        self.widget_align_pad = 30
        self.widget_thick_pad = 7
        self.label_output_width = 7
        self.label_output_name_width = 23
        self.label_output_pady = 1
        self.menu_item_font = ('Arial', 14, 'normal')

        self.style.configure(
            'headers.TLabel',
            fg='black',
            padding=3,
            font=('Arial', 16, 'bold'),
        )

        self.style.configure(
            'right_headers.TLabel',
            fg='black',
            padding=0,
            font=('Arial', 14, 'bold'),
        )

        self.style.configure(
            'file.TLabel',
            padding=self.widget_thick_pad,
            font=('Arial', 12, 'normal')

        )

        self.style.configure(
            'outputs.TLabel',
            font=('Arial', 14, 'normal'),
        )

        self.style.configure(
            'right_outputs.TLabel',
            font=('Arial', 12, 'normal'),
        )

        self.style.configure(
            'TButton',
            foreground='white',
            background=self.accent_colour,
            highlightthickness=0,
            bordercolor='white',
            relief='flat',
            padding=self.widget_thick_pad,
            font=('Arial', 14, 'bold')
        )

        self.style.configure(
            'TMenubutton',
            padding=self.widget_thick_pad,
            relief='flat',
            arrowcolor=self.accent_colour,
            background='white',
            arrowsize=10,
            fieldbackground='white',
            font=('Arial', 12, 'normal')
        )

        self.style.configure(
            'clusters.TMenubutton',
            padding=self.widget_thick_pad,
            relief='flat',
            arrowcolor=self.accent_colour,
            background=self.left_panel_bg,
            arrowsize=10,
            fieldbackground='red',
            font=('Arial', 12, 'normal')
        )

        # docs says some options for configure method only available for certain themes. Some are licensed.
        self.style.configure(
            'left_check.TCheckbutton',
            background=self.left_panel_bg,
            font=('Arial', 14, 'normal'),
            relief='flat'
        )

        self.style.configure(
            'right_check.TCheckbutton',
            background=self.right_panel_bg,
            font=('Arial', 14, 'normal'),
            relief='flat'
        )

        self.style.configure(
            'TNotebook',
            background='white',
            bordercolor='white',
            tabmargins=1,
            tabposition='n',
        )

        self.style.configure(
            'TNotebook.Tab',
            bordercolor='white',
            background=self.accent_colour,
            relief='flat',
            padding=[43, 15],
            foreground='white',
            font=('Arial', 14, 'normal'),
            width=15
        )
        self.style.map(
            'TNotebook.Tab',
            background=[('selected', '#2D4059')],
            padding=[('selected', [43, 15])],
            font=[('selected', ('Arial', 14, 'bold'))]
        )

        self.style.configure(
            'TEntry',
            padding=self.widget_thick_pad,
            relief='flat'
        )

        self.style.configure(
            'unchecked.TEntry',
            padding=self.widget_thick_pad,
            relief='flat',
            fieldbackground=self.left_panel_bg
        )