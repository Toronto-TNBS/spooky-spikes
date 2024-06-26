from PySide6.QtWidgets import *
from slots.backend import data
from slots.backend import files


FILEDATA = None

def button_file_clicked(app):
    global FILEDATA

    filetypes = 'Spike2 Datafiles (*.smr *.smrx);;MAT-Files (*.mat);;Comma Separated Values (*.csv)'
    filepath = QFileDialog.getOpenFileName(caption='Open File', dir='.', filter=filetypes)[0]
    if filepath == '':
        return

    file_data = files.load_spike2(filepath)
    FILEDATA = file_data

    app.label_file.setText(file_data.filepath.split('/')[-1])

    app.dropdown_channel.addItems([f'Ch {ch_id}' for ch_id in FILEDATA.channels.keys()])
    