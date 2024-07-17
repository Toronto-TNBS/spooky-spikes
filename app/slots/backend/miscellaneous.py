from PySide6 import QtGui
import pyqtgraph as pg


def check_convert_string(string, type):
        try:
            return type(string)
        except:
            return None


def remove_plot_items(plotitem, itemtype):
    for item in plotitem.listDataItems():
        if type(item) == itemtype:
            plotitem.removeItem(item)


def vertical_line_symbol():
        symbol = QtGui.QPainterPath()
        symbol.moveTo(0, -1)
        symbol.lineTo(0, 1)
        return symbol


def reset_spikesorting_clusters_data(app):
     for key in app.channeldata.spikesorting_clusters.keys():
          app.channeldata.spikesorting_clusters[key] = None
        

def get_existing_clusters(app):
     return [(i, key) for i, key in enumerate(app.channeldata.spikesorting_clusters.keys()) if app.channeldata.spikesorting_clusters[key] is not None]