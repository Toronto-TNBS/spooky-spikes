from PySide6 import QtGui
import pyqtgraph as pg


def check_convert_string(string, type):
        try:
            return type(string)
        except:
            return None


def remove_plot_item(plotitem, itemtype):
    for item in plotitem.listDataItems():
        if type(item) == itemtype:
            plotitem.removeItem(item)


def vertical_line_symbol():
        symbol = QtGui.QPainterPath()
        symbol.moveTo(0, -1)
        symbol.lineTo(0, 1)
        return symbol