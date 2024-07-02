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


def plot_spike_events(plotitem, spike_times):
    for event in spike_times:
        event_item = pg.InfiniteLine(pos=event, angle=90, pen=pg.mkPen('r', width=1))
        plotitem.addItem(event_item)    