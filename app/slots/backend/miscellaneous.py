

def check_convert_string(string, type):
        try:
            return type(string)
        except:
            return None


def remove_plot_item(plotitem, itemtype):
    for item in plotitem.listDataItems():
        if type(item) == itemtype:
            plotitem.removeItem(item)