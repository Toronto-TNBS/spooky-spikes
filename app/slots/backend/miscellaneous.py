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


def update_unit_features_display(app, ss=''):
     app.label_tab_features_patterned_firingrate_status.setText(str(round(app.channeldata.features_patterned['FR'], 2)))
     app.label_tab_features_patterned_burstindex_status.setText(str(round(app.channeldata.features_patterned['BI'], 2)))
     app.label_tab_features_patterned_cov_status.setText(str(round(app.channeldata.features_patterned['CV'], 2)))

     app.label_tab_features_qualitymetrics_snr_status.setText(str(round(app.channeldata.features_qualitymetrics['SNR'], 2)))
     app.label_tab_features_qualitymetrics_isi_status.setText(str(round(app.channeldata.features_qualitymetrics['ISI'], 2)))
     if ss != '':
          ss = str(round(app.channeldata.features_qualitymetrics['SS'], 2))
     app.label_tab_features_qualitymetrics_silhouette_status.setText(ss)

     app.label_tab_features_spiketrain_thetapower_status.setText(str(round(app.channeldata.features_spiketrain['Power']['Theta'], 2)))
     app.label_tab_features_spiketrain_alphapower_status.setText(str(round(app.channeldata.features_spiketrain['Power']['Alpha'], 2)))
     app.label_tab_features_spiketrain_lowbetapower_status.setText(str(round(app.channeldata.features_spiketrain['Power']['LowBeta'], 2)))
     app.label_tab_features_spiketrain_highbetapower_status.setText(str(round(app.channeldata.features_spiketrain['Power']['HighBeta'], 2)))

     app.label_tab_features_spiketrain_thetaburstduration_status.setText(str(round(app.channeldata.features_spiketrain['BurstDur']['Theta'], 2)))
     app.label_tab_features_spiketrain_alphaburstduration_status.setText(str(round(app.channeldata.features_spiketrain['BurstDur']['Alpha'], 2)))
     app.label_tab_features_spiketrain_lowbetaburstduration_status.setText(str(round(app.channeldata.features_spiketrain['BurstDur']['LowBeta'], 2)))
     app.label_tab_features_spiketrain_highbetaburstduration_status.setText(str(round(app.channeldata.features_spiketrain['BurstDur']['HighBeta'], 2)))


def update_lfp_features_display(app):
     app.label_tab_features_lfp_thetapower_status.setText(str(round(app.channeldata.features_lfp['Power']['Theta'], 2)))
     app.label_tab_features_lfp_alphapower_status.setText(str(round(app.channeldata.features_lfp['Power']['Alpha'], 2)))
     app.label_tab_features_lfp_lowbetapower_status.setText(str(round(app.channeldata.features_lfp['Power']['LowBeta'], 2)))
     app.label_tab_features_lfp_highbetapower_status.setText(str(round(app.channeldata.features_lfp['Power']['HighBeta'], 2)))