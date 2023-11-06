from sonpy import lib as sp
import numpy as np
import math
import pandas as pd
import scipy.io as spio
import mat73


class ReadSpike:

    def __init__(self):
        self.read_filepath = ''
        self.channel = 0
        self.save_filepath = ''

        self.snr = 0
        self.percent_isi_violations = 0
        self.firing_rate = 0
        self.burst_index = 0
        self.cov = 0
        self.silhouette = 0
        self.spiketrain_theta_power = 0
        self.spiketrain_alpha_power = 0
        self.spiketrain_low_beta_power = 0
        self.spiketrain_high_beta_power = 0
        self.spiketrain_theta_burst = 0
        self.spiketrain_alpha_burst = 0
        self.spiketrain_low_beta_burst = 0
        self.spiketrain_high_beta_burst = 0
        self.lfp_theta_power = 0
        self.lfp_alpha_power = 0
        self.lfp_low_beta_power = 0
        self.lfp_high_beta_power = 0

        self.threshold = None
        self.threshold_factor = None
        self.filtering = False
        self.highpass = None
        self.lowpass = None
        self.segment_inverted = False
        self.number_clusters = None
        self.lag_time = 0.5
        self.time_interval = 0.01

        self.supported_filetypes = ['smr', 'smrx', 'mat', 'csv']
        self.file_ext = None
        self.available_channels = None


    def get_channel_list(self, ext):
        if ext == 'smr' or ext == 'smrx':
            data = sp.SonFile(sName=self.read_filepath, bReadOnly=True)
            channel_list = [f'{i + 1} ({str(data.ChannelType(i)).split(".")[-1]})' for i in range(data.MaxChannels())]

        elif ext == 'mat':
            try:
                data = spio.loadmat(self.read_filepath)
            except:
                data = mat73.loadmat(self.read_filepath)

            channel_list = []
            for var in data.keys():
                if var[:6] == 'values':
                    try:
                        type_test = type(int(var[6:]))
                    except:
                        continue
                    channel_list.append(var)
            self.channel = 0

        elif ext == 'csv':
            data = pd.read_csv(self.read_filepath)
            channel_list = []
            for var in data.columns:
                if var[:6] == 'values':
                    try:
                        type_test = type(int(var[6:]))
                    except:
                        continue
                    channel_list.append(var)
            self.channel = 0

        return channel_list


    def get_lowest_wave_channel(self):
        file = sp.SonFile(self.read_filepath, True)
        wave_channels = [i for i in range(file.MaxChannels()) if file.ChannelType(i) == sp.DataType.Adc]
        for i in range(file.MaxChannels()):
            if file.ChannelType(i) == sp.DataType.RealWave:
                wave_channels.append(i)
        if len(wave_channels) == 0:
            print('Error: no wave spikedata found in this file.')
        else:
            self.channel = min(wave_channels)


    def read_smr_waves(self):
        file = sp.SonFile(sName=self.read_filepath, bReadOnly=True)

        if file.ChannelType(self.channel) == sp.DataType.Adc or file.ChannelType(self.channel) == sp.DataType.RealWave:
            read_max_time = file.ChannelMaxTime(self.channel) * file.GetTimeBase()

            period = file.ChannelDivide(self.channel) * file.GetTimeBase()
            num_points = math.floor(read_max_time / period)
            fs = 1 / period

            wave = file.ReadFloats(self.channel, num_points, 0)

            if len(wave) == 0:
                print('There was no spikedata read.')
            elif len(wave) != num_points:
                print(f'Mismatched number of points. Expected {num_points} points, but got {len(wave)} instead.')

            return wave, fs
        else:
            return -1


    def read_mat_wave(self):
        try:
            data = spio.loadmat(self.read_filepath)
        except:
            data = mat73.loadmat(self.read_filepath)
        var = self.available_channels[self.channel]
        id = var[6:]
        wave = np.ndarray.flatten(data[var])
        fs = float(np.ndarray.flatten(data['fs' + id])[0])
        return wave, fs


    def read_csv_wave(self):
        data = pd.read_csv(self.read_filepath)
        var = self.available_channels[self.channel]
        id = var[6:]
        wave = data[var].values
        fs = float(data['fs' + id][0])
        return wave, fs


    def read_wave(self, ext):
        if ext == 'smr' or ext == 'smrx':
            wave_data = self.read_smr_waves()
        elif ext == 'mat':
            wave_data = self.read_mat_wave()
        elif ext == 'csv':
            wave_data = self.read_csv_wave()

        if wave_data == -1:
            return -1
        return wave_data[0], wave_data[1]


    def save_segment(self, events):
        file = sp.SonFile(sName=self.read_filepath, bReadOnly=True)

        read_max_time = file.ChannelMaxTime(self.channel) * file.GetTimeBase()

        period = file.ChannelDivide(self.channel) * file.GetTimeBase()
        num_points = math.floor(read_max_time / period)

        time_points = np.arange(0, num_points * period, period)
        wave = file.ReadInts(self.channel, num_points, 0)

        if len(wave) == 0:
            print('There was no data read.')
            quit()
        elif len(wave) != num_points:
            print(f'Mismatched number of points. Expected {num_points} points, but got {len(wave)} instead.')
            quit()

        filename = self.save_filepath
        new = sp.SonFile(filename)

        wave_offset = file.GetChannelOffset(self.channel)
        wave_scale = file.GetChannelScale(self.channel)
        wave_time_base = file.GetTimeBase()
        wave_Y_range = file.GetChannelYRange(self.channel)
        wave_Y_low = wave_Y_range[0]
        wave_Y_high = wave_Y_range[1]
        wave_units = file.GetChannelUnits(self.channel)
        wave_title = file.GetChannelTitle(self.channel)
        wave_channel_type = file.ChannelType(self.channel)
        wave_sample_ticks = file.ChannelDivide(self.channel)
        wave_Fs = 1 / period
        time_start = 0

        new.SetTimeBase(wave_time_base)
        new_wave_channel = 0

        new.SetWaveChannel(new_wave_channel, wave_sample_ticks, sp.DataType.Adc, wave_Fs)
        new.SetChannelTitle(new_wave_channel, wave_title)
        new.SetChannelUnits(new_wave_channel, wave_units)
        new.SetChannelScale(new_wave_channel, wave_scale)
        new.SetChannelOffset(new_wave_channel, wave_offset)
        new.SetChannelYRange(new_wave_channel, wave_Y_low, wave_Y_high)

        spike_channel = new.GetFreeChannel()
        max_event_rate = 1 / (wave_time_base * wave_sample_ticks)

        new.SetEventChannel(spike_channel, max_event_rate, sp.DataType.EventFall)
        new.SetChannelTitle(spike_channel, 'ST-1')
        new.SetChannelUnits(spike_channel, wave_units)
        new.SetChannelScale(spike_channel, wave_scale)
        new.SetChannelOffset(spike_channel, wave_offset)
        new.SetChannelYRange(spike_channel, wave_Y_low, wave_Y_high)
        del new

        reload = sp.SonFile(filename, False)
        if self.segment_inverted:
            wave = [-i for i in wave]
        write_wave = reload.WriteInts(new_wave_channel, wave, time_start)

        if write_wave < 0:
            print(f'Error code: {write_wave}')
            print(sp.GetErrorString(write_wave))
        else:
            print('Waveform written successfully.')

        spike_ticks = [int(i / wave_time_base) for i in events]
        spike_times_formatted = np.array(spike_ticks, dtype=np.int64)
        write_train = reload.WriteEvents(spike_channel, spike_times_formatted)

        if write_train < 0:
            print(f'Error Code: {write_train}')
            print(sp.GetErrorString(write_train))
            status = -1
        else:
            print('Spike train written successfully.')
            status = 0

        del reload

        return status

    def save_properties(self):
        properties = pd.DataFrame(
            {
                'Channel': self.channel,
                'Threshold': self.threshold,
                'Threshold_Factor': self.threshold_factor,
                'Filtering': self.filtering,
                'Highpass': self.highpass,
                'Lowpass': self.lowpass,
                'Inverted': self.segment_inverted,
                'Number_Clusters': self.number_clusters,
                'SNR': self.snr,
                'Percent_ISI_Violations': self.percent_isi_violations,
                'Firing_Rate': self.firing_rate,
                'Burst_Index': self.burst_index,
                'CoV': self.cov,
                'Silhouette': self.silhouette,
                'Spiketrain_Theta_Power': self.spiketrain_theta_power,
                'Spiketrain_Alpha_Power': self.spiketrain_alpha_power,
                'Spiketrain_Low_Beta_Power': self.spiketrain_low_beta_power,
                'Spiketrain_High_Beta_Power': self.spiketrain_high_beta_power,
                'Theta_Mean_Burst_Duration': self.spiketrain_theta_burst,
                'Alpha_Mean_Burst_Duration': self.spiketrain_alpha_burst,
                'Low_Beta_Mean_Burst_Duration': self.spiketrain_low_beta_burst,
                'High_Beta_Mean_Burst_Duration': self.spiketrain_high_beta_burst,
                'LFP_Theta_Power': self.lfp_theta_power,
                'LFP_Alpha_Power': self.lfp_alpha_power,
                'LFP_Low_Beta_Power': self.lfp_low_beta_power,
                'LFP_High_Beta_Power': self.lfp_high_beta_power
            },
            index=[0]
        )
        properties.to_csv(self.save_filepath)

        return 0
