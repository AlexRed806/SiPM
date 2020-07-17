import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import argrelextrema
from pandas import DataFrame
import time
from time import sleep
import progressbar
import sys
import os
# import statistics
# from statistics import mean

class Waveform ():

    number_of_events = 0

    def __init__(self, timestamp_filename="R00029_Dark_OV2/TimeTrend_R00029_Cold_26-06-2020_10ev.csv",\
     waveform_filename="R00029_Dark_OV2/WaveForms_R00029_Cold_26-06-2020_10ev.csv",\
     wf_data_points=6250):
        # inport input files
        self.table_timestamp = pd.read_csv(timestamp_filename)

        # changes the header of the timestamp file
        #if self.table_timestamp.columns[0] == "X: (s)":
        self.table_timestamp.rename(columns={'X: (s)':'Event','Y: (Hits)':'DeltaT'},inplace=True)
        self.number_of_events = len(self.table_timestamp)

        # searches the waveform file for header
        waveform_file = open(waveform_filename)
        line_counter = 0
        n_line = -1
        while n_line == -1:
            line = waveform_file.readline()
            if line.startswith("TIME"): n_line = line_counter
            line_counter += 1
            if line_counter == 100:
                print("ERROR: header not found in", waveform_filename, "after the first 100 lines!")
                sys.exit(1)
        waveform_file.close()
        # then stores the waveform in a table
        self.table_waveform = pd.read_csv(waveform_filename,header=n_line if n_line==0 else n_line-1)

        #add column for absolute timestamp
        self.table_timestamp.at[0,'Timestamp'] = self.table_timestamp.iloc[0]['DeltaT']

        #calculating the timptamp
        for index in self.table_timestamp.index[1:]:
            self.table_timestamp.at[index,"Timestamp"] = self.table_timestamp.at[index-1,"Timestamp"] + self.table_timestamp.at[index,"DeltaT"]
        #defining the number of data points in a waveform
        self.wf_data_points = wf_data_points
        # wf_data_points = 1250 #First test
        # wf_data_points = 6250 #RoomT || Dark
        # wf_data_points = 1563 #DAq Room
        self.table_global = pd.DataFrame()

        #plt.plot(table_timestamp.index, table_timestamp['DeltaT'])
        #plt.axis([-100, 1100, -0.0000001, 0.000005]) # Set axis limits
        #plt.show()

    def event_minima(self, event, bsl_method, bsl_n_points, min_method, min_search_range, min_back_shift, min_n_points, min_gap, min_n_close, show_plot, save_plot, folder_name):
        ''' find minima of a single event (a snippet of the waveform table) '''
        do_plot = show_plot or save_plot
        if do_plot: fig, ax = plt.subplots()

        ev_waveform = self.table_waveform.loc[event*self.wf_data_points : event*self.wf_data_points+self.wf_data_points-1]
        #ev_waveform["CH1"] = ev_waveform["CH1"].replace(to_replace="-inf",value=ev_waveform.min())
        ev_waveform["CH1"].replace({(float)("-inf"): np.nanmin(ev_waveform[ev_waveform != -np.inf])}, inplace=True)

        minimum_list = argrelextrema(ev_waveform.CH1.values, np.less_equal, order=min_search_range)[0]
        minimum_list = list(minimum_list)
        #THIS IS WHERE THE WARNING IS
        #ev_waveform.at[minimum_list, 'min'] = ev_waveform.loc[minimum_list, 'CH1']
        #ev_waveform.loc[minimum_list].at['min'] = ev_waveform.loc[minimum_list, 'CH1']
        #print(ev_waveform.at[minimum_list,"CH1"])
        ev_waveform.loc[:,'min'] = ev_waveform.iloc[minimum_list]['CH1']
        #ev_waveform.at[:,'min'] = ev_waveform.loc[minimum_list,'CH1']
        #ev_waveform.insert(loc=minimum_list, column="min", value="CH1")

        baseline = 0
        if bsl_method==0: baseline = find_baseline_method_0(ev_waveform,bsl_n_points)
        elif bsl_method==1: baseline = find_baseline_method_1(ev_waveform,bsl_n_points)

        clean_minimum_list = []
        if min_method==0: clean_minimum_list = clean_minima_method_0(ev_waveform,minimum_list,min_back_shift,min_n_points,min_gap)
        elif min_method==1: clean_minimum_list = clean_minima_method_1(ev_waveform,minimum_list,baseline,min_gap,min_n_close)

        ev_waveform.loc[:,'clean_min'] = ev_waveform.iloc[clean_minimum_list]['CH1']

        #produce the waveform plots
        if do_plot:
            #ax.plot(ev_waveform.index, ev_waveform['CH1'], marker="2", linestyle="-", linewidth=1)
            ax.plot(ev_waveform.index, ev_waveform['CH1'], linestyle="-", linewidth=1)
            ax.scatter(ev_waveform.index, ev_waveform['min'], marker='o', s=10, c='r')
            ax.scatter(ev_waveform.index, ev_waveform['clean_min'], marker='o', s=50, c='g')
            ax.axhline(baseline, c='b')
            if show_plot: plt.show()

            if save_plot:
                name = '{0}/ev{1}.png'.format(folder_name,event)
                fig.savefig(name)
            plt.close(fig)

        #replace nana with a number and search for those events
        ev_waveform = ev_waveform.fillna(-12345)
        ev_waveform = ev_waveform[ ev_waveform['clean_min'] != -12345 ]
        ev_waveform = ev_waveform.drop_duplicates(subset='clean_min',keep='first')

        #adding the event number as a cloumn
        ev_waveform["Event"] = event
        #adding the Timestamp
        ev_waveform["Timestamp"] = float(self.table_timestamp[ self.table_timestamp['Event'] == event ]["Timestamp"]) + ev_waveform["TIME"]
        #adding the baseline
        ev_waveform["Baseline"] = baseline
        #adding the amplitude
        ev_waveform["Amplitude"] = -(ev_waveform["clean_min"] - baseline)

        return ev_waveform


    def find_all_minima(self, ev_list, bsl_method, bsl_n_points, min_method, min_search_range=50, min_back_shift=100, min_n_points=100, min_gap=0.003, min_n_close = 100, show_plot=False, save_plot=False, folder_name='plots'):

        if ev_list == []:
            ev_list = [int(i) for i in range(len(self.table_timestamp))]
        n_ev = len(ev_list)

        if save_plot and not os.path.exists(folder_name): os.mkdir(folder_name)

        print("-------------------------------------")
        print("Analyzing the waveforms to get minima")
        print("-------------------------------------")
        counter = 0
        bar = progressbar.ProgressBar(maxval=n_ev, \
        widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
        bar.start()
        #for ev in self.table_timestamp["Event"]:
        for event in ev_list:

            ev_waveform = self.event_minima(event, bsl_method, bsl_n_points, min_method, min_search_range, min_back_shift, min_n_points, min_gap, min_n_close, show_plot, save_plot, folder_name)
            self.table_global = pd.concat([self.table_global,ev_waveform],ignore_index=True)
            bar.update(counter+1)
            sleep(0.1)
            #if counter == n_ev-1: break
            counter+=1
        bar.finish()
        print("-------------------------------------\n")


    def save_minimum_table(self, table_name):
        pd.options.display.precision = 12
        '''name of the table without extention (wanna save in csv and txt)'''
        self.table_global = self.table_global[['Timestamp','Baseline','Amplitude']]
        print("-------------------------------------")
        print("           Table of minima           ")
        print("-------------------------------------")
        print(self.table_global)
        print("-------------------------------------")
        self.table_global.to_csv(table_name+".csv", header=True, index=False, sep=',')
        np.savetxt(table_name+".txt", self.table_global.values, fmt='%4.12f', delimiter=' ')



def clean_minima_method_0(waveform,index_list,n_previous=50,n_check=100,gap=0.003):
    '''
    waveform [DataFrame]: waveform to analyze
    index_list [list]: list of the index minima to clean
    '''
    idx_list_clean = []
    for index in index_list:
        if index-n_previous < 0:
            continue
        idx_start = 0 if index-n_previous-n_check < 0 else index-n_previous-n_check
        waveform_subset = waveform.iloc[idx_start+1:index-n_previous]

        waveform_subset["CH1"] = waveform_subset["CH1"].astype(float)

        if ( abs(waveform_subset["CH1"].mean() - (float)(waveform["CH1"].iloc[index])) > gap ):
            idx_list_clean.append(index)
    return idx_list_clean

def clean_minima_method_1(waveform,index_list,baseline,gap=0.005,n_close=100):
    idx_list_above_bsl = []
    idx_list_clean = []
    for index in index_list:
        if baseline - waveform["CH1"].iloc[index] > gap:# and len(waveform[index-1:index])-1 > n_close:
            idx_list_above_bsl.append(index)
    return idx_list_above_bsl


def find_baseline_method_0(waveform,n_firsts=100):
    waveform_start = waveform.iloc[0:n_firsts]
    return np.polyfit(waveform_start["TIME"],waveform_start["CH1"],0)[0]

def find_baseline_method_1(waveform,n_before=100):
    maxch1 = waveform["CH1"].loc[n_before:].idxmax()
    waveform_near_max = waveform.loc[maxch1-n_before:maxch1]
    return np.polyfit(waveform_near_max["TIME"],waveform_near_max["CH1"],0)[0]
