import scipy
from scipy import optimize
from scipy.stats import norm
from scipy.signal import argrelextrema
from scipy.signal import find_peaks

import matplotlib
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import matplotlib.ticker

import pandas as pd

import numpy as np

#import time
#import sys

class Peak ():

    def __init__(self, minima_filename="R00029_Dark_OV2/Results.csv", verbose=False):
        # reads file and creates a table of minima
        self.table_minima = pd.read_csv(minima_filename)
        # adds a column DeltaT..
        self.table_minima.at[0,'DeltaT'] = self.table_minima.iloc[0]['Timestamp']
        #.. and calculates the DeltaT from the timestamp
        for index in self.table_minima.index[1:]:
            self.table_minima.at[index,"DeltaT"] = self.table_minima.at[index,"Timestamp"] - self.table_minima.at[index-1,"Timestamp"]
        # in a very few cases (1 case up to now), the DeltaT is < 0  (have to check why) -> set it to an infinitesimal
        self.table_minima["DeltaT"][self.table_minima["DeltaT"] < 0] = 1e-10
        # print the table
        if verbose:
            print("Table of peaks")
            print(self.table_minima.describe(),"\n")


    def exclude_bursts(self, saturating_events_only=False, n_min_burst_ev=50, verbose=False):

        #check if we want to exclude bursts after saturating events or all events,..
        #.. and define the list of events to loop over accordingly
        bursts_startsing_events = []
        if saturating_events_only:
            bursts_startsing_events = self.table_minima[self.table_minima["Amplitude"] > self.table_minima["Amplitude"].max()*0.99].index
        else:
            bursts_startsing_events = self.table_minima.index[:len(self.table_minima)-2]

        #.. then proceed to exclude the bursts if asked to
        # first creates a list of indexes of events in bursts
        burst_idx = []
        #.. and a global variable for the integrated total time of these bursts
        bursts_time = 0
            # then for every saturating peak (or all events if exclude_bursts = 1)..
        for event_idx in bursts_startsing_events:
            #for event_idx in table_minima.index[:len(table_minima)-1]:
            #.. checks if the peak index is not already in a burst..
            if event_idx in burst_idx: continue
            #.. creates a table of all following peaks..
            after_saturation = self.table_minima.loc[event_idx:]
            #.. checks what is the first peak in the table of following peaks that is more then 0.1 sec away
            if any(after_saturation["DeltaT"].loc[event_idx+1:] > 0.1):
                end_of_burst = after_saturation[after_saturation["DeltaT"] > 0.1].loc[event_idx+1:].index[0]
            else:
                end_of_burst = after_saturation.loc[event_idx+1:].index[-1]
            #.. and puts all events between the former and the latter in a table called burst
            single_burst = after_saturation.loc[event_idx:end_of_burst]
            # now if this burst has more than 10 events..
            if len(single_burst) > n_min_burst_ev:
                # it calculates the time, and adds it up to the total time of bursts
                single_burst_time = single_burst["Timestamp"].loc[end_of_burst] - single_burst["Timestamp"].loc[event_idx]
                bursts_time += single_burst_time
                # then adds the events in the list of indexes of bursts
                burst_idx += list(single_burst.index)

                print("Found a burst of",single_burst_time,"s, with",len(single_burst),"events in it -> excluded from the list of peaks")
        # print the table(s)
        if verbose:
            print("\nTable of burst peaks")
            print(self.table_minima.loc[burst_idx],"\n")
        # once the indexes of all peaks in bursts are stored, they are dropped from the global table of events
        self.table_minima = self.table_minima.drop(burst_idx)
        # now the job is done, burst peaks are removed from the stat
        if verbose:
            print("New table of peaks")
            print(self.table_minima.describe(),"\n")
        return bursts_time


    def plot_dark_count(self,ampl_bin_range,show_plot=False,save_plot=False,sipm_name="R00029",ov="2"):

        # first we need to create a dummy plot to get the binning
        fig_dummy,axes_dummy = plt.subplots()
        bin_content, bins, _ = axes_dummy.hist(self.table_minima["DeltaT"],alpha=0,bins=100)
        plt.close(fig_dummy)

        # now we do the split plot for the dark count
        fig,axes = plt.subplots(2, sharex=True, gridspec_kw={'height_ratios': [2, 1]})
        plt.xscale("log")
        figure_title = sipm_name + "  -  OV = " + str(ov) + "V"
        axes[0].set_title(figure_title)
        axes[0].scatter(self.table_minima["DeltaT"],self.table_minima["Amplitude"],marker="o",facecolors="none",edgecolors="black")
        axes[0].set_xlim([1e-8,30])
        axes[0].set_ylim(ampl_bin_range[0],ampl_bin_range[1])
        entries = r"$N_{entries}$ = " + str(len(self.table_minima))
        axes[0].text(0.02,0.88,entries,horizontalalignment="left",transform=axes[0].transAxes)
        axes[0].set_ylabel("Amplitude (V)")

        log_bins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
        bin_content, _, _ = axes[1].hist(self.table_minima["DeltaT"],bins=log_bins,histtype="step")
        bin_max = max(bin_content)*1.5
        axes[1].set_yscale("log")

        axes[1].set_ylabel("counts")
        axes[1].set_xlabel(r"$\Delta$t (s)")

        locmaj = matplotlib.ticker.LogLocator(base=10.0, numticks=10)
        axes[1].xaxis.set_major_locator(locmaj)
        axes[1].set_ylim(ymax=bin_max)
        locmin = matplotlib.ticker.LogLocator(base=10.0, subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9), numticks=10)
        axes[1].xaxis.set_minor_locator(locmin)
        axes[1].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

        plt.tight_layout()
        plt.subplots_adjust(hspace=0)

        if show_plot: plt.show()
        if save_plot:
            fig_name = "Amplitude_vs_Dt_" + sipm_name + "_OV" + ov + ".pdf"
            fig.savefig(fig_name)

    def plot_times(self,time_n_bins,show_plot=False,save_plot=False,sipm_name="R00029",ov="2"):

        fig, axes = plt.subplots(2, sharex=False, gridspec_kw={'height_ratios': [1, 1]})
        figure_title = sipm_name + "  -  OV = " + str(ov) + "V"
        axes[0].set_title(figure_title)
        axes[0].set_yscale("log")
        axes[0].set_ylabel("counts")
        axes[0].set_xlabel("Timestamp (s)")
        axes[0].hist(self.table_minima["Timestamp"],histtype="step",bins=time_n_bins)

        axes[1].set_xlabel(r"$N_{event}$")
        axes[1].set_ylabel(r"$\Delta$t (s)")
        self.table_minima = self.table_minima.reset_index()
        axes[1].step(self.table_minima.index,self.table_minima["DeltaT"])#,marker=".")#,facecolors="none",edgecolors="black")

        plt.tight_layout()
        if show_plot: plt.show()
        if save_plot:
            fig_name = "Timestamp_Dt_" + sipm_name + "_OV" + ov + ".pdf"
            fig.savefig(fig_name)

    def plot_amplitude(self,ampl_n_bins=1000,show_plot=False,save_plot=False,sipm_name="R00029",ov="2",do_fit=True,save_fit_results=False):
        save_fit_results = do_fit and save_fit_results

        fig, axes = plt.subplots()
        figure_title = sipm_name + "  -  OV = " + str(ov) + "V"
        axes.set_title(figure_title)
        axes.set_yscale("log")
        axes.set_ylabel("counts")
        axes.set_xlabel("Amplitude (V)")
        # gets the bin content and center to find maxima..
        bin_content, bin_edge, _ = axes.hist(self.table_minima["Amplitude"],histtype="step",bins=ampl_n_bins)
        # calculate bin center from bin edge
        bin_width = bin_edge[1]-bin_edge[0]
        bin_center = bin_edge+bin_width/2.
        # adjust length (bin center was taken from edges so it's n+1)
        bin_center = bin_center[:-1]
        bin_center_dictionary = {}
        # creates a dictionary to retrieve bin index from its central calue
        for idx in range(len(bin_center)):
            bin_center_dictionary[bin_center[idx]] = idx

        # define a function that recursively looks for maxima in sliding slices of the histogram
        def find_hist_peaks(bin_center, bin_content):
            if len(bin_content) == 0 or max(bin_content) < 3:
                return

            max_bin = bin_center[bin_content==max(bin_content[0:ampl_n_bins/10])][0]
            max_bin_center.append(max_bin)

            start_searching_from = max_bin + max_bin_center[0]/2.
            new_bin_center = bin_center[bin_center>start_searching_from]
            new_bin_content = bin_content[bin_center>start_searching_from]

            find_hist_peaks(new_bin_center, new_bin_content)

        # now create a list of bins, call the function to fill it, and retrieve the list of indexes using the dictionary
        max_bin_center = []
        find_hist_peaks(bin_center, bin_content)
        max_list = [bin_center_dictionary[bc] for bc in max_bin_center]

        #.. and finaly plots the peaks as a scatter plot
        axes.scatter(bin_center[max_list],bin_content[max_list], marker='x', s=10, c='g')
        # before plotting fit on top of the histogram, retrievew y axis range
        my_ylim = axes.get_ylim()

        # define here our gaussian function
        def gaussian(x, amplitude, mean, stdev):
            return amplitude / stdev * np.sqrt(2*np.pi) * np.exp( - ( (x - mean) / (stdev*np.sqrt(2)) )**2 )

        # define the index width in which to fit every peak
        index_width = int((max_list[1]-max_list[0])/2)

        if do_fit:
            # open a file to store fit DataFrame
            if save_fit_results:
                ofile_name =  "pe_fit_" + sipm_name + "_OV" + ov + ".txt"
                my_ofile = open(ofile_name,"w")
                print >> my_ofile, "pe mu sigma"

            # now loop over the maxima and fit them
            fit_counter = 1
            for max_index in max_list:
                # define a range of bins to fit for every maximum
                bin_range = [max(max_index-index_width,0),min(max_index+index_width,len(bin_center)-1)]
                fit_bin_range = bin_center[bin_range[0]:bin_range[1]]
                fit_bin_values = bin_content[bin_range[0]:bin_range[1]]
                # then fit them and plot the fit
                try:
                    #fit_bounds = [ [0,np.inf],[min(fit_bin_range),max(fit_bin_range)],[0,np.inf] ]
                    fit_bounds = [ [0,min(fit_bin_range),0],[0.1,max(fit_bin_range),np.inf] ]
                    popt, _ = optimize.curve_fit(gaussian, fit_bin_range, fit_bin_values, maxfev=10000, bounds=fit_bounds)
                    if save_fit_results: print >> my_ofile, fit_counter, popt[1], popt[2]
                    axes.plot(fit_bin_range, gaussian(fit_bin_range, *popt), linewidth=1, color="orange")
                except: print("Failed at fitting peak on bin",max_index,"!")
                fit_counter += 1

            if save_fit_results:
                print("Fitting results for different pe peaks saved in ",ofile_name,"\n")
                my_ofile.close()
        # here it reapplies the range found before the gaussian fits were made
        axes.set_ylim(my_ylim)
        # and finally shows and saves the plot
        plt.tight_layout()
        if show_plot: plt.show()
        if save_plot:
            fig_name = "Amplitude_" + sipm_name + "_OV" + ov + ".pdf"
            fig.savefig(fig_name)
