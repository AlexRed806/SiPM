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

import time

import sys

if len(sys.argv) != 5:
    print "\nUSAGE: python plotter.py file_path simpm_name ov exclude_bursts"
    print "\nARGUMENTS:"
    print "file_path:\t\t path to .csv file with timestamps and amplitude of peaks (e.g. \"R00029_Dark/OV2/Results.csv\")"
    print "sipm_name:\t\t SiPM name, used for plot title (e.g. \"R00029\")"
    print "ov:\t\t\t over voltage value, used for plot title (2-5)"
    print "exclude_bursts:\t\t 0 - plot all events; 1 - exclude event bursts after saturation; 2 - exclude all bursts\n"
    sys.exit()

file_path = sys.argv[1] # path to file with timestamps and amplitude of peaks (e.g. "R00029_Dark/OV2/Results.txt")
sipm_name = sys.argv[2] # SiPM name, used for pplot title (e.g. "R00029")
ov = sys.argv[3] # over voltage value (2-5)

# reads file and creates a table of minima
table_minima = pd.read_csv(file_path)

# adds a column DeltaT..
table_minima.at[0,'DeltaT'] = table_minima.iloc[0]['Timestamp']
#.. and calculates the DeltaT from the timestamp
for index in table_minima.index[1:]:
    table_minima.at[index,"DeltaT"] = table_minima.at[index,"Timestamp"] - table_minima.at[index-1,"Timestamp"]

# in a very few cases (1 case up to now), the DeltaT is < 0  (have to check why) -> set it to an infinitesimal
table_minima["DeltaT"][table_minima["DeltaT"] < 0] = 1e-10

# check if we have to exclude bursts from the table..
exclude_bursts = 0
saturating_events = []
if (int)(sys.argv[4]) == 1:
    saturating_events = table_minima[table_minima["Amplitude"] > table_minima["Amplitude"].max()*0.99].index
    exclude_bursts = 1
elif (int)(sys.argv[4]) == 2:
    saturating_events = table_minima.index[:len(table_minima)-2]
    exclude_bursts = 1
elif (int)(sys.argv[4]) == 0:
    exclude_bursts = 0
else:
    print "ERROR: invalid value for exclude_bursts parameter, chose 0, 1, or 2"
    sys.exit()

#.. then proceed to exclude them if asked to

# first creates a list of indexes of events in bursts
burst_idx = []
#.. and a global variable for the integrated total time of these bursts
bursts_time = 0

if exclude_bursts:
    # then for every saturating peak (or all events if exclude_bursts = 1)..
    for my_event_idx in saturating_events:
        #for my_event_idx in table_minima.index[:len(table_minima)-1]:
        #.. checks if the peak index is not already in a burst..
        if my_event_idx in burst_idx: continue
        #.. creates a table of all following peaks..
        after_saturation = table_minima.loc[my_event_idx:]
        #.. checks what is the first peak in the table of following peaks that is more then 0.1 sec away
        end_of_burst = after_saturation[after_saturation["DeltaT"] > 0.1].loc[my_event_idx+1:].index[0]
        #.. and puts all events between the former and the latter in a table called burst
        single_burst = after_saturation.loc[my_event_idx:end_of_burst]
        # now if this burst has more than 10 events..
        if len(single_burst) > 50:
            # it calculates the time, and adds it up to the total time of bursts
            single_burst_time = single_burst["Timestamp"].loc[end_of_burst] - single_burst["Timestamp"].loc[my_event_idx]
            bursts_time += single_burst_time
            # then adds the events in the list of indexes of bursts
            burst_idx += list(single_burst.index)

            print "Found a burst of",single_burst_time,"s, with",len(single_burst),"events in it -> excluded from the list of peaks"

    print "\nTable of burst peaks"
    print table_minima.loc[burst_idx],"\n"

    print "Total integrated time of bursts:",bursts_time,"\n"

    # once the indexes of all peaks in bursts are stored, they are dripped from the global table of events
    table_minima = table_minima.drop(burst_idx)

# now the job is done, burst peaks are removed from the stat
#table_minima = table_minima.astype(float)
print table_minima.describe(),"\n"

# now we calculate the dark count rate
dark_count_rate = len(table_minima[ (table_minima["DeltaT"]>1e-6) ])# & (table_minima["Amplitude"]<0.0016) ])#0.012
dark_count_error = np.sqrt(dark_count_rate)
print "Dark counts = ",dark_count_rate
run_time = table_minima["Timestamp"].iloc[len(table_minima)-1]
run_time -= bursts_time
print "Total run time = ",run_time
dark_count_rate /= run_time
dark_count_error /= run_time
print "Dark count rate = ",dark_count_rate,"+/-",dark_count_error,"\n"

# now that we are done with calculation, it's time to plot!

#####################
## Dark count plot ##
#####################
# we need to create a dummy plot to get the binning
fig_dummy,axes_dummy = plt.subplots()
bin_content, bins, _ = axes_dummy.hist(table_minima["DeltaT"],alpha=0,bins=100)
plt.close(fig_dummy)

# now we do the split plot for the dark count
fig,axes = plt.subplots(2, sharex=True, gridspec_kw={'height_ratios': [2, 1]})
plt.xscale("log")
figure_title = sipm_name + "  -  OV = " + str(ov) + "V"
axes[0].set_title(figure_title)
axes[0].scatter(table_minima["DeltaT"],table_minima["Amplitude"],marker="o",facecolors="none",edgecolors="black")
axes[0].set_xlim([1e-8,30])
axes[0].set_ylim([0,0.01])
entries = r"$N_{entries}$ = " + str(len(table_minima))
axes[0].text(0.02,0.88,entries,horizontalalignment="left",transform=axes[0].transAxes)
axes[0].set_ylabel("Amplitude (V)")

log_bins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
bin_content, _, _ = axes[1].hist(table_minima["DeltaT"],bins=log_bins,histtype="step")
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

#plt.show()
fig_name = "Amplitude_vs_Dt_" + sipm_name + "_OV" + ov + ".pdf"
fig.savefig(fig_name)


#############################
## Timestamp & DeltaT plot ##
#############################
fig, axes = plt.subplots(2, sharex=False, gridspec_kw={'height_ratios': [1, 1]})
figure_title = sipm_name + "  -  OV = " + str(ov) + "V"
axes[0].set_title(figure_title)
axes[0].set_yscale("log")
axes[0].set_ylabel("counts")
axes[0].set_xlabel("Timestamp (s)")
axes[0].hist(table_minima["Timestamp"],histtype="step",bins=1000)

axes[1].set_xlabel(r"$N_{event}$")
axes[1].set_ylabel(r"$\Delta$t (s)")
table_minima = table_minima.reset_index()
axes[1].step(table_minima.index,table_minima["DeltaT"])#,marker=".")#,facecolors="none",edgecolors="black")

plt.tight_layout()
#plt.show()
fig_name = "Timestamp_Dt_" + sipm_name + "_OV" + ov + ".pdf"
fig.savefig(fig_name)


###################
## Amplitue plot ##
###################
fig, axes = plt.subplots()
figure_title = sipm_name + "  -  OV = " + str(ov) + "V"
axes.set_title(figure_title)
axes.set_yscale("log")
axes.set_ylabel("counts")
axes.set_xlabel("Amplitude (V)")
# gets the bin content and center to find maxima..
bin_content, bin_edge, _ = axes.hist(table_minima["Amplitude"],histtype="step",bins=10000)
# calculate bin center from bin edge
bin_width = bin_edge[1]-bin_edge[0]
bin_center = bin_edge+bin_width/2.
# adjust length (bin center was taken from edges so it's n+1)
bin_center = bin_center[:-1]
bin_center_dictionary = {}
for idx in range(len(bin_center)):
    bin_center_dictionary[bin_center[idx]] = idx

def find_hist_peaks(bin_center, bin_content):
    if max(bin_content) < 2:
        print max_bin_center
        return
        #return max_bin_list

    max_bin = bin_center[bin_content==max(bin_content)][0]
    max_bin_center.append(max_bin)
    print max_bin
    print max_bin_center

    start_searching_from = max_bin + max_bin_center[0]/2.
    new_bin_center = bin_center[bin_center>start_searching_from]
    new_bin_content = bin_content[bin_center>start_searching_from]

    find_hist_peaks(new_bin_center, new_bin_content)

max_bin_center = []
find_hist_peaks(bin_center, bin_content)
print max_bin_center
max_list = [bin_center_dictionary[bc] for bc in max_bin_center]
print max_list

#sys.exit()
#.. then find those maxima, excluding bins where the content is 1..
# (max list is a list of indexes showing where the maxima are in bin_content)
#max_list = argrelextrema(bin_content[:len(bin_content)-10], np.greater_equal, order=30)[0]
# max_list, _ = scipy.signal.find_peaks(bin_content[:len(bin_content)-10],prominence=3)
# print max_list
# max_list = max_list[bin_content[max_list]>1]
# print max_list
#.. and finaly plots the peaks as a scatter plot
axes.scatter(bin_center[max_list],bin_content[max_list], marker='o', s=100, c='g')
# now get y axis range before plotting fit on top
my_ylim = axes.get_ylim()

# define here our gaussian function
def gaussian(x, amplitude, mean, stdev):
    return amplitude / stdev * np.sqrt(2*np.pi) * np.exp( - ( (x - mean) / (stdev*np.sqrt(2)) )**2 )

# now loop over the maxima and fit them
for max_index in max_list:
    print max_index
    # define a range of bins to fit for every maximum
    bin_range = [max(max_index-20,0),min(max_index+20,len(bin_center)-1)]
    fit_bin_range = bin_center[bin_range[0]:bin_range[1]]
    print fit_bin_range
    fit_bin_values = bin_content[bin_range[0]:bin_range[1]]
    # then fit them and plot the fit
    #fit_bin_sigmas = np.where(np.sqrt(fit_bin_values)==0,1,np.sqrt(fit_bin_values))
    # popt, _ = optimize.curve_fit(gaussian, fit_bin_range, fit_bin_values, sigma=fit_bin_sigmas, maxfev=10000)
    try:
        popt, _ = optimize.curve_fit(gaussian, fit_bin_range, fit_bin_values, maxfev=10000)
        axes.plot(fit_bin_range, gaussian(fit_bin_range, *popt), linewidth=1, color="orange")
    except: print "Failed at fitting peak on bin",max_index,"!"
# here it reapplies the range found before the gaussian fits were made
axes.set_ylim(my_ylim)
# and finally shows the plot
plt.tight_layout()
plt.show()
