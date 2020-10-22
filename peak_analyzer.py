import sys
from peak_class import *

# print usage
if len(sys.argv) != 5 and len(sys.argv)!=8 and len(sys.argv)!=9:
    print("\nUSAGE: python peak_analyzer.py file_path simpm_name ov exclude_bursts (n_min_burst_ev max_burst_delay min_burst_frequency (max_threshold))")
    print("\nARGUMENTS:")
    print("file_path:\t\t path to .csv file with timestamps and amplitude of peaks (e.g. \"R00029_Dark/OV2/Results.csv\")")
    print("sipm_name:\t\t SiPM name, used for plot title (e.g. \"R00029\")")
    print("ov:\t\t\t over voltage value, used for plot title (2-5)")
    print("exclude_bursts:\t\t 0 - plot all events; 1 - exclude event bursts after saturation; 2 - exclude all bursts")
    print("(n_min_burst_ev):\t if exclude_bursts == 1 or 2 - can pecify minimum number of events in a burst")
    print("(max_burst_delay):\t if exclude_bursts == 1 or 2 - can specify maximum time delay between consecutive events in a burst (default 20)")
    print("(min_burst_frequency):\t if exclude_bursts == 1 or 2 - can specify minimum frequencey of events in a burst (default 0.2)")
    print("(max_threshold):\t if exclude_bursts == 1 - can specify what fraction of the maximum amplitude is considered a saturating event (default 0.8)")

    sys.exit()
print("\n----------  PEAK ANALYZER LAUNCHED  ----------\n")

# from user input, get path to file to access table of peaks, and sipm name and ov for plot titles
file_path = sys.argv[1] # path to file with timestamps and amplitude of peaks (e.g. "R00029_Dark/OV2/Results.txt")
sipm_name = sys.argv[2] # SiPM name, used for pplot title (e.g. "R00029")
ov = sys.argv[3] # over voltage value (2-5)

# from user input with default values, get parameter of burst search
n_min_burst_ev = 20 if len(sys.argv)!=8 and len(sys.argv)!=9 else int(sys.argv[5])
max_burst_delay = 0.2 if len(sys.argv)!=8 and len(sys.argv)!=9 else float(sys.argv[6])
min_burst_delay = 4e-6 # to avoid counting after pulses
min_burst_frequency = 5.0 if len(sys.argv)!=8 and len(sys.argv)!=9 else float(sys.argv[7])
max_threshold = 0.8 if len(sys.argv)!=9 else float(sys.argv[8])

# initialize the class, and define the burst time
_peak_ = Peak(sys.argv[1],False)
n_bursts = 0
bursts_time = 0

#calculate run time (to compute rates)
run_time = _peak_.table_minima["Timestamp"].iloc[len(_peak_.table_minima)-1]

# check if we have to exclude bursts from the table..
if (int)(sys.argv[4]) == 0:
    print("Event bursts will not be removed")
elif (int)(sys.argv[4]) == 1:
    print("Removing event bursts after a saturating event")
    n_bursts, bursts_time = _peak_.exclude_bursts(n_min_burst_ev,max_burst_delay,min_burst_delay,min_burst_frequency,True,max_threshold)
    print("Number of bursts:",n_bursts,"; burst rate",n_bursts/run_time)
    print("Total integrated time of bursts:",bursts_time,"\n")
elif (int)(sys.argv[4]) == 2:
    print("Removing all event bursts")
    n_bursts, bursts_time = _peak_.exclude_bursts(n_min_burst_ev,max_burst_delay,min_burst_delay,min_burst_frequency,False,max_threshold)
    print("Number of bursts:",n_bursts,"; burst rate",n_bursts/run_time)
    print("Total integrated time of bursts:",bursts_time,"\n")
else:
    print("ERROR: invalid value for exclude_bursts parameter, chose 0, 1, or 2")
    sys.exit()

# Now, time to plot!
# First we plot amplitude, and obtain the value of the fist pe in V
ampl_n_bins = 400
n_ampl_peaks = 1
first_pe =_peak_.plot_amplitude(ampl_n_bins,n_ampl_peaks,False,False,sipm_name,ov,True,True,False)
print("First photoelectron peak found at",first_pe,"V\n")
first_pe = first_pe*1.35
after_pulse_end = 4e-6
# now use such value to calculate the dark count and cross talk rates
# definition of dark count based on discussion during the meetings
dark_count_rate = len(_peak_.table_minima[ (_peak_.table_minima["DeltaT"]>after_pulse_end) ])
dark_count_error = np.sqrt(dark_count_rate)
cross_talk_rate = len(_peak_.table_minima[ (_peak_.table_minima["Amplitude"]>first_pe) ])
cross_talk_error = np.sqrt(cross_talk_rate)
after_pulse_rate = len(_peak_.table_minima[ (_peak_.table_minima["DeltaT"]<after_pulse_end) & (_peak_.table_minima["Amplitude"]<first_pe) ])
total_rate=len(_peak_.table_minima)
total_rate_error=np.sqrt(total_rate)
after_pulse_error = np.sqrt(after_pulse_rate)
print("Dark counts =",dark_count_rate)
print("Cross talk events =",cross_talk_rate)
print("After pulses =",after_pulse_rate)
run_time -= bursts_time
print("Total run time = ",run_time)
total_rate/=run_time
total_rate_error/=run_time
dark_count_error /= run_time
dark_count_rate /= run_time
cross_talk_rate /= run_time
cross_talk_error /= run_time
after_pulse_rate /= run_time
after_pulse_error /= run_time
print("Total count rate =",total_rate,"+/-",total_rate_error)
print("Dark count rate =",dark_count_rate,"+/-",dark_count_error)
print("Cross talk rate =",cross_talk_rate,"+/-",cross_talk_error)
print("After pulse rate =",after_pulse_rate,"+/-",after_pulse_error,"\n")

# now that we are done with calculation, it's time to plot the DCR and time
ampl_bin_range = [0,0.035]
_peak_.plot_dark_count(ampl_bin_range,False,False,sipm_name,ov,True,first_pe,after_pulse_end)

bins_per_sec = 1
_peak_.plot_times(bins_per_sec,False,False,sipm_name,ov)

# finnaly, we show all plots
plt.show()
