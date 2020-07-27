import sys
from peak_class import *

# print usage
if len(sys.argv) != 5:
    print("\nUSAGE: python peak_analyzer.py file_path simpm_name ov exclude_bursts")
    print("\nARGUMENTS:")
    print("file_path:\t\t path to .csv file with timestamps and amplitude of peaks (e.g. \"R00029_Dark/OV2/Results.csv\")")
    print("sipm_name:\t\t SiPM name, used for plot title (e.g. \"R00029\")")
    print("ov:\t\t\t over voltage value, used for plot title (2-5)")
    print("exclude_bursts:\t\t 0 - plot all events; 1 - exclude event bursts after saturation; 2 - exclude all bursts\n")
    sys.exit()
print("\n----------  PEAK ANALYZER LAUNCHED  ----------\n")

# from user input, get path to file to access table of peaks, and sipm name and ov for plot titles
file_path = sys.argv[1] # path to file with timestamps and amplitude of peaks (e.g. "R00029_Dark/OV2/Results.txt")
sipm_name = sys.argv[2] # SiPM name, used for pplot title (e.g. "R00029")
ov = sys.argv[3] # over voltage value (2-5)

# initialize the class, and define the burst time
_peak_ = Peak(sys.argv[1],False)
n_bursts = 0
bursts_time = 0

#calculate run time (to compute rates)
run_time = _peak_.table_minima["Timestamp"].iloc[len(_peak_.table_minima)-1]

# check if we have to exclude bursts from the table..
n_min_burst_ev = 20
max_hreshold = 0.9

if (int)(sys.argv[4]) == 0:
    print("Event bursts will not be removed")
elif (int)(sys.argv[4]) == 1:
    print("Removing event bursts after a saturating event")
    n_bursts, bursts_time = _peak_.exclude_bursts(True,n_min_burst_ev)
    print("Number of bursts:",n_bursts,"; burst rate",n_bursts/run_time)
    print("Total integrated time of bursts:",bursts_time,"\n")
elif (int)(sys.argv[4]) == 2:
    print("Removing all event bursts")
    n_bursts, bursts_time = _peak_.exclude_bursts(False,n_min_burst_ev)
    print("Number of bursts:",n_bursts,"; burst rate",n_bursts/run_time)
    print("Total integrated time of bursts:",bursts_time,"\n")
else:
    print("ERROR: invalid value for exclude_bursts parameter, chose 0, 1, or 2")
    sys.exit()

# Now, time to plot!
# First we plot amplitude, and obtain the value of the fist pe in V
ampl_n_bins = 400
first_pe =_peak_.plot_amplitude(ampl_n_bins,False,False,sipm_name,ov,True,False)
print("First photoelectron peak found at",first_pe,"V")
first_pe = first_pe*1.35
after_pulse_end = 1e-6
# now use such value to calculate the dark count and cross talk rates
# definition of dark count based on discussion during the meetings
dark_count_rate = len(_peak_.table_minima[ (_peak_.table_minima["DeltaT"]>after_pulse_end) ])
dark_count_error = np.sqrt(dark_count_rate)
cross_talk_rate = len(_peak_.table_minima[ (_peak_.table_minima["Amplitude"]>first_pe) ])
cross_talk_error = np.sqrt(cross_talk_rate)
after_pulse_rate = len(_peak_.table_minima[ (_peak_.table_minima["DeltaT"]<after_pulse_end) & (_peak_.table_minima["Amplitude"]<first_pe) ])
after_pulse_error = np.sqrt(after_pulse_rate)
print("Dark counts =",dark_count_rate)
print("Cross talk events =",cross_talk_rate)
print("After pulses =",after_pulse_rate)
run_time -= bursts_time
print("Total run time = ",run_time)
dark_count_rate /= run_time
dark_count_error /= run_time
cross_talk_rate /= run_time
cross_talk_error /= run_time
after_pulse_rate /= run_time
after_pulse_error /= run_time
print("Dark count rate =",dark_count_rate,"+/-",dark_count_error)
print("Cross talk rate =",cross_talk_rate,"+/-",cross_talk_error)
print("After pulse rate =",after_pulse_rate,"+/-",after_pulse_error,"\n")

# now that we are done with calculation, it's time to plot the DCR and time
ampl_bin_range = [0,0.05]
_peak_.plot_dark_count(ampl_bin_range,True,False,sipm_name,ov,True,first_pe,after_pulse_end)

time_n_bins = 1000
_peak_.plot_times(time_n_bins,True,False,sipm_name,ov)

# finnaly, we show all plots
#plt.show()
