import sys
from waveform_class import *

#creates an object waveform of the waveform class
#_waveform_ = Waveform(sys.argv[1],sys.argv[2],sys.argv[5] if len(sys.argv)==6 else 6250)
_waveform_ = Waveform(sys.argv[1],sys.argv[2],(int)(sys.argv[3]))

#creates a dictionary of command options to be passed to the class methods
def wf_dict():
   config_file = open(sys.argv[4])

   lines = config_file.readlines()
   mydict = {}
   for line in lines:
       if line[0] == "-": continue
       var,val = line.split('=')
       mydict[var] = val.strip()
       if mydict[var]=="True" or mydict[var]=="true" or mydict[var]==1: mydict[var] = True
       elif mydict[var]=="False" or mydict[var]=="false" or mydict[var]==0: mydict[var] = False

       if mydict[var]=="first_points": mydict[var] = 0
       if mydict[var]=="max_range": mydict[var] = 1

       if mydict[var]=="previous_points": mydict[var] = 0
       if mydict[var]=="baseline_gap": mydict[var] = 1

   return mydict

#creates on object dictionary
_wf_dict_ = wf_dict()


# ---- Baseline ----
# baseline_method = "first_point"
# baseline_n_points = 100
# ---- Minima ----
# minimum_search_range=50 OK
# minimum_back_shift=100
# minimum_n_points=100
# minimum_gap=0.03


#runs the algorithm
n_ev = sys.argv[6] if (len(sys.argv)==7 and sys.argv[6]>0) else 1000000
_waveform_.find_all_minima(n_ev,_wf_dict_["baseline_method"],(int)(_wf_dict_["baseline_n_points"]),(int)(_wf_dict_["minimum_method"]),(int)(_wf_dict_["minimum_search_range"]),(int)(_wf_dict_["minimum_back_shift"]),(int)(_wf_dict_["minimum_n_points"]),(float)(_wf_dict_["minimum_gap"]),(int)(_wf_dict_["minimum_n_close"]),_wf_dict_["show_plot"],_wf_dict_["save_plot"],sys.argv[5])

_waveform_.save_minimum_table(sys.argv[5])
