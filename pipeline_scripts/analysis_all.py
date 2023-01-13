# analysis script for plotting and stats for all the diff engines
# anna

from inspect import BoundArguments
import warnings
import BioSimSpace as BSS
from BioSimSpace import _Exceptions
from BioSimSpace import Units as _Units
import sys
import csv
import pickle
import os as _os
import numpy as _np
import pandas as _pd
import math as _math
import itertools as it
from scipy.stats import sem

try:
    import pipeline
except:
    print("adding code to the pythonpath...")
    code = '/home/anna/Documents/code/python'
    if code not in sys.path:
        sys.path.insert(1, code)
    import pipeline

from pipeline.analysis import *


# location of repeat files
# results folder for output
main_dir = _os.environ["MAINDIRECTORY"]

if _os.path.exists(f"{main_dir}/outputs/results"):
    final_results_folder = f"{main_dir}/outputs/results"
elif _os.path.exists(f"{main_dir}/outputs_extracted/results"):
    final_results_folder = f"{main_dir}/outputs_extracted/results"
else:
    raise ValueError(f"results directory not found in the {main_dir}")

# options that should be done in dict form
# options
analysis_options = {"engines": "all",
                    "bar":True,
                    "scatter":True,
                    "cycle_closure": True,
                    "convergence": True,
                    "XXX":"XXX",
                    "try_pickle":True,'save_pickle':True}


all_analysis_object = analysis_engine(final_results_folder)
all_analysis_object.set_options(analysis_options)
all_analysis_object.set_experimental() # add experimental files
# add any other results files

all_analysis_object.plot_graphs() 

# plot all graphs plots all the below if it is in the options
# plotting individually will always set it to true for plotting?
all_analysis_object.plot_bar_pert_single()
all_analysis_object.plot_bar_pert_all()
all_analysis_object.plot_bar_lig_single()
all_analysis_object.plot_bar_lig_all()
all_analysis_object.plot_scatter_pert_single()
all_analysis_object.plot_scatter_pert_all()
all_analysis_object.plot_scatter_lig_single() # cinnabar
all_analysis_object.plot_scatter_lig_all()
all_analysis_object.plot_convergence()
all_analysis_object.statistical_analysis()
all_analysis_object.cycle_closure()