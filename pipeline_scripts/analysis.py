# analysis script for all the diff engines
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
from pipeline.utils import write_analysis_file

# someway that this is for the transformations listed in the csv file made at the start for the setup
# or alternatively do this in bash
trans = sys.argv[1].rstrip()
engine = sys.argv[2].rstrip()

# TODO add the different analysis options

# options
analysis_options = {'estimator': "MBAR", "method":"alchemlyb",
                    "check_overlap":True,
                    "try_pickle":True, 'save_pickle':True,
                    "auto_equilibration": False,
                    "truncate_percentage": 0,
                    "truncate_keep":"start"}

main_dir = _os.environ["MAINDIRECTORY"]

# find correct path, use extracted if it exists
if _os.path.exists(f"{main_dir}/outputs/{engine}_extracted/{trans}"):
    path_to_dir = f"{main_dir}/outputs/{engine}_extracted/{trans}"
    final_results_folder = f"{main_dir}/outputs"
elif _os.path.exists(f"{main_dir}/outputs_extracted/{engine}/{trans}"):
    path_to_dir = f"{main_dir}/outputs_extracted/{engine}/{trans}"
    final_results_folder = f"{main_dir}/outputs_extracted"
else:
    path_to_dir = f"{main_dir}/outputs/{engine}/{trans}"
    final_results_folder = f"{main_dir}/outputs"

if not _os.path.exists(path_to_dir):
    raise OSError(f"{path_to_dir} does not exist.")

print(f'analysing results for {path_to_dir}')
print(f"using {analysis_options} for analysis")

# using the pipeline module for analysis
analysis = analyse(path_to_dir)
analysis.set_options(analysis_options)
avg, error, repeats_tuple_list = analysis.analyse_all_repeats()
analysis.plot_graphs()

# write the final result
write_analysis_file(analysis, final_results_folder)
