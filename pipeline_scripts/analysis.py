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

# someway that this is for the transformations listed in the csv file made at the start for the setup
# or alternatively do this in bash
trans = sys.argv[1].rstrip()
lig_1 = trans.split('~')[0]
lig_2 = trans.split('~')[1]
engine = sys.argv[2].rstrip()

# options
extra_options = {'save_graphs': True, 'save_pickle':True}
chosen_estimator = "MBAR" # MBAR or TI
chosen_method = "alchemlyb" # native or alchemlyb
est_meth = f"{chosen_estimator}_{chosen_method}"
method = f"benchmark" # this is for writing in the output folder
file_ext = f"{engine}_{est_meth}_{method}" # for all files written

main_dir = _os.environ["MAINDIRECTORY"]
# main_dir = "/home/anna/Documents/benchmark/test_tyk2_benchmark_sage"
final_results_file_path = f"{main_dir}/outputs/final_summary_{file_ext}.csv"

# find correct path, use extracted if it exists
if _os.path.exists(f"{main_dir}/outputs/{engine}_extracted/{trans}"):
    path_to_dir = f"{main_dir}/outputs/{engine}_extracted/{trans}"
else:
    path_to_dir = f"{main_dir}/outputs/{engine}/{trans}"

if not _os.path.exists(path_to_dir):
    raise OSError(f"{path_to_dir} does not exist.")

print(f'analysing results for {path_to_dir}')
print(f"using {chosen_method} and {chosen_estimator} for analysis")

analysis = analyse(path_to_dir)
analysis.set_options(extra_options)

avg, error, repeats_tuple_list = analysis.analyse_all_repeats()

# ####### WRITING DATA for the final result
# data point for average
data_point_avg = [lig_1, lig_2, str(avg), str(error), engine, est_meth, method]
print(data_point_avg)

# use csv to open the results file.
with open(final_results_file_path, "a") as freenrg_writefile:
    writer = csv.writer(freenrg_writefile)

    # first, write a header if the file is created for the first time.
    if _os.path.getsize(final_results_file_path) == 0:
        print(f"Starting {final_results_file_path} file.")
        writer.writerow(["lig_1", "lig_2", "freenrg",
                        "SEM", "engine", "estimator", "method"])


with open(final_results_file_path, "r") as freenrg_readfile:
    # then, grab all of the data that is already in the file.
    reader = csv.reader(freenrg_readfile)
    data_entries = [row for row in reader]

# check if our data entry is not already in the results file. Raise an error if is.
if data_point_avg in data_entries:
    warnings.warn(
        f"Results for in {trans}, {engine} are already in {final_results_file_path}.")

else:
    # at this point we know that we are writing a new entry in the results file. Append the line to the file.
    # use csv to open the results file.
    with open(final_results_file_path, "a") as freenrg_writefile:
        writer = csv.writer(freenrg_writefile)
        print(
            f"Writing results. Average free energy of binding is {avg} and the SEM is {error} for {trans}, {engine}.")
        writer.writerow(data_point_avg)


# write results for each repeat also
no_repeats = list(range(len(repeats_tuple_list)))
# use csv to open the results file.
for r in no_repeats:
    data_point = [lig_1, lig_2, repeats_tuple_list[r][1], repeats_tuple_list[r][2], engine, est_meth, method]
    results_file_path = f"{main_dir}/outputs/repeat_{no_repeats.index(r)}_{file_ext}.csv"
    with open(results_file_path, "a") as freenrg_writefile:
        writer = csv.writer(freenrg_writefile)

        # first, write a header if the file is created for the first time.
        if _os.path.getsize(results_file_path) == 0:
            print(f"Starting {results_file_path} file.")
            writer.writerow(["lig_1", "lig_2", "freenrg",
                            "error", "engine", "estimator", "method"])


    with open(results_file_path, "r") as freenrg_readfile:
        # then, grab all of the data that is already in the file.
        reader = csv.reader(freenrg_readfile)
        data_entries = [row for row in reader]

    # check if our data entry is not already in the results file. Raise an error if is.
    if data_point in data_entries:
        warnings.warn(
            f"Results for in {trans}, {engine} are already in {results_file_path}.")

    else:
        # at this point we know that we are writing a new entry in the results file. Append the line to the file.
        # use csv to open the results file.
        with open(results_file_path, "a") as freenrg_writefile:
            writer = csv.writer(freenrg_writefile)
            print(
                f"Writing results. For repeat {r}, free energy of binding is {repeats_tuple_list[r][1]} and the error is {repeats_tuple_list[r][2]} for {trans}, {engine}.")
            writer.writerow(data_point)
