
# import libraries
import BioSimSpace as BSS
import os
import glob
import csv
import numpy as np
import math
import pandas as pd
import networkx as nx
import yaml
from scipy.stats import sem as sem
from scipy.stats import bootstrap
from sklearn.metrics import mean_absolute_error as mae
import pickle
import tempfile
import itertools

import matplotlib.pyplot as plt
import csv
import numpy as np
import pandas as pd 

from ..utils import *


def convert_yml_into_freenrgworkflows(exp_file, exp_file_dat):
    # get the experimental data into a useable format (from yml to csv)
    # for freenergworkflows, want to save as lig, Ki
    # experimental values (e.g. ic50/ki) for all ligands in our set.

    exp_file = validate.file_path(exp_file)

    with open(exp_file, "r") as file:
        data = yaml.safe_load(file) # loads as dictionary

    with open(exp_file_dat, "w") as file:
        writer = csv.writer(file, delimiter=",")
        writer.writerow(["ligand","value","error"])

        # the data needs to be IC50, uM
        # am assuming that ki and IC50 are the same
        
        for key in data.keys(): # write for each ligand that was in yaml file
            if data[key]['measurement']['unit'] == 'uM':
                writer.writerow([key, data[key]['measurement']['value'], data[key]['measurement']['error']])
            elif data[key]['measurement']['unit'] == 'nM':
                writer.writerow([key, "{:.4f}".format(data[key]['measurement']['value']/1000), data[key]['measurement']['error']/1000])



def convert_cinnabar_file(results_files, exper_val_dict, output_file):
    # files is a list of files
    # output file

    # write to a csv file
    with open(f"{output_file}.csv", "w") as cinnabar_data_file:
        writer = csv.writer(cinnabar_data_file, delimiter=",")

        # first, write the experimental data
        writer.writerow(["# Experimental block"])
        writer.writerow(["# Ligand","expt_DDG","expt_dDDG"])


        # TODO write function to convert experimental values instead of freenergworkflows (take from other ana)
        for lig in exper_val_dict.keys():
            writer.writerow([lig,f"{exper_val_dict[lig][0]}",f"{exper_val_dict[lig][1]}"])


        # second write the perturbation data
        writer.writerow([" "])
        writer.writerow(["# Calculated block"])
        writer.writerow(["# Ligand1","Ligand2","calc_DDG","calc_dDDG(MBAR)", "calc_dDDG(additional)"])

        for file in results_files:
            with open(file, "r") as res_file:
                for line in res_file:
                    if "freenrg" in line: # avoid the header
                        pass
                    else:                           # write each perturbation and repeat to the file
                        lig_0 = line.split(",")[0]
                        lig_1 = line.split(",")[1]
                        comp_ddG = line.split(",")[2]
                        comp_err = line.split(",")[3]
            
                        writer.writerow([lig_0, lig_1, comp_ddG, comp_err, "0.0"])
