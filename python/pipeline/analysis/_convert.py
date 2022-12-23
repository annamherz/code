
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

class convert:
    """class of static methods for converting data 
    """

    def __init__(self):
        pass

    # TODO more robust yml file conversion
    @staticmethod
    def yml_into_freenrgworkflows(exp_file, exp_file_dat, data_format=None):
        # get the experimental data into a useable format (from yml to csv)
        # for freenergworkflows, want to save as lig, Ki
        # experimental values (e.g. ic50/ki) for all ligands in our set.

        exp_file = validate.file_path(exp_file)
        exp_file_dat = validate.string(exp_file_dat)

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

    @staticmethod
    def convert_M_kcal(value, magnitude = "uM"):
        temp = 300
        if magnitude == "uM":
            power = 10**-6
        # gas constant in kcal per Kelvin per mol, exp val converted into M
        kcal_val = 0.0019872041*temp*np.log(value*(power))

        return kcal_val

    @staticmethod
    def yml_into_exper_dict(exp_file, exp_file_dat, data_format=None):
        
        # format of the data from yml into freenrgworkflows is uM
        convert.yml_into_freenrgworkflows(exp_file, exp_file_dat)

        exper_val_dict = {}

        with open(exp_file_dat, "r") as file:
            for line in file:
                if not "ligand" in line: # dont use header
                    lig = line.split(",")[0]

                    exp_val = float(line.split(",")[1].strip())
                    # convert into kcal mol
                    exp_kcal = convert.convert_M_kcal(exp_val)

                    # convert both upper and lower error bounds for this too
                    # get average and keep this as the error
                    err = float(line.split(",")[2].strip())
                    exp_upper = exp_val + err
                    exp_lower = exp_val - err
                    exp_upper_kcal = convert.convert_M_kcal(exp_upper)
                    exp_lower_kcal = convert.convert_M_kcal(exp_lower)
                    err_kcal = abs(exp_upper_kcal - exp_lower_kcal)/2

                    # add to dict
                    exper_val_dict.update({lig:(exp_kcal, err_kcal)})            

        return exper_val_dict

    @staticmethod
    def cinnabar_file(results_files, exper_val, output_file):
        # files is a list of files
        # output file

        # first check if the experimental values are a dict or a file
        try:
            exper_val_dict = validate.dictionary(exper_val)
        except:
            validate.file_path(exper_val)
            is_file = True
        
        if is_file:
            print("input is a file, will convert this into a dict...")
            print("please check that the conversion of values is okay.")
            exper_val_dict = convert.yml_into_exper_dict(exper_val)

        # write to a csv file
        with open(f"{output_file}.csv", "w") as cinnabar_data_file:
            writer = csv.writer(cinnabar_data_file, delimiter=",")

            # first, write the experimental data
            writer.writerow(["# Experimental block"])
            writer.writerow(["# Ligand","expt_DDG","expt_dDDG"])


            # TODO calc exp from yml, follwed by conversion into dict
            for lig in exper_val_dict.keys():
                writer.writerow([lig,f"{exper_val[lig][0]}",f"{exper_val[lig][1]}"])


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
