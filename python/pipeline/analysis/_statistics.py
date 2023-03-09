import BioSimSpace as BSS
import os
import itertools as it
import sys
import re
import glob
import csv
import numpy as np
import math
import pandas as pd
import networkx as nx
import yaml
import pickle
import tempfile
import itertools as it
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import csv
import numpy as np
import pandas as pd 
from rdkit import Chem

from ..utils import *
from ._network import *
from ._analysis import *
from ._plotting import *
from ._convergence import *
from ._dictionaries import *

import cinnabar
from cinnabar import wrangle,plotting,stats

class stats_engines(plotting_engines):


    def __init__(self, analysis_object=None, output_folder=None):
        # inherit the init from other protocol too
        super().__init__(analysis_object, output_folder)

        self._set_statistic_dicts()
    
    def _set_statistic_dicts(self):

        self.statistics  = ['RMSE', 'MUE', 'R2', 'rho','RAE','KTAU']

        # for statistics compared to experimental value
        self.statistics_dict_exper = {}
        for statistic in self.statistics:
            self.statistics_dict_exper[statistic] = None
        
        # TODO also set MAE plain


    def calc_mae(self, pert_val=None, engines=None):
        # calc mae for a provided dictionary in the format wanted

        pv = validate.pert_val(pert_val)

        values_dict = self.values_dict
        if engines:
            engines = validate.engines(engines)
        else:
            engines = self.engines

        mae_pert_df = pd.DataFrame(columns=engines,index=engines)
        mae_pert_df_err = pd.DataFrame(columns=engines,index=engines)

        # iterate over all possible combinations
        for combo in it.product(engines, engines):
            eng1 = combo[0]
            eng2 = combo[1]

            eng1_vals = []
            eng2_vals = []

            #TODO also need to fix this to contain the correct matched values!!
            # first create df of values
            # make sure the values exist!
            for pert in values_dict[eng1][f"{pv}_results"]:
                if pert in values_dict[eng2][f"{pv}_results"]:
                    if values_dict[eng1][f"{pv}_results"][pert][0] != None:
                        if values_dict[eng2][f"{pv}_results"][pert][0] != None:
                            eng1_vals.append(values_dict[eng1][f"{pv}_results"][pert][0])
                            eng2_vals.append(values_dict[eng2][f"{pv}_results"][pert][0])

            # double check only values w corresponding values were used
            if len(eng1_vals) == len(eng2_vals):
                mean_absolute_error = mae(eng1_vals,eng2_vals)  
                data_for_df = {"eng1":eng1_vals,"eng2":eng2_vals}
                data_df= pd.DataFrame(data_for_df)
            else:
                print("cant calc")

            boots = []
            n_boots = 10000

            for n in range(n_boots):
                sample_df = data_df.sample(n=len(eng1_vals), replace=True)
                mae_sample = (abs(sample_df['eng1'] - sample_df['eng2']).sum())/len(eng1_vals)
                boots.append(mae_sample)
            
            mae_err = (np.std(boots))

            mae_pert_df.loc[eng1,eng2]=mean_absolute_error
            mae_pert_df_err.loc[eng1,eng2]=mae_err

        mae_pert_df.to_csv(f"{self.output_folder}/mae_pert_{self.file_ext}.csv", sep=" ")
        mae_pert_df_err.to_csv(f"{self.output_folder}/mae_pert_err_{self.file_ext}.csv", sep=" ")

        return mae_pert_df, mae_pert_df_err

    def _get_x_y(self, pert_val=None, data_x=None, data_y=None, x=None, y=None, xerr=None, yerr=None):

        # default is getting compared to experimental
        # this should be so can match any combo of a and b based on what is available in values dict

        # have dict to df in plotting and also have match engine to other results - dont wanna use dfs!
                    
        pv = validate.pert_val(pert_val)

        if data_x and data_y:

            data_x = self.values_dict[data_x][f"{pv}_results"]
            data_y = self.values_dict[data_x][f"{pv}_results"]

            # need the same length of the data
            # check that 
            len_x = len(data_x.values())
            len_y = len(data_y.values())
            print(len_x)
            print(len_y)
            if len_x != len_y:
                print("data points missing, discarding some for stats analysis...")
                if len_x < len_y:
                    data_used = data_x
                elif len_y < len_x:
                    data_used = data_y

                for _data in [data_x, data_y]:
                    for key in _data.keys():
                        if key not in data_used.keys():
                            del _data[key]
            
            # make sure all the values are paired together
            for value in data_x:
                print(value)
            #     try:
            #         exp_ddG = self.values_dict["experimental"][f"{pv}_results"][value][0]
            #         exp_err = self.values_dict["experimental"][f"{pv}_results"][value][1]
            #         comp_ddG = self.values_dict[eng][f"{pv}_results"][value][0]
            #         comp_err = self.values_dict[eng][f"{pv}_results"][value][1]
            #         freenrg_pert_dict[value] = [exp_ddG, exp_err, comp_ddG, comp_err]
            #     except Exception as e:
            #         # print(e)
            #         print(f"could not convert analysis object {value}, {eng}, {pv}, into dataframe. Was it able to be computed earlier?")
            # freenrg_df = pd.DataFrame(freenrg_pert_dict, index=["freenrg_exp", "err_exp", "freenrg_calc", "err_calc"]).transpose()



            x = [val[0] for val in self.values_dict[data_x][f"{pv}_results"].values()]
            y = [val[0] for val in self.values_dict[data_y][f"{pv}_results"].values()]
            xerr = np.asarray([val[1] for val in self.values_dict[data_x][f"{pv}_results"].values()])
            yerr = np.asarray([val[1] for val in self.values_dict[data_y][f"{pv}_results"].values()])
  
        else:
            try:
                x = x
                y = y
                xerr = xerr
                yerr = yerr
            except:
                print("if not providing data_x and data_y (which should be a name in the self.names_list),\
                      please provide x,y,xerr,yerr values")
        
        return x,y,xerr,yerr


    def _compute_stats(self, pert_val=None, data_x=None, data_y=None, statistic=None, x=None, y=None, xerr=None, yerr=None):

        if statistic not in self.statistics:
            raise ValueError(f"please use one of the statistics in {self.statistics}")

        # get the x y values from the dictionaries, also validates pert val and engine
        x,y,xerr,yerr = self._get_x_y(pert_val, data_x, data_y, x, y, xerr, yerr)

        # using cinnabar function
        s = stats.bootstrap_statistic(x, y, xerr, yerr, statistic=statistic)
        values = (s['mle'], s['low'], s['high'])
        # string = f"{statistic}:   {s['mle']:.2f} [95%: {s['low']:.2f}, {s['high']:.2f}] " + "\n"
            
        return values

    def compute_statistics(self):

        for pv in ["pert","val"]:
            self.compute_mue(pv, self.engines)
            # all compute functions

    def compute_mue(self, pert_val=None, engines=None):

        # validate from other that it is in names list

        pert_val = validate.pert_val(pert_val)
        if engines:
            engines = validate.engines(engines)
        else:
            engines = self.engines

        for eng in engines:
            values = self._compute_stats(pert_val, "experimental", eng, statistic="MUE")
            print(values)
            self.statistics_dict_exper["MUE"][eng] = values

