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



def analyse_all_repeats(work_dir, estimator='MBAR', method="alchemlyb", extra_options=None):
    """Analyse all existing free-energy data from a simulation working directory.

        Parameters
        ----------

        work_dir : str
            The working directory for all the repeats of the simulation.
            The repeats must have either bound and/or free in the name to be recognised.

        estimator : str
            The estimator ('MBAR' or 'TI') used. MBAR is the default.

        extra_options : dict
            Extra options can include the following as keys:
            average_method : 'mean' (default is mean)
            error_method : 'SEM' (default is SEM)
            save_graphs : (bool) - these will be saved in work_dir/graphs

        Returns
        -------

        free_energy : (:class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)
            The average BSS.relative free-energy difference and its associated error.
    """

    if not isinstance(work_dir, str):
        raise TypeError("'work_dir' must be of type 'str'.")
    if not _os.path.isdir(work_dir):
        raise ValueError("'work_dir' doesn't exist!")

    if estimator not in ['MBAR', 'TI']:
        raise ValueError("'estimator' must be either 'MBAR' or 'TI'.")

    # default calculation methods
    error_method = 'SEM'
    average_method = 'mean'
    if estimator == "MBAR":
        check_overlap = True
    else:
        check_overlap = False
    save_graphs = True
    save_pickle = True

    pickle_dir = path_to_dir + '/pickle'
    if not _os.path.exists(pickle_dir):
        _os.mkdir(f"{pickle_dir}")
        extra_options['try_pickle'] =  False
    else:
        extra_options['try_pickle'] =  True

    # if there are extra options, override the default calculation methods
    if extra_options:
        if not isinstance(extra_options, dict):
            # add so it has to be std or smth else
            raise TypeError("'extra_options' must be of type 'dict'.")
        if "error_method" in extra_options:
            if not isinstance(extra_options["error_method"], str):
                raise TypeError(
                    "'error_method' value in the extra_options dictionary must be of type 'str'.")
            else:
                error_method = extra_options["error_method"]
        if "average_method" in extra_options:
            if not isinstance(extra_options["average_method"], str):
                raise TypeError(
                    "'average_method' value in the extra_options dictionary must be of type 'str'.")
            else:
                average_method = extra_options["average_method"]
        if "save_graphs" in extra_options:
            if not isinstance(extra_options["save_graphs"], bool):
                raise TypeError(
                    "'save_graphs' value in the extra_options dictionary must be of type 'bool'.")
            else:
                save_graphs = extra_options["save_graphs"]
        if "save_pickle" in extra_options:
            if not isinstance(extra_options["save_pickle"], bool):
                raise TypeError(
                    "'save_pickle' value in the extra_options dictionary must be of type 'bool'.")
            else:
                save_pickle = extra_options["save_pickle"]
        if "try_pickle" in extra_options:
            if not isinstance(extra_options["try_pickle"], bool):
                raise TypeError(
                    "'try_pickle' value in the extra_options dictionary must be of type 'bool'.")
            else:
                try_pickle = extra_options["try_pickle"]

    # Dictionaries
    bound_pmf_dict = {}  # for the intial results
    free_pmf_dict = {}
    bound_matrix_dict = {}  # for the energy matrix
    free_matrix_dict = {}
    bound_val_dict = {}  # for the actual value for all the windows
    free_val_dict = {}
    bound_err_dict = {}  # for the final error defined
    free_err_dict = {}

    # list for all the repeats
    repeats_tuple_list = []
    # list for the successful calculations
    bound_calculated = []
    free_calculated = []

    # Read how many repeats are in the directory.
    folders = (next(_os.walk(work_dir))[1])
    b_folders = []
    f_folders = []
    for f in folders:
        if 'bound' in f:
            b_folders.append(f'{f}')
        elif 'free' in f:
            f_folders.append(f'{f}')
        else:
            continue
    
    # sort the folders
    b_folders.sort()
    f_folders.sort()

    if not b_folders:
        raise ValueError(
            "Couldn't find any folders with bound in the specified directory?")
    elif not f_folders:
        raise ValueError(
            "Couldn't find any folders with free in the specified directory?")
    else:
        print(f'For bound values, analysing folder(s) {b_folders}')
        print(f'For free values, analysing folder(s) {f_folders}')

    no_of_b_repeats = len(b_folders)
    b_repeats = list(range(no_of_b_repeats))
    no_of_f_repeats = len(f_folders)
    f_repeats = list(range(no_of_f_repeats))

    if no_of_b_repeats != no_of_f_repeats:
        print(
            f"There are a different number of repeats for bound ({no_of_b_repeats}) and free ({no_of_f_repeats}) for {work_dir}.")
    else:
        print(
            f"There are {no_of_b_repeats} repeats for each the bound and the free for {work_dir}.")

    # try loading in if previously calculated
    if try_pickle:
        try:
            print("trying to locate pickles in folder...")
            # TODO fix so pickle name is in function as otherwise cant run independently of this script
            with open(f"{pickle_dir}/bound_pmf_{trans}_{engine}_{estimator}_{method}.pickle", "rb") as file:
                bound_pmf_dict = pickle.load(file)
            with open(f"{pickle_dir}/free_pmf_{trans}_{engine}_{estimator}_{method}.pickle", "rb") as file:
                free_pmf_dict = pickle.load(file)
            print("pickles found!")
        except:
            print("loading pickle failed. Calculating normally.")
            try_pickle = False

    # Analyse the results for each leg of the transformation.
    for b in b_repeats:
        try:
            name = str(b) + '_bound'
            if try_pickle:
                if name in bound_pmf_dict.keys():
                    bound_matrix_dict.update({name: None})
                    bound_calculated.append(name)
                else:
                    print("could not find pickle for that name, rerunning...")
                    pmf_bound, overlap_matrix_bound = BSS.FreeEnergy.Relative.analyse(
                        f'{work_dir}/{b_folders[b]}', estimator=estimator, method=method)
                    bound_pmf_dict.update({name: pmf_bound})
                    bound_matrix_dict.update({name: overlap_matrix_bound})
                    bound_calculated.append(name)                    
            else:
                pmf_bound, overlap_matrix_bound = BSS.FreeEnergy.Relative.analyse(
                    f'{work_dir}/{b_folders[b]}', estimator=estimator, method=method)
                bound_pmf_dict.update({name: pmf_bound})
                bound_matrix_dict.update({name: overlap_matrix_bound})
                bound_calculated.append(name)
        except Exception as e:
            print(e)
            print(
                f'Unable to analyse values for {name}, which is repeat {b_folders[b]} in {work_dir}.')

    for f in f_repeats:
        try:
            name = str(f) + '_free'
            if try_pickle:
                if name in free_pmf_dict.keys():
                    free_matrix_dict.update({name: None})
                    free_calculated.append(name)
                else:
                    pmf_free, overlap_matrix_free = BSS.FreeEnergy.Relative.analyse(
                        f'{work_dir}/{f_folders[f]}', estimator=estimator, method=method)
                    free_pmf_dict.update({name: pmf_free})
                    free_matrix_dict.update({name: overlap_matrix_free})
                    free_calculated.append(name)
            else:
                pmf_free, overlap_matrix_free = BSS.FreeEnergy.Relative.analyse(
                    f'{work_dir}/{f_folders[f]}', estimator=estimator, method=method)
                free_pmf_dict.update({name: pmf_free})
                free_matrix_dict.update({name: overlap_matrix_free})
                free_calculated.append(name)
        except Exception as e:
            print(e)
            print(
                f'Unable to analyse values for {name}, which is repeat {f_folders[f]} in {work_dir}.')

    # create a dictionary of the calculated values, for each bound and free
    for r in b_repeats:
        try:
            bound_name = str(r) + '_bound'
            bound_val = (bound_pmf_dict[bound_name])[-1][1] - \
                (bound_pmf_dict[bound_name])[0][1]
            bound_err = (bound_pmf_dict[bound_name])[-1][2] # TODO change this so as in diff?
            bound_val_dict.update({bound_name: bound_val})
            bound_err_dict.update({bound_name: bound_err})
        except:
            print(f'''Unable to compute values for {bound_name} in {work_dir}.\
                Check earlier error message if these values could be analysed.''')
    for r in f_repeats:
        try:
            free_name = str(r) + '_free'
            free_val = (free_pmf_dict[free_name])[-1][1] - \
                (free_pmf_dict[free_name])[0][1]
            free_err = (free_pmf_dict[free_name])[-1][2] # TODO change this so as in diff?
            free_val_dict.update({free_name: free_val})
            free_err_dict.update({free_name: free_err})
        except:
            print(f'''Unable to compute values for {free_name} in {work_dir}.\
                Check earlier error message if these values could be analysed.''')

    # get average of free energy values just calculated
    # if there is only one repeat use the BSS difference function
    if len(free_val_dict.values()) == 1 and len(bound_val_dict.values()) == 1:
        freenrg_rel = BSS.FreeEnergy.Relative.difference(
            list(bound_pmf_dict.items())[0][1], list(free_pmf_dict.items())[0][1])
        freenrg_val = freenrg_rel[0].value()
        freenrg_err = freenrg_rel[1].value()

    # otherwise, calculate the average and the SEM
    else:
        free_vals = list(
            val/_Units.Energy.kcal_per_mol for val in free_val_dict.values())
        free_avg = _np.mean(free_vals)
        free_sem = sem(free_vals)
        bound_vals = list(
            val/_Units.Energy.kcal_per_mol for val in bound_val_dict.values())
        bound_avg = _np.mean(bound_vals)
        bound_sem = sem(bound_vals)
        freenrg_val = (bound_avg-free_avg)
        freenrg_err = (_math.sqrt(
            _math.pow(bound_sem, 2)+_math.pow(free_sem, 2)))
        freenrg_rel = (freenrg_val * _Units.Energy.kcal_per_mol,
                        freenrg_err * _Units.Energy.kcal_per_mol)
    
    # create tuple list of each repeat that was calculated
    # first check the length of the calculated values and check if this is also the length of the folders found
    if len(bound_calculated) != no_of_b_repeats:
        print("the number of calculated values for bound does not match the number of found bound folders earlier.\
            Check previous messages to see which folder(s) couldn't be analysed.")
    if len(free_calculated) != no_of_f_repeats:
        print("the number of calculated values for free does not match the number of found free folders earlier.\
            Check previous messages to see which folder(s) couldn't be analysed.")   

    # if the numebr of calculated values is the same, match these evenly
    if len(bound_calculated) == len(free_calculated):
        print(f"There are {len(bound_calculated)} calculated values for each the bound and the free leg for the folders in {work_dir}.")
        no_of_repeats = len(bound_calculated)
        repeats = list(range(no_of_repeats))
        for r in repeats:
            freenrg_rel = BSS.FreeEnergy.Relative.difference(
                bound_pmf_dict[bound_calculated[r]], free_pmf_dict[free_calculated[r]])
            freenrg_val = freenrg_rel[0].value()
            freenrg_err = freenrg_rel[1].value()
            repeats_tuple_list.append((f"{str(r)}_repeat", freenrg_val, freenrg_err))

    elif len(bound_calculated) != len(free_calculated):
        print(f"There are {len(bound_calculated)} calculated values for the bound and \
        {len(free_calculated)} calculated values for the free leg for the folders in {work_dir}.")
        # use the shorter calculated values as the number of complete repeats
        if len(bound_calculated) < len(free_calculated):
            no_of_repeats = len(bound_calculated)
        else:
            no_of_repeats = len(free_calculated)
        print(f"The number of calculated values do not match. {no_of_repeats} repeats will be calculated.")
        r = 0
        for b,f in zip(bound_calculated, free_calculated):
            print(f"calculating repeat {r} as {b} and {f}.")
            freenrg_rel = BSS.FreeEnergy.Relative.difference(bound_pmf_dict[b], free_pmf_dict[f])
            freenrg_val = freenrg_rel[0].value()
            freenrg_err = freenrg_rel[1].value()
            repeats_tuple_list.append((f"{str(r)}_repeat", freenrg_val, freenrg_err))
            r += 1
    
    graph_dir = work_dir + '/graphs'
    if not _os.path.exists(graph_dir):
        _os.mkdir(f"{graph_dir}")
    else:
        pass

    if check_overlap:
        # check overlap matrix if okay
        for b in b_repeats:
            try:
                name = str(b) + '_bound'
                overlap = bound_matrix_dict[name]
                overlap_okay = BSS.FreeEnergy.Relative.checkOverlap(
                    overlap, estimator=estimator)
            except Exception as e:
                print(e)
                print(f"could not check overlap matrix for {name}")
            if save_graphs:
                try:
                    ax = BSS.FreeEnergy.Relative.plot(overlap, work_dir=graph_dir, plot_name=f"{name}_overlap_MBAR")
                except:
                    print("could not plt overlap matrix")

        for f in f_repeats:
            try:
                name = str(f) + '_free'
                overlap = free_matrix_dict[name]
                overlap_okay = BSS.FreeEnergy.Relative.checkOverlap(
                    overlap, estimator=estimator)
            except Exception as e:
                print(e)
                print(f"could not check overlap matrix for {name}")
            if save_graphs:
                try:
                    ax = BSS.FreeEnergy.Relative.plot(overlap, work_dir=graph_dir, plot_name=f"{name}_overlap_MBAR")
                except:
                    print("could not plt overlap matrix")


    if estimator == "TI" and save_graphs:
        
        for b in b_repeats:
            name = str(b) + '_bound'
            overlap = bound_matrix_dict[name]
            try:
                ax = BSS.FreeEnergy.Relative.plot(overlap, work_dir=graph_dir, plot_name=f"{name}_dHdl_TI")
            except Exception as e:
                print(e)
                print("could not plt dhdl")

        for f in f_repeats:
            name = str(f) + '_free'
            overlap = free_matrix_dict[name]
            try:
                ax = BSS.FreeEnergy.Relative.plot(overlap, work_dir=graph_dir, plot_name=f"{name}_dHdl_TI")
            except Exception as e:
                print(e)
                print("could not plt dhdl")

    if save_pickle:
        print("saving the pmf dictionaries for bound and free as pickles.")
        # write the pmf as a pickle
        with open (f"{pickle_dir}/bound_pmf_{trans}_{engine}_{estimator}_{method}.pickle", 'wb') as handle:
            pickle.dump(bound_pmf_dict, handle)
        with open (f"{pickle_dir}/free_pmf_{trans}_{engine}_{estimator}_{method}.pickle", 'wb') as handle:
            pickle.dump(free_pmf_dict, handle)
    else:
        pass

    return (freenrg_rel[0], freenrg_rel[1], repeats_tuple_list)





def calc_mae(values_dict=None, perts_ligs = None):
    # calc mae for a provided dictionary in the format wanted

    values_dict_list = []
    for key in values_dict.keys():
        values_dict_list.append(key)

    mae_pert_df = pd.DataFrame(columns=values_dict_list,index=values_dict_list)
    mae_pert_df_err = pd.DataFrame(columns=values_dict_list,index=values_dict_list)

    # iterate over all possible combinations
    for combo in itertools.product(values_dict_list,values_dict_list):
        eng1 = combo[0]
        eng2 = combo[1]

        eng1_vals = []
        eng2_vals = []

        # first create df of values
        if perts_ligs == "perts":
            for pert in values_dict[eng1]["pert_results"]:
                if pert in values_dict[eng2]["pert_results"]:
                    if values_dict[eng1]["pert_results"][pert][0] != None:
                        if values_dict[eng2]["pert_results"][pert][0] != None:
                            eng1_vals.append(values_dict[eng1]["pert_results"][pert][0])
                            eng2_vals.append(values_dict[eng2]["pert_results"][pert][0])
        elif perts_ligs == "ligs":
            for pert in values_dict[eng1]["val_results"]:
                if pert in values_dict[eng2]["val_results"]:
                    if values_dict[eng1]["val_results"][pert][0] != None:
                        if values_dict[eng2]["val_results"][pert][0] != None:
                            eng1_vals.append(values_dict[eng1]["val_results"][pert][0])
                            eng2_vals.append(values_dict[eng2]["val_results"][pert][0])
        else:
            raise ValueError("'perts_ligs' must be either 'perts' or 'ligs'!")


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

    mae_pert_df.to_csv(f"{res_folder}/mae_pert_{file_ext_out}.csv", sep=" ")
    mae_pert_df_err.to_csv(f"{res_folder}/mae_pert_err_{file_ext_out}.csv", sep=" ")

    return mae_pert_df, mae_pert_df_err

