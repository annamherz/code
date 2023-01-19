
# import libraries
import BioSimSpace as BSS
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
import itertools

import matplotlib.pyplot as plt
import csv
import numpy as np
import pandas as pd 

from ..utils import *
from ._convert import *

# functions
# TODO clean up and make sure have description at start

class make_dict():
    """class of static methods for making dicts to analyse the results
    """

    def __init__():
        pass

    
    @staticmethod
    def comp_results(results_files=None, perturbations=None, engine=None, output_file=None):
        
        results_files = validate.is_list(results_files)
        for file in results_files:
            validate.file_path(file)
        
        if perturbations:
            perturbations = validate.is_list(perturbations)
            make_pert_list = False
        else:
            make_pert_list = True
            perturbations = []

        # TODO some way to get engine from results files? - or only read for certain engines
        # make a dictionary with the results of the files
        comp_dict_list = {}
        comp_err_dict_list = {}

        # append for results file
        for res_file in results_files:
            res_df = pd.read_csv(res_file)
            for index,row in res_df.iterrows():
                lig_0 = row[0]
                lig_1 = row[1]
                pert = f"{lig_0}~{lig_1}"

                if make_pert_list:
                    if pert not in perturbations:
                        perturbations.append(pert)

                if not isinstance(row[2], float):
                    # to convert ?
                    ddG = BSS.Types.Energy(float(row[2].split()[0]),row[2].split()[-1])
                else:
                    ddG = row[2]
                if not isinstance(row[3], float):
                    ddG_err = BSS.Types.Energy(float(row[3].split()[0]),row[3].split()[-1])
                else:
                    ddG_err = row[3]
                    
                if pert in comp_dict_list:
                    # Key exist in dict, check if is a list
                    if not isinstance(comp_dict_list[pert], list):
                        # If type is not list then make it list
                        comp_dict_list[pert] = [comp_dict_list[pert]]
                    if not isinstance(comp_err_dict_list[pert], list):
                        # If type is not list then make it list
                        comp_err_dict_list[pert] = [comp_err_dict_list[pert]]
                    # Append the value in list
                    comp_dict_list[pert].append(ddG)
                    comp_err_dict_list[pert].append(ddG_err)
                else:
                    # As key is not in dict,
                    # so, add key-value pair
                    comp_dict_list[pert] = [ddG]
                    comp_err_dict_list[pert] = [ddG_err]

        # now calculate all the avg and SEM for the network perturbations
        # put these into a dictionary
        comp_diff_dict = {}

        # TODO combine this with file writer for network in _files
        if output_file:
            # write these to a csv file
            with open(f"{output_file}.csv", "w") as comp_pert_file:
                writer = csv.writer(comp_pert_file, delimiter=",")
                if engine:
                    writer.writerow(["lig_1","lig_2","freenrg","error","engine"])
                else:
                    writer.writerow(["lig_1","lig_2","freenrg","error","source"])
                
                for pert in perturbations:
                    lig_0 = pert.split("~")[0]
                    lig_1 = pert.split("~")[1]
                    
                    # check if the perturbations calculated are also those in the network file and if any are missing
                    try:
                        # find the values in the dictionary
                        ddGs = comp_dict_list[pert]
                        ddGs_error = comp_err_dict_list[pert]
                        # calculate the average and the error
                        comp_ddG = np.average(ddGs)
                        # comp_ddG = np.average([ddG.value() for ddG in ddGs])
                        if len(ddGs) == 1:
                            # comp_err = ddGs_error.value()
                            comp_err = ddGs_error[0]
                        else:
                            comp_err = sem(ddGs)
                            # comp_err = sem([ddG.value() for ddG in ddGs])
                
                    # if unable to calculate one of the perturbations, this is a None value.
                    except:
                        comp_ddG = None
                        comp_err = None

                    #update the dictionary for plotting later
                    comp_diff_dict.update({pert:(comp_ddG, comp_err)})

                    if engine:
                        writer.writerow([lig_0, lig_1, comp_ddG, comp_err, engine])
                    else:
                        writer.writerow([lig_0, lig_1, comp_ddG, comp_err, "source"])


        for pert in perturbations:
            lig_0 = pert.split("~")[0]
            lig_1 = pert.split("~")[1]
            
            # check if the perturbations calculated are also those in the network file and if any are missing
            try:
                # find the values in the dictionary
                ddGs = comp_dict_list[pert]
                ddGs_error = comp_err_dict_list[pert]
                # calculate the average and the error
                comp_ddG = np.average(ddGs)
                # comp_ddG = np.average([ddG.value() for ddG in ddGs])
                if len(ddGs) == 1:
                    comp_err = ddGs_error.value()
                else:
                    comp_err = sem(ddGs)
                    # comp_err = sem([ddG.value() for ddG in ddGs])
        
            # if unable to calculate one of the perturbations, this is a None value.
            except:
                comp_ddG = None
                comp_err = None
            
            # TODO some way to incl units - attributes?

            #update the dictionary for plotting later
            comp_diff_dict.update({pert:(comp_ddG, comp_err)})
        

        return comp_diff_dict


    @staticmethod
    def experimental_from_freenrgworkflows(experimental_DDGs, ligands, perturbations):
        
        ligands = validate.is_list(ligands)
        perturbations = validate.is_list(perturbations)

        exper_val_dict = make_dict._from_freenrgworkflows_experimental_val(experimental_DDGs, ligands)
        exper_diff_dict = make_dict._from_freenrgworkflows_experimental_diff(exper_val_dict, perturbations)

        return exper_diff_dict, exper_val_dict
    
    @staticmethod
    def _from_freenrgworkflows_experimental_val(experimental_DDGs, ligands):

        # create a dictionary for the experimental values
        exper_val_dict = {}

        # convert the list of dicitonaries from freenrgworkflows into a single dictionary
        for lig_dict in experimental_DDGs:
            lig_name = list(lig_dict.keys())[0]
            exper = lig_dict[lig_name]
            exper_err = lig_dict["error"]
            exper_val_dict.update({lig_name:(exper, exper_err)})

        # add any ligands that are in the ligands file but dont have experimental values for
        for lig_name in ligands:
            if lig_name in exper_val_dict:
                pass
            else:
                exper_val_dict.update({lig_name:(None, None)})
        
        return exper_val_dict

    @staticmethod
    def _from_freenrgworkflows_experimental_diff(exper_val_dict, perturbations):

        # we can also create a dictionary with all the experimental values for the perturbations
        exper_diff_dict = {}

        # calculate the experimental RBFEs

        for pert in perturbations:
            lig_0 = pert.split("~")[0]
            lig_1 = pert.split("~")[1]
            # exclude from calculating if one of the ligands is not available
            if exper_val_dict[lig_0][0] is None or exper_val_dict[lig_1][0] is None:
                exper_ddG =None
                exper_err = None
                exper_diff_dict.update({pert:(None, None)})
            # if experimental data is available, calculate experimental perturbation and propagate
            else:
                exper_ddG = exper_val_dict[lig_1][0] - exper_val_dict[lig_0][0]
                exper_err = math.sqrt(math.pow(exper_val_dict[lig_0][1], 2) + math.pow(exper_val_dict[lig_1][1], 2))
                exper_diff_dict.update({pert:(exper_ddG, exper_err)})
        
        return exper_diff_dict

    @staticmethod
    def from_freenrgworkflows_network_analyser(computed_relative_DDGs):

        freenrg_dict = {}

        # append computed freenrg and error.
        for item in computed_relative_DDGs:
            ligand = list(item.keys())[0]
            freenrg = list(item.values())[0]
            error = list(item.values())[1]

            freenrg_dict.update({ligand:(freenrg, error)})

        return freenrg_dict

    @staticmethod
    def from_cinnabar_network_edges(network, calc_exp, perturbations):

        if calc_exp not in ["calc","exp"]:
            raise ValueError("calc_exp must be either 'calc' or 'exp'")

        freenrg_dict = {}

        name_dict = {}
        for node in network.graph.nodes(data=True):
            name_dict.update({node[0]: node[1]["name"]})

        for edge in network.graph.edges(data=True):
            lig_0 = name_dict[edge[0]]
            lig_1 = name_dict[edge[1]]
            pert = f"{lig_0}~{lig_1}"
            anti_pert = f"{lig_1}~{lig_0}"
            if pert in perturbations:
                pert_name = pert
                add_dict = True
            elif anti_pert in perturbations:
                pert_name = anti_pert
                add_dict = True
            else:
                add_dict = False
            
            if add_dict:
                freenrg_dict.update({pert_name:(edge[2][f"{calc_exp}_DDG"], edge[2][f"{calc_exp}_dDDG"])})
        
        return freenrg_dict

    @staticmethod
    def from_cinnabar_network_node(network, calc_exp, normalise=False):

        if calc_exp not in ["calc","exp"]:
            raise ValueError("calc_exp must be either 'calc' or 'exp'")

        freenrg_dict = {}

        for node in network.graph.nodes(data=True):
            freenrg_dict.update({node[1]["name"]:(node[1][f"{calc_exp}_DG"], node[1][f"{calc_exp}_dDG"])})
        
        if normalise:
            normalised_freenrg_dict = make_dict._normalise_data(freenrg_dict)
            return normalised_freenrg_dict
        else:
            return freenrg_dict

    @staticmethod
    def experimental_for_network(exper_dict, ligands, perturbations):
        
        ligands = validate.is_list(ligands)
        perturbations = validate.is_list(perturbations)

        exper_val_dict = make_dict._exper_from_ligands(exper_dict, ligands)
        exper_diff_dict = make_dict._exper_from_perturbations(exper_val_dict, perturbations)

        return exper_diff_dict, exper_val_dict

    @staticmethod
    def _exper_from_ligands(exper_val_dict, ligands, normalise=False):

        exper_val_dict = validate.dictionary(exper_val_dict)
        ligands = validate.is_list(ligands)
        normalise = validate.boolean(normalise)

        new_exper_val_dict ={}

        for lig in ligands:
            exper_dG = exper_val_dict[lig][0]
            exper_err = exper_val_dict[lig][1]
            new_exper_val_dict.update({lig:(exper_dG, exper_err)})
        
        if normalise:
            normalised_exper_val_dict = make_dict._normalise_data(new_exper_val_dict)
            return normalised_exper_val_dict
        else:
            return new_exper_val_dict


    @staticmethod
    def _exper_from_perturbations(exper_val_dict, perturbations):

        exper_diff_dict = {}

        # calculate the experimental RBFEs
        for pert in perturbations:
            lig_0 = pert.split("~")[0]
            lig_1 = pert.split("~")[1]
            exper_ddG = exper_val_dict[lig_1][0] - exper_val_dict[lig_0][0]
            exper_err = math.sqrt(math.pow(exper_val_dict[lig_0][1], 2) + math.pow(exper_val_dict[lig_1][1], 2))
            exper_diff_dict.update({pert:(exper_ddG, exper_err)})

        return exper_diff_dict

    @staticmethod
    def _normalise_data(data):
        # normalise either a list / np / dict of data

        try:
            data_dict = validate.dictionary(data)
            is_dict = True
        except:
            data = validate.is_list(data)
            is_dict = False
        
        if is_dict:
            normalised_dict = {}
            values = []
            for val in data_dict.values():
                values.append(float(val[0])) # incase it is a tuple of vals
            avg_val = np.mean(values)
            for key in data_dict:
                normalised_dict.update({key:(float(data_dict[key][0]) - avg_val, data_dict[key][1])}) # also incl error

            return normalised_dict
        
        else:
            avg_val = np.mean(data)
            normalised_data = []
            for val in data:
                normalised_data.append(val - avg_val)
            
            return normalised_data


    @staticmethod
    def cycle_closures(pert_dict, cycle_closures):

        pert_dict = validate.dictionary(pert_dict)
        cycle_closures = validate.is_list(cycle_closures)

        cycles_dict = {}
        cycle_vals = []

        for cycle in cycle_closures:

            cycle_dict = {}
            # print(cycle)
            cycle_val = []
            cycle_val_err = []
            for pert in cycle:

                liga = pert.split("~")[0]
                ligb = pert.split("~")[1]
                rev_pert = f"{ligb}~{liga}"
            
                if pert in pert_dict:
                    if pert_dict[pert][0] is not None:
                        cycle_val.append(+pert_dict[pert][0])
                        cycle_val_err.append(pert_dict[pert][1])
                    else:
                        print(f"{pert} or {rev_pert} does not exist in the results for {cycle}. This cycle is not included.")
                        break
                elif rev_pert in pert_dict:
                    if pert_dict[rev_pert][0] is not None:
                        cycle_val.append(-pert_dict[rev_pert][0])
                        cycle_val_err.append(pert_dict[rev_pert][1])
                    else:
                        print(f"{pert} or {rev_pert} does not exist in the results for {cycle}. This cycle is not included.")
                        break

            if not all(i is None for i in cycle_val):
                cycle_vals.append(sum(cycle_val))
            else:
                pass

            cycles_dict.update({"_".join(cycle):(sum(cycle_val), cycle_val_err)}) # TODO also incl the error for a cycle depending on whether or not std or other error

        return (cycles_dict, cycle_vals, np.mean(cycle_vals), np.std(cycle_vals))

