
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

# functions
# TODO clean up and make sure have description at start

def convert_yml_into_freenrgworkflows(exp_file, exp_file_dat):
    # get the experimental data into a useable format (from yml to csv)
    # for freenergworkflows, want to save as lig, Ki
    # experimental values (e.g. ic50/ki) for all ligands in our set.
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



def get_info_network(engine=None, results_files=None, net_file=None, output_folder=None, extra_options=None):
    # get info from a network file for engine
    # For the network that we are considering,
    # we want to get results files with these perturbations from the overall file that contains the large network results (if this is the case).
    # This is just good to have a consistent format of results to analyse.
    if not engine:
        raise ValueError("must specify an engine") # TODO make so can be BSS engines or all engines allowed
    # TODO if not out_folder is a file path, raise issue
    # if not engine:
    #     raise ValueError("must specify an engine")    
    # TODO if not temp folder make

    # get all the indivisual results file for that engine
    # TODO if results file is empty, raise issue
    # TODO print how many results files
    # TODO extra options if want to change file name etc, also put if only want a specific engine here?
    # otherwise will find all engines in network file and use this

    # We also want to create a list of the perturbations in our network.
    # create a list of the perturbations
    perturbations = []

    # create a list of ligands
    ligands = []

    # use the network file to find the ligands and perturbations
    for line in open(f"{net_file}", "r"):
        if line.split()[-1] == engine:
            lig_0 = line.split()[0]
            lig_1 = line.split()[1]
            pert = f"{lig_0}~{lig_1}"
            perturbations.append(pert)
            if lig_0 not in ligands:
                ligands.append(lig_0)
            elif lig_1 not in ligands:
                ligands.append(lig_1)
            else:
                pass

    mod_results_files = []

    for file in results_files:
        new_file_name = f"{output_folder}/results_{results_files.index(file)}_{engine}.csv"
        with open(new_file_name, "w") as result_file:

            writer = csv.writer(result_file, delimiter=",")
            writer.writerow(["lig_1","lig_2","freenrg","error","engine"])

            for row, index in pd.read_csv(file).iterrows():
                pert = f"{index['lig_1']}~{index['lig_2']}"
                if pert in perturbations:
                    writer.writerow([index['lig_1'], index['lig_2'], index['freenrg'], index['error'], index['engine']])    

            mod_results_files.append(new_file_name)

    return perturbations, ligands, mod_results_files


def gen_graph(ligands=None, perturbations=None, file_dir=None):

    # TODO error check if not lists, etc
    # Generate the graph.
    graph = nx.Graph()

    # Loop over the nligands and add as nodes to the graph.
    for lig in ligands:
        graph.add_node(lig, label=lig, labelloc="t")

    # Loop over the edges in the dictionary and add to the graph.
    for edge in perturbations:
        lig_0 = edge.split("~")[0]
        lig_1 = edge.split("~")[1]
        graph.add_edge(lig_0, lig_1)

    # Plot the networkX graph.
    pos = nx.kamada_kawai_layout(graph)
    plt.figure(figsize=(8,8), dpi=150)
    nx.draw(
        graph, pos, edge_color='black', width=1, linewidths=1,
        node_size=2100, node_color='darkblue', font_size = 9.5,
        labels={node: node for node in graph.nodes()},
        font_color = "white")

    plt.savefig(f"{file_dir}/analysis_network.png", dpi=300)
    plt.show()


# make dictionary of values of results
def make_dict_comp_results(results_files=None, perturbations=None, file_name=None, engine=None):

    # TODO some way to get engine from results files?
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
            if not isinstance(row[2], float):
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
                comp_dict_list[pert] = ddG
                comp_err_dict_list[pert] = ddG_err

    # now calculate all the avg and SEM for the network perturbations
    # put these into a dictionary
    comp_diff_dict = {}

    # TODO if none, don't write to csv file
    # write these to a csv file
    with open(f"{file_name}.csv", "w") as comp_pert_file:
        writer = csv.writer(comp_pert_file, delimiter=",")
        writer.writerow(["lig_1","lig_2","freenrg","error","engine"])
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
            
            # TODO some way to incl units

            #update the dictionary for plotting later
            comp_diff_dict.update({pert:(comp_ddG, comp_err)})

            writer.writerow([lig_0, lig_1, comp_ddG, comp_err, engine])
    
    return comp_diff_dict

def freenrgworkflows_into_dict(experimental_DDGs, ligands, perturbations):
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

    # now that we have our dictionary, 
    # we can also create a dictionary with all the experimental values for the perturbations
    exper_diff_dict = {}

    # calculate the experimental RBFEs
    # write these to a csv file
    with open("experimental_perturbations.csv", "w") as exp_pert_file:
        writer = csv.writer(exp_pert_file, delimiter=",")
        writer.writerow(["lig_1","lig_2","freenrg","error","engine"])

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

            writer.writerow([lig_0, lig_1, exper_ddG, exper_err, "experimental"])

    return exper_diff_dict, exper_val_dict


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

def gen_graph(ligands=None, perturbations=None, file_dir=None):

    # TODO error check if not lists, etc
    # Generate the graph.
    graph = nx.Graph()

    # Loop over the nligands and add as nodes to the graph.
    for lig in ligands:
        graph.add_node(lig, label=lig, labelloc="t")

    # Loop over the edges in the dictionary and add to the graph.
    for edge in perturbations:
        lig_0 = edge.split("~")[0]
        lig_1 = edge.split("~")[1]
        graph.add_edge(lig_0, lig_1)

    # Plot the networkX graph.
    pos = nx.kamada_kawai_layout(graph)
    plt.figure(figsize=(8,8), dpi=150)
    nx.draw(
        graph, pos, edge_color='black', width=1, linewidths=1,
        node_size=2100, node_color='darkblue', font_size = 9.5,
        labels={node: node for node in graph.nodes()},
        font_color = "white")

    plt.savefig(f"{file_dir}/analysis_network.png", dpi=300)
    plt.show()

    return graph

