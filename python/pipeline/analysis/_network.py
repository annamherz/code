
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



# from alvaro
def get_average_weighted_simple_paths(G):
    '''Calculate the average number of connection between each pair of nodes. 
    '''
    paths_per_nodepair_combination = []
    for node_i in G.nodes:
        for node_j in G.nodes:
            if node_i == node_j: break
            possible_paths = nx.all_simple_edge_paths(G, node_i, node_j)
            sum_of_weighted_averaged_paths = sum([np.average([G.get_edge_data(*edge)['weight']
                                                                                      for edge in path])
                                                                                      for path in possible_paths])
            paths_per_nodepair_combination.append(sum_of_weighted_averaged_paths)
    return np.average(paths_per_nodepair_combination)