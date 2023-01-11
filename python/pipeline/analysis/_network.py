
# import libraries
import BioSimSpace as BSS
import csv
import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import sem as sem
from scipy.stats import bootstrap
from sklearn.metrics import mean_absolute_error as mae

import matplotlib.pyplot as plt
import csv
import numpy as np
import pandas as pd 

from ..utils import *

def get_info_network(net_file=None, results_files=None, extra_options=None):
    # get info from a network file for engine
    # For the network that we are considering,
    # we want to get results files with these perturbations from the overall file that contains the large network results (if this is the case).
    # This is just good to have a consistent format of results to analyse.
    
    use_net_file = False

    if net_file:
        try:
            net_file = validate.file_path(net_file)
            use_net_file = True
        except:
            print("can't use net_file, will use results files instead.")

    if not use_net_file:
        try:
            results_files = validate.is_list(results_files)
            for file in results_files:
                validate.file_path(file)
        except Exception as e:
            print(e)
            print("cant use network or results files, please provide one or the other.")
            raise
    
    # set extra_options variables as defaults
    engines = [eng.upper() for eng in BSS.FreeEnergy.engines()] # use all

    if extra_options:
        extra_options = validate.dictionary(extra_options)

        if "engine" in extra_options.keys():
            engines = [validate.engine(extra_options["engine"])]
        if "engines" in extra_options.keys():
            engines = validate.is_list(extra_options["engines"])
            for engine in engines:
                engine_val = validate.engine(engine)
                engines = [engine_val if i == engine else i for i in engines]
       
    
    # We also want to create a list of the perturbations in our network.
    perturbations = []
    # create a list of ligands
    ligands = []

    if use_net_file:
        # use the network file to find the ligands and perturbations
        with open(f"{net_file}", "r") as file:
            for line in file:
                for engine in engines:
                    if line.split()[-1] == engine:
                        lig_0 = line.split()[0]
                        lig_1 = line.split()[1]
                        pert = f"{lig_0}~{lig_1}"
                        if pert not in perturbations:
                            perturbations.append(pert)
                        if lig_0 not in ligands:
                            ligands.append(lig_0)
                        if lig_1 not in ligands:
                            ligands.append(lig_1)
                        else:
                            pass
    
    else:
        for res_file in results_files:
            # use the network file to find the ligands and perturbations
            with open(f"{res_file}", "r") as file:
                for line in file:
                    for engine in engines:
                        if line.split(",")[4] == engine:
                            lig_0 = line.split(",")[0]
                            lig_1 = line.split(",")[1]
                            pert = f"{lig_0}~{lig_1}"
                            if pert not in perturbations:
                                perturbations.append(pert)
                            if lig_0 not in ligands:
                                ligands.append(lig_0)
                            if lig_1 not in ligands:
                                ligands.append(lig_1)
                            else:
                                pass     

    return (perturbations, ligands)
            
class net_graph():
    
    def __init__(self, ligands, perturbations, file_dir=None):

        self.ligands = validate.is_list(ligands)
        self.perturbations = validate.is_list(perturbations)

        if file_dir:
            self.file_dir = validate.folder_path(file_dir, create=True)
            self._save_image = True
        else:
            self._save_image = False

        net_graph._gen_graph(self)

    def _gen_graph(self):

        # Generate the graph.
        graph = nx.Graph()

        # Loop over the nligands and add as nodes to the graph.
        for lig in self.ligands:
            graph.add_node(lig, label=lig, labelloc="t")

        # Loop over the edges in the nxgraph and add to the graph.
        for edge in self.perturbations:
            lig_0 = edge.split("~")[0]
            lig_1 = edge.split("~")[1]
            graph.add_edge(lig_0, lig_1)
        
        self.graph = graph


    def draw_graph(self, file_dir=None):
        
        graph = validate.nxgraph(self.graph)

        # Plot the networkX graph.
        pos = nx.kamada_kawai_layout(graph)
        plt.figure(figsize=(8,8), dpi=150)
        nx.draw(
            graph, pos, edge_color='black', width=1, linewidths=1,
            node_size=2100, node_color='darkblue', font_size = 9.5,
            labels={node: node for node in graph.nodes()},
            font_color = "white")

        if self._save_image:
            plt.savefig(f"{self.file_dir}/analysis_network.png", dpi=300)
        if file_dir:
            file_dir = validate.folder_path(file_dir, create = True)
            plt.savefig(f"{file_dir}/analysis_network.png", dpi=300)

        plt.show()


    def cycle_closures(self):

        cycles = nx.cycle_basis(self.graph)

        # get list of all cycle closures
        cycle_closures = []

        for cycle in cycles:
            length = len(cycle)
            ligas = []
            ligbs = []
            for i in range(0, length-1):
                ligas.append(cycle[i])
                ligbs.append(cycle[i+1])
            # add final cycle closure
            ligas.append(cycle[-1])
            ligbs.append(cycle[0])
            # make list for cycle closure
            cycle_closure = []
            for liga, ligb in zip(ligas,ligbs):
                cycle_closure.append(f"{liga}~{ligb}")
            
            cycle_closures.append(cycle_closure)

            self.cycles_list = cycle_closures    
        
        return cycle_closures


    def add_weight(self, input_data=None):
        
        if not input_data:
            raise ValueError("need some kind of input to add weights from! dict or file.")

        else:
            try:
                weight_dict = validate.dictionary(input_data)
                for key in weight_dict.keys():
                    if not isinstance(key, tuple):
                        raise TypeError("dict entry must be of the format {(lig_0, lig_1): weight}")
                use_file = False
                # print("using dict to add weights")
            except:
                use_file = True
                weight_dict = {}

        if use_file:
            # print("using input file to get the weights")
            scores_file = validate.file_path(input_data)
            
            with open(scores_file, "r") as file:
                for line in file:
                    lig_0 = line.split(",")[0]
                    lig_1 = line.split(",")[1]
                    weight = line.split(",")[-1].strip()
                    weight_dict[(lig_0, lig_1)] = float(weight)            

        nx.set_edge_attributes(self.graph, weight_dict, name="weight")


    # from alvaro
    def get_average_weighted_simple_paths(self):
        '''Calculate the average number of connection between each pair of nodes. 
        '''

        G = self.graph

        paths_per_nodepair_combination = []
        for node_i in G.nodes:
            for node_j in G.nodes:
                if node_i == node_j: break
                possible_paths = nx.all_simple_edge_paths(G, node_i, node_j)
                sum_of_weighted_averaged_paths = sum([np.average([G.get_edge_data(*edge)['weight']
                                                                                        for edge in path])
                                                                                        for path in possible_paths])
                paths_per_nodepair_combination.append(sum_of_weighted_averaged_paths)
        
        # no of possible paths between each two nodes on average
        return np.average(paths_per_nodepair_combination)


# TODO in new - use functions from network and dictionary to make the results object
# this would also incl the plotting which uses the network as well