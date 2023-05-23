
# import libraries
import BioSimSpace as BSS
import csv
import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import sem as sem
from scipy.stats import bootstrap
from sklearn.metrics import mean_absolute_error as mae
from rdkit import Chem

import matplotlib.pyplot as plt
import csv
import numpy as np
import pandas as pd 

from ..utils import *

def get_ligands_from_perts(perturbations):

    perturbations = validate.is_list(perturbations)
    ligands = []

    for pert in perturbations:
        lig_0 = pert.split("~")[0]
        lig_1 = pert.split("~")[1]        
        if lig_0 not in ligands:
            ligands.append(lig_0)
        if lig_1 not in ligands:
            ligands.append(lig_1)
    
    return ligands

def get_info_network(net_file=None, results_files=None, extra_options=None):
    """get information about the network from the network file

    Args:
        net_file (str, optional): network file. Defaults to None.
        results_files (list, optional): list of results files. Defaults to None.
        extra_options (dict, optional): extra options (engine or engines). Defaults to None.

    Returns:
        tuple: (perturbations, ligands)
    """
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
            results_files = validate.is_list(results_files, make_list=True)
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
            try:
                engines = [validate.engine(extra_options["engine"])]
            except Exception as e:
                print(e)
        if "engines" in extra_options.keys():
            try:
                engines = validate.engines(extra_options["engines"])
            except Exception as e:
                print(e)
    
    # We also want to create a list of the perturbations in our network.
    perturbations = []

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
    
    else:
        for res_file in results_files:
            # use the network file to find the ligands and perturbations
            with open(f"{res_file}", "r") as file:
                lines = (line.rstrip() for line in file) # All lines including the blank ones
                lines = (line for line in lines if line) # Non-blank lines
                for line in lines:
                    for engine in engines:
                        if line.split(",")[4] == engine:
                            lig_0 = line.split(",")[0]
                            lig_1 = line.split(",")[1]
                            pert = f"{lig_0}~{lig_1}"
                            if pert not in perturbations:
                                perturbations.append(pert)

    ligands = get_ligands_from_perts(perturbations)    

    return (perturbations, ligands)


def get_info_network_from_dict(res_dict):
    """get list of perturbations and ligands from a dictionary

    Args:
        res_dict (dict): dictionary of free energy results

    Returns:
        tuple: (perturbations, ligands)
    """
    # get info for the network from a perturbation results dictionary

    res_dict = validate.dictionary(res_dict)

    perturbations = []

    for key in res_dict.keys():
        if key not in perturbations:
            perturbations.append(key)
    
    ligands = get_ligands_from_perts(perturbations)

    return (perturbations, ligands)


class net_graph():
    
    def __init__(self, ligands, perturbations, file_dir=None, ligands_folder=None):
        """_summary_

        Args:
            ligands (list): list of ligands
            perturbations (list): list of perturbations
            file_dir (str, optional): folder path to save graph image in. Defaults to None.
        """

        self.ligands = validate.is_list(ligands)
        self.perturbations = validate.is_list(perturbations)

        if file_dir:
            self.file_dir = validate.folder_path(file_dir, create=True)
            self._save_image = True
        else:
            self._save_image = False

        if ligands_folder:
            self.ligands_folder = validate.folder_path(ligands_folder)
        else:
            self.ligands_folder = None

        net_graph._gen_graph(self)

    def _gen_graph(self):
        """generate a network x graph from the perturbations and ligands.
        """

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
        """draw the network x graph.

        Args:
            file_dir (str, optional): where to save the image. Defaults to None, no image is saved.
        """
        
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
            plt.savefig(f"{self.file_dir}/network.png", dpi=300)
        if file_dir:
            file_dir = validate.folder_path(file_dir, create = True)
            plt.savefig(f"{file_dir}/network.png", dpi=300)

        plt.show()

    def draw_ligand(self, ligand):

        if not self.ligands_folder:
            raise ValueError("please provide a ligands dir w the files inside.")
        
        m = Chem.SDMolSupplier(f"{self.ligands_folder}/{ligand}.sdf")[0]
        smi = Chem.MolToSmiles(m)
        m2 = Chem.MolFromSmiles(smi)
        img = Chem.Draw.MolToImage(m2)

        return img

    def draw_perturbation(self, pert):

        lig_0 = pert.split("~")[0]
        lig_1 = pert.split("~")[1]   

        img = self.draw_ligand(lig_0)
        img2 = self.draw_ligand(lig_1)

        # plt.rcParams['figure.figsize'] = 11 ,8
        fig, ax = plt.subplots(1,2)

        ax[0].imshow(img)
        ax[0].title.set_text(f"{lig_0}")
        ax[0].axis("off")
        ax[1].imshow(img2)
        ax[1].title.set_text(f"{lig_1}")
        ax[1].axis("off")

        return fig

    def disconnected_ligands(self):

        ligs = [lig for lig in nx.isolates(self.graph)]

        return ligs

    def cycle_closures(self):
        """get cycle closures in the network

        Returns:
            list: list of cycle closures
        """

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


    def add_weight(self, input_data):
        """add weights to the network x graph for the edges and each perturbation.

        Args:
            input_data (str or dict, optional): file path or dict of the weights for the perturbation edges.

        Raises:
            TypeError: if dict, must be of format {(lig_0, lig_1): weight}
        """
        # these weights are for the lomap or rbfenn score
        
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