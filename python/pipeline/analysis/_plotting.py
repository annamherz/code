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
import pickle
import tempfile
import itertools as it
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import csv
import numpy as np
import pandas as pd 
from rdkit import Chem

from scipy.stats import sem as sem
from scipy.stats import bootstrap, norm
from sklearn.metrics import mean_absolute_error as mae

from ..utils import *
from ._network import *

class plotting_engines():

    def __init__(self, analysis_object=None, output_folder=None, verbose=False):
        """for plotting analysis network results.

        Args:
            analysis_object (pipeline.analysis.analysis_network, optional): analysis object that is to be plotted for. Defaults to None.
            output_folder (str, optional): output folder for generated csv files and graph images. Defaults to None.

        Raises:
            ValueError: must ptovide an analysis network object.
        """

        self.is_verbose(verbose)
        
        if analysis_object:
            self._analysis_object = analysis_object
            # get info about things for plotting from analysis
            self.analysis_obj_into_format()
        else:
            raise ValueError("please provide an analysis object to be plotted for.")

        # place to write results to
        if not output_folder:
            # want to write to the graph directory
            self.output_folder = self._analysis_object.output_folder
            self.graph_folder = self._analysis_object.graph_dir
        else:
            self.output_folder = validate.folder_path(output_folder, create=True)
            self.graph_folder = validate.folder_path(f"{output_folder}/graphs", create=True)

        # set the colours and bar spacing
        self._set_style()

        # convert the dictionaries into dataframes for plotting
        self._analysis_dicts_to_df()

    def is_verbose(self, value):

        verbose = validate.boolean(value)
        self._is_verbose = verbose

        return verbose

    def analysis_obj_into_format(self):
        """turn the passed pipeline.analysis.analysis_network object into format for this class
        """

        ana_obj = self._analysis_object
        
        # analysis information
        self.engines = sorted(ana_obj.engines)
        self.ligands = ana_obj.ligands
        self.perturbations = ana_obj.perturbations

        # for other results
        self.other_results_names = ana_obj.other_results_names
    
        # name of all options
        self._eng_other_list()
        
        # file extension
        self.file_extension(self.default_file_ext())
        self.network_extension(self.default_net_ext())

        # dictionaries of engines for plotting from cinnabar
        if ana_obj.cinnabar_calc_pert_dict:
            self.calc_val_dict = ana_obj.cinnabar_calc_val_dict
            self.exper_val_dict = ana_obj.cinnabar_exper_val_dict
            self.calc_pert_dict = ana_obj.cinnabar_calc_pert_dict
        else:
            print("no cinnabar calculation has been performed. Can only plot 'pert' values.")
            self.calc_val_dict = {}
            self.exper_val_dict = {}
            self.calc_pert_dict = ana_obj.calc_pert_dict

        # experimental calculated directly from exp values (for bar)
        self.all_exper_pert_dict = ana_obj.exper_pert_dict
        self.all_exper_val_dict = ana_obj.normalised_exper_val_dict

    def default_file_ext(self):
        """file extension from the analysis object, if none is provided na.

        Returns:
            str: file extension
        """

        file_ext = self._analysis_object.file_ext

        if file_ext == ".+":
            file_ext = "na"
 
        return file_ext
    
    def default_net_ext(self):
        """the network extension

        Returns:
            str: network extension
        """

        net_ext = self._analysis_object.net_ext
        self.net_ext = net_ext    

        return net_ext    
    
    def file_extension(self, file_ext=None):
        """set or return the file extension

        Args:
            file_ext (str, optional): file extension to set. Defaults to None.

        Returns:
            str: the file extension
        """

        if file_ext:
            self.file_ext = validate.string(file_ext)
        else:
            pass

        return self.file_ext

    def network_extension(self, net_ext=None):
        """set or return the network extension

        Args:
            file_ext (str, optional): network extension to set. Defaults to None.

        Returns:
            str: the network extension
        """

        if net_ext:
            self.net_ext = validate.string(net_ext)
        else:
            pass

        return self.net_ext

    def _eng_other_list(self):
        """list of engines and any other results and the experimental.

        Returns:
            list: names list
        """

        names_list = []

        # all possible engines
        for eng in self.engines:
            names_list.append(eng)
        # all other results
        for name in self.other_results_names:
            names_list.append(name)
        # the experimental
        names_list.append("experimental")
        
        self.names_list = names_list
            
        return names_list

    def _validate_in_names_list(self, name):
        """validate if the name is in the names list

        Args:
            name (str): the name to validate

        Raises:
            ValueError: if not in names list

        Returns:
            str: the validated name
        """

        name = validate.string(name)
        if name not in self.names_list:
            raise ValueError(f"{name} must be in {self.names_list}")

        return name            

    def _analysis_dicts_to_df(self):
        """turn the dicts from the analysis network to ones for the plotting object
        """

        # create the overall dict from all the passed results
        self._overall_dict()
        self.freenrg_df_dict = {}
        for name in self.names_list:
            self.freenrg_df_dict.update({name:None})
        self._calc_all_dict_to_df()


    def _overall_dict(self):
        """create an overall dict with the information passed from the analysis network object.

        Returns:
            dict: values dict for all the engines, other results, and experimental. Each dict in the dict contains the 'perts', 'ligs', 'pert_results', and 'val_results'.
        """

        values_dict = {}
        for name in self.names_list:
            values_dict.update({name:{}})

        # run for all engines with selected network and populate the dictionary for plotting
        for eng in (self.engines + self.other_results_names):
            try:
                # get perts and ligands for each engine
                pert_lig = get_info_network_from_dict(self.calc_pert_dict[eng])
                values_dict[eng]["perts"] = pert_lig[0]
                values_dict[eng]["ligs"] = pert_lig[1]
                # put results into values dict
                values_dict[eng]["pert_results"] = self.calc_pert_dict[eng]
                
            except Exception as e:
                values_dict[eng]["perts"] = [None]
                values_dict[eng]["ligs"] = [None]
                values_dict[eng]["pert_results"] = [None]
                print(e)
                print(f"could not convert {eng} values for plotting. None will be used. Was earlier analysis okay?")

            try:
                values_dict[eng]["val_results"] = self.calc_val_dict[eng]
            except Exception as e:
                values_dict[eng]["val_results"] = [None]
                print(e)
                print(f"could not convert val {eng} values for plotting. None will be used. Was cinnabar analysis carried out correctly?")

        values_dict["experimental"]["perts"] = self.perturbations
        values_dict["experimental"]["ligs"] = self.ligands
        values_dict["experimental"]["pert_results"] = self.all_exper_pert_dict
        values_dict["experimental"]["val_results"] = self.all_exper_val_dict #normalised data

        self.values_dict = values_dict

        return values_dict

    def _calc_all_dict_to_df(self):
        """for all identified engines, other results, and experimental, convert the dictionary into a dataframe.
        """

        for name in self.names_list:
            self._dict_to_df(x_name=name)

    def _dict_to_df(self, x_name="experimental"):
        """turn a dictionary of results into a dataframe

        Args:
            x_name (str, optional): name of the dictionary to convert. Defaults to "experimental".

        Returns:
            dict: dictionary of pandas dataframe, which is a dictionary (other names) of the 'pv' (pert or val) to match the values found in the x_name dictionary.
        """

        # calculate df w respect to each other value

        x_name = self._validate_in_names_list(x_name)

        freenrg_df_dict = {}

        to_convert_list = [x for x in self.names_list]
        # to_convert_list.remove(x_name)
        for name in to_convert_list:
            freenrg_df_dict.update({name:{}})

        # construct dict with experimental freenrg and error and computed
        for name in to_convert_list: # will do this for engines and other results

            for pv in ["pert","val"]:

                if pv == "pert":
                    which_list = "perts"
                elif pv == "val":
                    which_list = "ligs"

                freenrg_df = self.match_dicts_to_df(self.values_dict[x_name][f"{pv}_results"],
                                                    self.values_dict[name][f"{pv}_results"],
                                                    x_name,
                                                    "calc",
                                                    self.values_dict[x_name][which_list],
                                                    verbose=self._is_verbose)

                freenrg_df_dict[name][pv] = freenrg_df

                # save our results to a file that can be opened in e.g. Excel.
                freenrg_df.to_csv(f"{self.output_folder}/{name}_vs_{x_name}_{pv}_results_table_{self.file_ext}_{self.net_ext}.csv")
        
        self.freenrg_df_dict[x_name] = freenrg_df_dict

        return freenrg_df_dict

    @staticmethod
    def match_dicts_to_df(dict_x, dict_y, x_name, y_name, values=None, verbose=False):
        """match two dictionaries into one dataframe (to be used for plotting)

        Args:
            dict_x (dict): dictionary of values for x, in format value : (freenerg, err)
            dict_y (dict): dictionary of values for y, in format value : (freenerg, err)
            x_name (str): name of the x values ( eg experimental )
            y_name (str): name of the y values
            values (list, optional): list of dict values to convert. Defaults to None.

        Returns:
            _type_: _description_
        """

        freenrg_dict = {}

        if values:
            values = validate.is_list(values)
        else:
            values = dict_x.keys()

        for value in values:
            try:
                x_ddG = dict_x[value][0]
                x_err = dict_x[value][1]
                y_ddG = dict_y[value][0]
                y_err = dict_y[value][1]
                freenrg_dict[value] = [x_ddG, x_err, y_ddG, y_err]
            except Exception as e:
                if verbose:
                    print(f"{value} not in both dicts, {x_name} and {y_name}")
        
        freenrg_df = pd.DataFrame(freenrg_dict, index=[f"freenrg_{x_name}", f"err_{x_name}", f"freenrg_{y_name}", f"err_{y_name}"]).transpose()

        return freenrg_df

    @staticmethod
    def _prune_perturbations(df, perturbations, remove=False):
        """keep or remove perturbations in list from the dataframe

        Args:
            df (pandas.dataframe): dataframe of results
            perturbations (list): list of perturbations that want to keep.
            remove (boolean): whether to keep or remove the perturbations in the list

        Returns:
            df: pruned dataframe
        """
        
        remove = validate.boolean(remove)
        perturbations = validate.is_list(perturbations, make_list=True)

        if remove:
            # keep only specified perturbations in the dataframe
            to_del = []
            for pert in df.index:
                if pert in perturbations:
                    to_del.append(pert)
        else:
            # keep only specified perturbations in the dataframe
            to_del = []
            for pert in df.index:
                if pert not in perturbations:
                    to_del.append(pert)

        for pert in to_del:
            df = df.drop(index=[pert])

        return df


    def _set_style(self):
        """set the style of the graphs.
        """

        self.set_colours()
        self._get_bar_spacing()


    def set_colours(self, colour_dict=None):
        """set the colours of the bars or scatter plots.

        Args:
            colour_dict (dict, optional): dicitonary of names and their colours. Defaults to None.

        Returns:
            dict: dictionary of new colours
        """
        
        set_colour_dict = self._set_colours(colour_dict)
        self.colours = set_colour_dict

        return set_colour_dict

    @staticmethod
    def _set_colours(colour_dict=None):
        """set colours to replace those in the default dictionary.

        Args:
            colour_dict (_type_, optional): _description_. Defaults to None.

        Returns:
            _type_: _description_
        """

        default_colour_dict = {"AMBER":"orange",
                    "SOMD":"darkturquoise",
                    "GROMACS":"orchid",
                    "experimental":"midnightblue"
                    }

        if colour_dict:
            colour_dict = validate.dictionary(colour_dict)
            for key in colour_dict:
                # replace default colour dict keys with those in the passed dictionary
                default_colour_dict[key] = colour_dict[key]
        
        return default_colour_dict


    def _get_bar_spacing(self):
        """_get the bar spacing based on how many engines.
        """

        bar_spacing, bar_width = plotting_engines.get_bar_spacing(names=self.engines)

        self._bar_spacing = bar_spacing
        self._bar_width = bar_width


    @staticmethod
    def get_bar_spacing(names=None):

        names = validate.is_list(names)

        # so experimental is the last thing plotted
        if "experimental" in names:
            names.remove("experimental")
            names.append("experimental")
        
        placement_dict = {}

        if len(names) == 6:
            width = 0.12  # set bar width
            placement = [-width*(5/2), -width*(3/2), -width*(1/2), width*(1/2), width*(3/2), width*(5/2)]
        elif len(names) == 5:
            width = 0.14  # set bar width
            placement = [-width*(4/2), -width*(2/2), 0, width*(2/2), width*(4/2)]
        elif len(names) == 4:
            width = 0.15  # set bar width
            placement = [-width*(3/2), -width*(1/2), width*(1/2), width*(3/2)]
        elif len(names) == 3:
            width = 0.23  # set bar width
            placement = [-width*(2/2), 0, width*(2/2)]
        elif len(names) == 2:
            width = 0.4  # set bar width
            placement = [-width*(1/2), width*(1/2)]
        elif len(names) == 1:
            width = 0.6  # set bar width
            placement = [0]
        else:
            raise ValueError("length of the engine list + exp cannot exceed 6? must have atleast 1 engine/exp.")

        for eng,place in zip(names, placement):
            placement_dict.update({eng:place}) # for each engine

        return placement_dict, width


    def _plotting_engines(self, engines_query):
        """engines to be used for plotting

        Args:
            engines_query (list): list of engines to be plotted for

        Returns:
            list: validated list of engines to be plotted for.
        """

        # get engines for analysis
        if not engines_query:
            engines = self.engines
        else:
            try:
                engines = validate.engines(engines_query)
            except:
                print("engine input not recognised. Will use all engines for which there is data.")
                engines = self.engines
        
        return engines

    def bar(self, pert_val=None, names=None, values=None, **kwargs):
        """plot a bar plot of the results

        Args:
            pert_val (str, optional): whether plotting 'pert' ie perturbations or 'val' ie values (per ligand result). Defaults to None.
            names (list, optional): engines and other results and experimental to plot for. Defaults to None.
            values (list, optional): list of values (perturbations or ligands) to plot for. Defaults to None.
        """

        pert_val = validate.pert_val(pert_val)

        if names:
            names = validate.is_list(names, make_list=True)
            for eng in names:
                self._validate_in_names_list(eng)
            bar_spacing, width = plotting_engines.get_bar_spacing(names=names)
        # if no engines provided, use the defaults that were set based on the analysis object
        else:
            names = self.engines
            bar_spacing = self._bar_spacing
            width = self._bar_width

        if not values:
            if pert_val == "pert":
                values = self.perturbations
            elif pert_val == "val":
                values = self.ligands
        else:
            values = validate.is_list(values)

        plt.rc('font', size=12)
        fig, ax = plt.subplots(figsize=(15,8))

        # df_dict = freenrg_df_plotting

        for eng in names:

            col = self.colours[eng]
            space = bar_spacing[eng]
            
            # just always compare to experimental for this
            freenrg_df_plotting = self.freenrg_df_dict["experimental"][eng][pert_val].fillna(0)

            # prune df to only have perturbations considered
            freenrg_df_plotting = self._prune_perturbations(freenrg_df_plotting, values)

            # determine positions for X axis labels.
            x_locs = np.arange(len(freenrg_df_plotting))

            # plot both our experimental and FEP free energies using an offset on the x position so bars don't overlap.
            ax.bar(x_locs + space, height=freenrg_df_plotting["freenrg_calc"], width=width, yerr=freenrg_df_plotting["err_calc"],
                            label=eng, color=col)

        #plt.xlabel('ΔΔG for experimental (kcal/mol)')
        #plt.ylabel('ΔΔG for calculated (kcal/mol)')
        # format the plot further.
        plt.axhline(color="black")
        plt.title(f"Freenrg for {self.file_ext.replace('_',',')}, {self.net_ext.replace('_',',')}")
        if pert_val == "pert":
            plt.ylabel("$\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
            plt.xlabel("perturbations")
        elif pert_val == "val":
            plt.ylabel("$\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
            plt.xlabel("ligands")
        plt.xticks(x_locs, freenrg_df_plotting.index, rotation=70, ha="right")
        plt.legend()

        eng_name = "_".join(str(eng) for eng in names)
        plt.savefig(f"{self.graph_folder}/barplot_{pert_val}_{self.file_ext}_{self.net_ext}_{eng_name}.png", dpi=300, bbox_inches='tight') # TODO fix eng name
        plt.show()


    def scatter(self, pert_val=None, engines=None, name="experimental", values=None, **kwargs): # TODO change x_name, y_name
        """plot scatter plot.

        Args:
            pert_val (str, optional): whether plotting 'pert' ie perturbations or 'val' ie values (per ligand result). Defaults to None.
            engines (list, optional): engines to plot for. Defaults to None.
            name (str, optional): what to plot against. Defaults to "experimental". #TODO fix so can plot multiple other results
            values (list, optional): list of values (perturbations or ligands) to plot for. Defaults to None.

        Raises:
            ValueError: the name must be available in the other names list ie have results assosciated with it.
        """

        pert_val = validate.pert_val(pert_val)
        name = self._validate_in_names_list(name)

        if engines:
            engines = validate.is_list(engines, make_list=True)
            for eng in engines:
                if eng not in self.names_list:
                    raise ValueError("name must be in calc names")

            # engines = self._plotting_engines(engines)
        # if no engines provided, use the defaults that were set based on the analysis object
        else:
            engines = self.engines

        if not values:
            if pert_val == "pert":
                values = self.perturbations
            elif pert_val == "val":
                values = self.ligands
        else:
            values = validate.is_list(values)

        # plot a scatter plot
        plt.rc('font', size=20)
        fig, ax = plt.subplots(figsize=(7,7))

        lines = []

        for eng in engines:

            col = self.colours[eng]

            freenrg_df_plotting = self.freenrg_df_dict[name][eng][pert_val].dropna()

            # prune df to only have perturbations considered
            freenrg_df_plotting = self._prune_perturbations(freenrg_df_plotting, values)

            x = freenrg_df_plotting[f"freenrg_{name}"]
            y = freenrg_df_plotting["freenrg_calc"]
            x_er = freenrg_df_plotting[f"err_{name}"]
            y_er = freenrg_df_plotting["err_calc"]               

            scatterplot = [plt.scatter(x, y, zorder=10, c=col)]    

            # diff shapes for diff proteins
            # scatterplot = [plt.scatter(x[:4], y[:4], zorder=10, c=col, label="TYK2"),
            #                plt.scatter(x[4:5], y[4:5], zorder=10, c=col, marker="D", label="p38"),
            #                plt.scatter(x[5:], y[5:], zorder=10, c=col, marker="s",label="MCL1")]
            
            lines += plt.plot(0,0,c=col, label=eng)

            plt.errorbar(x , y,
                        yerr=y_er,
                        xerr=x_er,   # comment this line to hide experimental error bars \
                                    # as this can sometimes overcrowd the plot.
                        ls="none",
                        lw=0.5, 
                        capsize=2,
                        color="black",
                        zorder=5
                        )

            #plotting lines - need to change this part so incl from the correct written equations if want a linear fit line
            # x_line = np.linspace(-2,2,20)
            # y_line = (slope)*(x_line) + (intercept)
            # ax.plot(x_line, y_line, label=eng)

        #plotting error bars
        plt.errorbar(x , y,
                    yerr=y_er,
                    # xerr=x_er,   # comment this line to hide experimental error bars \
                                # as this can sometimes overcrowd the plot.
                    ls="none",
                    lw=0.5, 
                    capsize=2,
                    color="black",
                    zorder=5
                    )

        # plot 1/2 kcal bounds:
        plt.fill_between(
                        x=[-100, 100], 
                        y2=[-100.25,99.75],
                        y1=[-99.75, 100.25],
                        lw=0, 
                        zorder=-10,
                        alpha=0.3,
                        color="grey")
        # upper bound:
        plt.fill_between(
                        x=[-100, 100], 
                        y2=[-99.5,100.5],
                        y1=[-99.75, 100.25],
                        lw=0, 
                        zorder=-10,
                        color="grey", 
                        alpha=0.2)
        # lower bound:
        plt.fill_between(
                        x=[-100, 100], 
                        y2=[-100.25,99.75],
                        y1=[-100.5, 99.5],
                        lw=0, 
                        zorder=-10,
                        color="grey", 
                        alpha=0.2)

        min_lim, max_lim = self._get_bounds_scatter(engines, self.freenrg_df_dict[name], pert_val, values, name)

        # for a scatterplot we want the axis ranges to be the same. 
        plt.xlim(min_lim*1.3, max_lim*1.3)
        plt.ylim(min_lim*1.3, max_lim*1.3)

        plt.axhline(color="black", zorder=1)
        plt.axvline(color="black", zorder=1)
        
        # default
        y_label = None
        x_label = None
        title = None
        include_key = True

        # check kwargs incase there is plotting info
        for key,value in kwargs.items():
            if key == "y label":
                y_label = value
            if key == "x label":
                x_label = value
            if key == "title":
                title = value
            if key == "key":
                include_key = validate.boolean(value)
            if key == "save":
                save_fig_location = validate.string(f"{value}.png")

        if include_key:
            labels = [l.get_label() for l in lines]
            plt.legend(lines, labels, loc='upper left')

        if title:
            title = title
        else:
            if name:
                title = f"Computed vs {name}\nfor {self.file_ext.replace('_',',')}, {self.net_ext.replace('_',',')}"
            else:
                title = f"Computed vs Experimental\nfor {self.file_ext.replace('_',',')}, {self.net_ext.replace('_',',')}"

        plt.title(title, fontsize=20)

        if y_label:
            plt.ylabel(f"{y_label}")
        else:
            if pert_val == "pert":
                plt.ylabel("Computed $\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
            elif pert_val == "val":
                plt.ylabel("Computed $\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")

        if x_label:
            plt.xlabel(f"{x_label}")
        else:
            if pert_val == "pert":
                if name:
                    plt.xlabel(f"{name} " + "$\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
                else:
                    plt.xlabel("Experimental $\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
            elif pert_val == "val":
                if name:
                    plt.xlabel(f"{name} " + "$\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
                else:
                    plt.xlabel("Experimental $\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")

        eng_name = self._get_eng_name(engines)
        save_fig_location = f"{self.graph_folder}/calc_vs_{name}_scatterplot_{pert_val}_{self.file_ext}_{self.net_ext}_{eng_name}.png"
        plt.savefig(save_fig_location, dpi=300, bbox_inches='tight')
        plt.show()


    def outlier(self, pert_val="pert", engines=None, outliers=3,  name="experimental", values=None, **kwargs):
        """plot scatter plot with annotated outliers.

        Args:
            pert_val (str, optional): whether plotting 'pert' ie perturbations or 'val' ie values (per ligand result). Defaults to None.
            engines (list, optional): engines to plot for. Defaults to None.
            outliers (int, optional): number of outliers to annotate. Defaults to 3.
            name (str, optional): what to plot against. Defaults to "experimental". #TODO fix so can plot multiple other results
            values (list, optional): list of values (perturbations or ligands) to plot for. Defaults to None.
            
        """

        pert_val = validate.pert_val(pert_val)
        name = self._validate_in_names_list(name)

        if engines:
            engines = self._plotting_engines(engines)
        # if no engines provided, use the defaults that were set based on the analysis object
        else:
            engines = self.engines

        number_outliers_to_annotate = validate.integer(outliers)

        if not values:
            if pert_val == "pert":
                values = self.perturbations
            elif pert_val == "val":
                values = self.ligands
        else:
            values = validate.is_list(values)

        # plot a scatter plot
        plt.rc('font', size=12)
        fig, ax = plt.subplots(figsize=(10,10))

        lines = []

        for eng in engines:

            col = self.colours[eng]

            freenrg_df_plotting = self.freenrg_df_dict[name][eng][pert_val].dropna()
            # prune df to only have perturbations considered
            freenrg_df_plotting = self._prune_perturbations(freenrg_df_plotting, values)

            x = freenrg_df_plotting[f"freenrg_{name}"]
            y = freenrg_df_plotting["freenrg_calc"]
            x_er = freenrg_df_plotting[f"err_{name}"]
            y_er = freenrg_df_plotting["err_calc"] 

            # get an array of the MUE values comparing experimental and FEP values. Take the absolute values.
            mue_values = abs(x - y)
 
            # find the n ligand names that are outliers.
            outlier_names = mue_values.nlargest(number_outliers_to_annotate).index.values.tolist()
            print(outlier_names)

            # construct a list of labels to annotate the scatterplot with.
            annot_labels = []
            colours = []
            for label in freenrg_df_plotting.index.values:
                # if the ligand is an outlier, append the name to the annotation labels list.
                if label in outlier_names:
                    if len(engines) > 1:
                        annot_labels.append(f"{label}, {eng}")
                    else:
                        annot_labels.append(f"{label}")
                    colours.append(self.colours["experimental"])
                else:
                    # if the ligand is not an outlier, append an empty string to the annotation labels list.
                    annot_labels.append("")
                    colours.append(self.colours[eng])

            scatterplot = [plt.scatter(x, y, zorder=10, c=colours)] 
            lines += plt.plot(0,0,c=col, label=eng)
            plt.errorbar(x , y,
                        yerr=y_er,
                        xerr=x_er,   # comment this line to hide experimental error bars \
                                    # as this can sometimes overcrowd the plot.
                        ls="none",
                        lw=0.5, 
                        capsize=2,
                        color="black",
                        zorder=5
                        )

            # then, after generating the figure, we can annotate:
            for i, txt in enumerate(annot_labels):
                plt.annotate(txt, 
                            (freenrg_df_plotting[f"freenrg_{name}"].values.tolist()[i]+0.1,     # x coords
                            freenrg_df_plotting["freenrg_calc"].values.tolist()[i]+0.1),    # y coords
                            size=15, color=self.colours["experimental"])

        # can plot a line for ideal
        # plt.plot((min_lim*1.3,max_lim*1.3),(min_lim*1.3,max_lim*1.3), color="teal")

        # or if want to plot 1/2 kcal bounds
        plt.fill_between(
                        x=[-100, 100], 
                        y2=[-100.25,99.75],
                        y1=[-99.75, 100.25],
                        lw=0, 
                        zorder=-10,
                        alpha=0.3,
                        color="grey")
        # upper bound:
        plt.fill_between(
                        x=[-100, 100], 
                        y2=[-99.5,100.5],
                        y1=[-99.75, 100.25],
                        lw=0, 
                        zorder=-10,
                        color="grey", 
                        alpha=0.2)
        # lower bound:
        plt.fill_between(
                        x=[-100, 100], 
                        y2=[-100.25,99.75],
                        y1=[-100.5, 99.5],
                        lw=0, 
                        zorder=-10,
                        color="grey", 
                        alpha=0.2)

        labels = [l.get_label() for l in lines]
        plt.legend(lines, labels, loc='upper left')
        
        # for a scatterplot we want the axis ranges to be the same.
        min_lim, max_lim = self._get_bounds_scatter(engines, self.freenrg_df_dict[name], pert_val, values, name)
        plt.xlim(min_lim*1.3, max_lim*1.3)
        plt.ylim(min_lim*1.3, max_lim*1.3)

        plt.axhline(color="black", zorder=1)
        plt.axvline(color="black", zorder=1)

        if name:
            plt.title(f"Computed vs {name} outliers\nfor {self.file_ext.replace('_',',')}, {self.net_ext.replace('_',',')}")
        else:
            plt.title(f"Computed vs Experimental outliers\nfor {self.file_ext.replace('_',',')}, {self.net_ext.replace('_',',')}")
        if pert_val == "pert":
            plt.ylabel("$\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
            if name:
                plt.xlabel(f"{name} " + "$\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
            else:
                plt.xlabel("Experimental $\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
        elif pert_val == "val":
            plt.ylabel("$\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
            if name:
                plt.xlabel(f"{name} " + "$\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
            else:
                plt.xlabel("Experimental $\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")

        eng_name = self._get_eng_name(engines)
        plt.savefig(f"{self.graph_folder}/calc_vs_{name}_outlierplot_{pert_val}_{self.file_ext}_{self.net_ext}_{eng_name}.png", dpi=300, bbox_inches='tight')
        plt.show()    

    # some functions so cleaner in plotting functions
    @staticmethod
    def _get_eng_name(engines):
        """get the engine names for writing the titles and file names

        Args:
            engines (list): list of engines

        Returns:
            str: name for the plotting and file name
        """
        #TODO adjust so also includes other results names
        if len(engines) == 3:
            eng_name = "all"
        else:
            eng_name = "_".join(str(eng) for eng in engines)

        return eng_name

    @staticmethod
    def _get_bounds_scatter(engines, freenrg_df_dict, pert_val, values, name):
        """get the upper and lower bounds of the scatter plot based on the results plotted.

        Args:
            engines (list): list of engines and other results.
            freenrg_df_dict (dict): dictionary of dataframes of results
            pert_val (str, optional): whether plotting 'pert' ie perturbations or 'val' ie values (per ligand result).
            values (list): list of values (perturbations or ligands) to plot for. Defaults to None.
            name (str): what is being plotted against.

        Returns:
            tuple: min and max limit for plotting.
        """

        # get the bounds. This can be done with min/max or simply by hand.
        all_freenrg_values_pre = []
        for eng in engines:
            freenrg_df_plotting = freenrg_df_dict[eng][pert_val].dropna()
            freenrg_df_plotting = plotting_engines._prune_perturbations(freenrg_df_plotting, values)
            x = np.array(freenrg_df_plotting[f"freenrg_{name}"]).tolist() # TODO also fix so does based on which df is passed - get default non calc name?
            y = np.array(freenrg_df_plotting["freenrg_calc"]).tolist()
            all_freenrg_values_pre.append(x)
            all_freenrg_values_pre.append(y)

        all_freenrg_values = []
        for sublist in all_freenrg_values_pre:
            for item in sublist:
                all_freenrg_values.append(item)

        min_lim = min(all_freenrg_values)   
        max_lim = max(all_freenrg_values)

        return min_lim, max_lim

    # TODO plot cycle closure things?

    def histogram(self, engines=None, pert_val=None, error_dict=None, file_ext=None, perturbations=None, name="experimental"):
        """default plots histogram of SEM, if error dict supplied in format {engine: error_list}, will plot these

        Args:
            engines (list, optional): list of engines. Defaults to None.
            pert_val (str, optional): whether plotting 'pert' ie perturbations or 'val' ie values (per ligand result). Defaults to None.
            error_dict (dict, optional): dictionary of errors if want to plot that instead. Defaults to None.
            file_ext (str, optional): file extension to be used for the plots. Defaults to None (object file extension).
            perturbations (list, optional): list of perturbations to plot for. Defaults to None.
            name (str): what is being plotted against. Defaults to "experimental".

        Returns:
            dict: dictionary of histograms (for names / engines, and the over all 'dist')
        """

        name = self._validate_in_names_list(name)

        if error_dict:
            type_error = "error"
            if file_ext:
                type_error = validate.string(file_ext)
        else:
            pert_val = validate.pert_val(pert_val)
            if pert_val == "pert":
                type_error = "perturbations"
            if pert_val == "val":
                type_error = "value"
                
        if engines:
            engines = self._plotting_engines(engines)
        # if no engines provided, use the defaults that were set based on the analysis object
        else:
            engines = self.engines

        if not perturbations:
            perturbations = self.perturbations
        else:
            perturbations = validate.is_list(perturbations)

        best_fit_dict = {}
        histogram_dict = {}

        for eng in engines:

            # set plot defaults
            plt.rc('font', size=12)
            fig, ax = plt.subplots(figsize=(8,8))
            col = self.colours[eng]
     
            if error_dict:
                x = error_dict[eng] 
            else:
                freenrg_df_plotting = [self.freenrg_df_dict][name][eng][pert_val].dropna()
                freenrg_df_plotting = self._prune_perturbations(freenrg_df_plotting, perturbations)
                x = freenrg_df_plotting["err_calc"]   

            # no_bins = int(len(freenrg_df_plotting["err_exp"])/8)
            no_bins = 6 # TODO calculate no of bins so min and max per bin

            # Fit a normal distribution to the data, mean and standard deviation
            mu, std = norm.fit(x)
    
            # plot histogram
            plt.hist(x, bins=no_bins, density=True, alpha=0.7, color=col, edgecolor="grey")
            
            # Plot the PDF.
            xmin, xmax = plt.xlim()
            x = np.linspace(xmin, xmax, 100)
            y = norm.pdf(x, mu, std)
            
            plt.plot(x, y, '--', linewidth=2, color=self.colours["experimental"])
            
            best_fit_dict.update({eng:((x, y), mu, std)})
            # TODO also save as pickle? and also save as tipe with mu and std

            #plot
            plt.xlabel('Error')
            plt.ylabel('Frequency')
            plt.title(f"Distribution of error for {type_error}, {eng}, {self.net_ext.replace('_',', ')}\n mu = {mu:.3f} , std = {std:.3f}")
            plt.savefig(f"{self.graph_folder}/fep_vs_exp_histogram_{type_error}_{self.file_ext}_{self.net_ext}_{eng}.png", dpi=300, bbox_inches='tight')
            plt.show()

            # add to histogram dict for shared plotting
            histogram_dict.update({eng: fig})


        # plot the distributions
        fig, ax = plt.subplots(figsize=(10,10))

        lines = []

        for eng in engines:
            col = self.colours[eng]
            plt.plot(best_fit_dict[eng][0][0], best_fit_dict[eng][0][1], 'k', linewidth=2, color=col)
            lines += plt.plot(0,0,c=col, label=eng)

        labels = [l.get_label() for l in lines]
        plt.legend(lines, labels, loc='upper right')

        plt.xlabel('Error')
        plt.ylabel('Frequency')  
        eng_name = self._get_eng_name(engines)
        plt.title(f"Distribution of error for {type_error}, {eng_name}, {self.net_ext.replace('_',', ')}")
        plt.savefig(f"{self.graph_folder}/fep_vs_exp_normal_dist_{type_error}_{self.file_ext}_{self.net_ext}_{eng_name}.png", dpi=300, bbox_inches='tight')
        plt.show()

        histogram_dict.update({"dist": fig})

        return histogram_dict
