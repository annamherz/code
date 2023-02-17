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

    def __init__(self, analysis_object=None, output_folder=None):

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


    def analysis_obj_into_format(self):

        ana_obj = self._analysis_object
        
        # analysis information
        self.engines = sorted(ana_obj.engines)
        self.ligands = ana_obj.ligands
        self.perturbations = ana_obj.perturbations

        # for other results
        self.other_results_names = ana_obj.other_results_names
        
        # file extension
        self._file_ext()
        self._net_ext()

        # dictionaries of engines for plotting from cinnabar
        self.calc_val_dict = ana_obj.cinnabar_calc_val_dict
        self.exper_val_dict = ana_obj.cinnabar_exper_val_dict
        self.calc_pert_dict = ana_obj.cinnabar_calc_pert_dict
        self.exper_pert_dict = ana_obj.cinnabar_exper_pert_dict
        # experimental calculated directly from exp values (for bar)
        self.all_exper_pert_dict = ana_obj.exper_pert_dict
        self.all_exper_val_dict = ana_obj.normalised_exper_val_dict

    def _file_ext(self):

        file_ext = self._analysis_object.file_ext

        if file_ext == ".+":
            self.file_ext = "na"
        else:
            self.file_ext = file_ext
 
        return file_ext
    
    def _net_ext(self):

        net_ext = self._analysis_object.net_ext
        self.net_ext = net_ext    

        return net_ext    
    
    def set_file_ext(self, file_ext):

        self.file_ext = validate.string(file_ext)

    def set_net_ext(self, net_ext):

        self.net_ext = validate.string(net_ext)

    def _analysis_dicts_to_df(self):

        self._overall_dict()
        self._dict_to_df()


    def _overall_dict(self):

        values_dict = {}
        for eng in self.engines:
            values_dict.update({eng:{}})
        values_dict.update({"experimental":{}})
        for name in self.other_results_names:
            values_dict.update({name:{}})

        # run for all engines with selected network and populate the dictionary for plotting
        for eng in self.engines:
            # get perts and ligands for each engine
            pert_lig = get_info_network_from_dict(self.calc_pert_dict[eng])
            values_dict[eng]["perts"] = pert_lig[0]
            values_dict[eng]["ligs"] = pert_lig[1]
            # put results into values dict
            values_dict[eng]["pert_results"] = self.calc_pert_dict[eng]
            values_dict[eng]["val_results"] = self.calc_val_dict[eng]


        values_dict["experimental"]["perts"] = self.perturbations
        values_dict["experimental"]["ligs"] = self.ligands
        values_dict["experimental"]["pert_results"] = self.all_exper_pert_dict
        values_dict["experimental"]["val_results"] = self.all_exper_val_dict #normalised data

        # add other values to dict 
        for name in self.other_results_names:
            # get perts and ligands for each engine
            pert_lig = get_info_network_from_dict(self.calc_pert_dict[name])
            values_dict[name]["perts"] = pert_lig[0]
            values_dict[name]["ligs"] = pert_lig[1]
            # put results into values dict
            values_dict[name]["pert_results"] = self.calc_pert_dict[name]
            values_dict[name]["val_results"] = self.calc_val_dict[name]

        self.values_dict = values_dict

        return values_dict

    def _dict_to_df(self):

        freenrg_df_dict = {}
        to_convert_list = []
        for eng in self.engines:
            freenrg_df_dict.update({eng:{}})
            to_convert_list.append(eng)
        for name in self.other_results_names:
            freenrg_df_dict.update({name:{}})
            to_convert_list.append(name)

        # construct dict with experimental freenrg and error and computed
        for eng in to_convert_list: # will do this for engines and other results

            for pv in ["pert","val"]:

                freenrg_pert_dict = {}

                if pv == "pert":
                    which_list = "perts"
                elif pv == "val":
                    which_list = "ligs"

                for value in self.values_dict[eng][which_list]:
                    exp_ddG = self.values_dict["experimental"][f"{pv}_results"][value][0]
                    exp_err = self.values_dict["experimental"][f"{pv}_results"][value][1]
                    comp_ddG = self.values_dict[eng][f"{pv}_results"][value][0]
                    comp_err = self.values_dict[eng][f"{pv}_results"][value][1]
                    freenrg_pert_dict[value] = [exp_ddG, exp_err, comp_ddG, comp_err]
                freenrg_df = pd.DataFrame(freenrg_pert_dict, index=["freenrg_exp", "err_exp", "freenrg_fep", "err_fep"]).transpose()
                
                freenrg_df_dict[eng][pv] = freenrg_df

                # save our results to a file that can be opened in e.g. Excel.
                freenrg_df.to_csv(f"{self.output_folder}/fep_{pv}_results_table_{self.file_ext}_{self.net_ext}_{eng}.csv")
        
        self.freenrg_df_dict = freenrg_df_dict

        return freenrg_df_dict
    
    def _match_engine_and_other_results(self, name, pv=None):

        pv = validate.pert_val(pv)

        if pv == "pert":
            which_list = "perts"
        elif pv == "val":
            which_list = "ligs"

        if name in self.other_results_names or name in self.engines:
            
            plotting_df_dict = {}

            for eng in self.engines:
                plotting_df_dict.update({eng:{}})

                plotting_dict = {}

                for value in self.values_dict[eng][which_list]:
                    try:
                        exp_ddG = self.values_dict[f"{name}"][f"{pv}_results"][value][0]
                        exp_err = self.values_dict[f"{name}"][f"{pv}_results"][value][1]
                        comp_ddG = self.values_dict[eng][f"{pv}_results"][value][0]
                        comp_err = self.values_dict[eng][f"{pv}_results"][value][1]
                        plotting_dict[value] = [exp_ddG, exp_err, comp_ddG, comp_err]
                    except:
                        warnings.warn(f"{value} is not available in {name}.")

                plotting_df = pd.DataFrame(plotting_dict, index=[f"freenrg_{name}", f"err_{name}", "freenrg_fep", "err_fep"]).transpose()
                plotting_df_dict[eng][pv] = plotting_df

            return plotting_df_dict


    def _set_style(self):

        self.set_colours()
        self._get_bar_spacing()


    def set_colours(self, colour_dict=None):
        
        set_colour_dict = self._set_colours(colour_dict)
        self.colours = set_colour_dict

        return set_colour_dict

    @staticmethod
    def _set_colours(colour_dict=None):

        default_colour_dict = {"AMBER":"orange",
                    "SOMD":"darkturquoise",
                    "GROMACS":"orchid",
                    "experimental":"midnightblue"
                    }
        allowed_keys = ["AMBER","SOMD","GROMACS","experimental"]

        if not colour_dict:
            colour_dict = default_colour_dict

        else:
            colour_dict = validate.dictionary(colour_dict)
            for key in colour_dict:
                if key not in allowed_keys:
                    raise ValueError(f"{colour_dict} may only have the keys in {allowed_keys}.")
                # replace in the default dict and have this as new colour dict
                default_colour_dict[key] = colour_dict[key]
                colour_dict = default_colour_dict
        
        return colour_dict


    def _get_bar_spacing(self):

        bar_spacing, bar_width = plotting_engines.get_bar_spacing(engines=self.engines, experimental=True)

        self._bar_spacing = bar_spacing
        self._bar_width = bar_width


    @staticmethod
    def get_bar_spacing(engines=None, experimental=True):
        
        engines = validate.is_list(engines)
        experimental = validate.boolean(experimental)
        
        placement_dict = {}

        if experimental:
            exp_len = 1
        else:
            exp_len = 0

        if (len(engines) + exp_len) == 4:
            width = 0.15  # set bar width
            placement = [-width*(3/2), -width*(1/2), width*(1/2), width*(3/2)]
        elif (len(engines) + exp_len) == 3:
            width = 0.23  # set bar width
            placement = [-width*(2/2), 0, width*(2/2)]
        elif (len(engines) + exp_len) == 2:
            width = 0.4  # set bar width
            placement = [-width*(1/2), width*(1/2)]
        elif (len(engines) + exp_len) == 1:
            width = 0.6  # set bar width
            placement = [0]
        else:
            raise ValueError("length of the engine list + exp cannot exceed 4? must have atleast 1 engine/exp.")

        for eng,place in zip(engines, placement):
            placement_dict.update({eng:place}) # for each engine
        if experimental:
            placement_dict.update({"experimental":placement[-1]}) # add experimental

        return placement_dict, width


    def _plotting_engines(self, engines_query):

        # get engines for analysis
        if not engines_query:
            engines = self.engines
        else:
            try:
                engines_query = validate.is_list(engines_query)
                val_engines = []
                for engine in engines_query:
                    engine_val = validate.engine(engine)
                    val_engines.append(engine_val)
                engines = val_engines
            # if single engine string, put into list
            except:
                try:
                    engines_query = validate.string(engines_query)
                    engines_query = validate.engine(engines_query)
                    engines = [engines_query]
                except:
                    print("engine input not recognised. Will use all engines for which there is data.")
                    engines = self.engines
        
        return engines

    def bar(self, pert_val=None, engines=None):

        pert_val = validate.pert_val(pert_val)

        if engines:
            engines = self._plotting_engines(engines)
            bar_spacing, width = plotting_engines.get_bar_spacing(engines=engines, experimental=True)
        # if no engines provided, use the defaults that were set based on the analysis object
        else:
            engines = self.engines
            bar_spacing = self._bar_spacing
            width = self._bar_width

        plt.rc('font', size=12)
        fig, ax = plt.subplots(figsize=(15,8))

        for eng in engines:

            col = self.colours[eng]
            space = bar_spacing[eng]
            
            freenrg_df_plotting = self.freenrg_df_dict[eng][pert_val].fillna(0)

            # determine positions for X axis labels.
            x_locs = np.arange(len(freenrg_df_plotting))

            # plot both our experimental and FEP free energies using an offset on the x position so bars don't overlap.
            ax.bar(x_locs + space, height=freenrg_df_plotting["freenrg_fep"], width=width, yerr=freenrg_df_plotting["err_fep"],
                            label=eng, color=col)

        # plot experimental
        ax.bar(x_locs + bar_spacing["experimental"], height=freenrg_df_plotting["freenrg_exp"], width=width, yerr=freenrg_df_plotting["err_exp"],
                        label='Experimental', color=self.colours["experimental"]) 

        #plt.xlabel('ΔΔG for experimental (kcal/mol)')
        #plt.ylabel('ΔΔG for calculated (kcal/mol)')
        # format the plot further.
        plt.axhline(color="black")
        plt.title(f"Computed vs Experimental for {self.file_ext.replace('_',',')}, {self.net_ext.replace('_',',')}")
        if pert_val == "pert":
            plt.ylabel("$\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
            plt.xlabel("perturbations")
        elif pert_val == "val":
            plt.ylabel("$\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
            plt.xlabel("ligands")
        plt.xticks(x_locs, freenrg_df_plotting.index, rotation=70, ha="right")
        plt.legend()

        eng_name = self._get_eng_name(engines)
        plt.savefig(f"{self.graph_folder}/fep_vs_exp_barplot_{pert_val}_{self.file_ext}_{self.net_ext}_{eng_name}.png", dpi=300, bbox_inches='tight')
        plt.show()


    def scatter(self, pert_val=None, engines=None, name=None, **kwargs):

        pert_val = validate.pert_val(pert_val)

        if engines:
            engines = self._plotting_engines(engines)
        # if no engines provided, use the defaults that were set based on the analysis object
        else:
            engines = self.engines

        if name:
            exp_name = name
        else:
            exp_name = "exp"

        # plot a scatter plot
        plt.rc('font', size=12)
        fig, ax = plt.subplots(figsize=(10,10))

        lines = []

        for eng in engines:

            col = self.colours[eng]

            if name:
                plotting_df = self._match_engine_and_other_results(name, pv=pert_val)
                freenrg_df_plotting = plotting_df[eng][pert_val].dropna()
            else:
                freenrg_df_plotting = self.freenrg_df_dict[eng][pert_val].dropna()

            x = freenrg_df_plotting[f"freenrg_{exp_name}"]
            y = freenrg_df_plotting["freenrg_fep"]
            x_er = freenrg_df_plotting[f"err_{exp_name}"]
            y_er = freenrg_df_plotting["err_fep"]               

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

        labels = [l.get_label() for l in lines]
        plt.legend(lines, labels, loc='upper left')
        # plt.legend(scatterplot, ["TYK2","p38","MCL1"])

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

        min_lim, max_lim = self._get_bounds_scatter(engines, self.freenrg_df_dict, pert_val)

        # for a scatterplot we want the axis ranges to be the same. 
        plt.xlim(min_lim*1.3, max_lim*1.3)
        plt.ylim(min_lim*1.3, max_lim*1.3)

        plt.axhline(color="black", zorder=1)
        plt.axvline(color="black", zorder=1)
        
        # default
        y_label = None
        x_label = None
        title = None

        # check kwargs incase there is plotting info
        for key,value in kwargs.items():
            if key == "y label":
                y_label = value
            if key == "x label":
                x_label = value
            if key == "title":
                title = value

        if title:
            plt.title(f"{title}")
        else:
            if name:
                plt.title(f"Computed vs {name}\nfor {self.file_ext.replace('_',',')}, {self.net_ext.replace('_',',')}")
            else:
                plt.title(f"Computed vs Experimental\nfor {self.file_ext.replace('_',',')}, {self.net_ext.replace('_',',')}")
        
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
        plt.savefig(f"{self.graph_folder}/fep_vs_{exp_name}_scatterplot_{pert_val}_{self.file_ext}_{self.net_ext}_{eng_name}.png", dpi=300, bbox_inches='tight')
        plt.show()


    def outlier(self, pert_val="pert", engines=None, outliers=3, name=None):

        pert_val = validate.pert_val(pert_val)

        if engines:
            engines = self._plotting_engines(engines)
        # if no engines provided, use the defaults that were set based on the analysis object
        else:
            engines = self.engines

        number_outliers_to_annotate = validate.integer(outliers)

        if name:
            exp_name = name
        else:
            exp_name = "exp"

        # plot a scatter plot
        plt.rc('font', size=12)
        fig, ax = plt.subplots(figsize=(10,10))

        lines = []

        for eng in engines:

            col = self.colours[eng]

            if name:
                plotting_df = self._match_engine_and_other_results(name, pv=pert_val)
                freenrg_df_plotting = plotting_df[eng][pert_val].dropna()
            else:
                freenrg_df_plotting = self.freenrg_df_dict[eng][pert_val].dropna()

            x = freenrg_df_plotting[f"freenrg_{exp_name}"]
            y = freenrg_df_plotting["freenrg_fep"]
            x_er = freenrg_df_plotting[f"err_{exp_name}"]
            y_er = freenrg_df_plotting["err_fep"] 

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
                            (freenrg_df_plotting[f"freenrg_{exp_name}"].values.tolist()[i]+0.1,     # x coords
                            freenrg_df_plotting["freenrg_fep"].values.tolist()[i]+0.1),    # y coords
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
        min_lim, max_lim = self._get_bounds_scatter(engines, self.freenrg_df_dict, pert_val)
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
        plt.savefig(f"{self.graph_folder}/fep_vs_{exp_name}_outlierplot_{pert_val}_{self.file_ext}_{self.net_ext}_{eng_name}.png", dpi=300, bbox_inches='tight')
        plt.show()    

    # some functions so cleaner in plotting functions
    @staticmethod
    def _get_eng_name(engines):

        if len(engines) == 3:
            eng_name = "all"
        else:
            eng_name = "_".join(str(eng) for eng in engines)

        return eng_name
    
    # TODO method to get the

    @staticmethod
    def _get_bounds_scatter(engines, freenrg_df_dict, pert_val):

        # get the bounds. This can be done with min/max or simply by hand.
        all_freenrg_values_pre = []
        for eng in engines:
            freenrg_df_plotting = freenrg_df_dict[eng][pert_val].dropna()
            x = np.array(freenrg_df_plotting["freenrg_exp"]).tolist()
            y = np.array(freenrg_df_plotting["freenrg_fep"]).tolist()
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

    def histogram(self, engines=None, pert_val=None, error_dict=None, file_ext=None):
        """default plots histogram of SEM, if error dict supplied in format {engine: error_list}, will plot these
        """
        # TODO this is currently plotting the standard error of the 

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
                freenrg_df_plotting = self.freenrg_df_dict[eng][pert_val].dropna()
                x = freenrg_df_plotting["err_fep"]   

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
            # TODO also save as pickle, and also save as tipe with mu and std

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


    def calc_mae(self, pert_val=None):
        # calc mae for a provided dictionary in the format wanted

        pv = validate.pert_val(pert_val)

        values_dict = self.values_dict
        engines = self.engines

        mae_pert_df = pd.DataFrame(columns=engines,index=engines)
        mae_pert_df_err = pd.DataFrame(columns=engines,index=engines)

        # iterate over all possible combinations
        for combo in it.product(engines, engines):
            eng1 = combo[0]
            eng2 = combo[1]

            eng1_vals = []
            eng2_vals = []

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