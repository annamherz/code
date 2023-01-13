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
from rdkit import Chem

from ..utils import *
from ._network import *

class plotting_engines():

    def __init__(self, analysis_object=None, res_folder=None):

        if analysis_object:
            self._analysis_object = analysis_object
            # get info about things for plotting from analysis
            self.analysis_obj_into_format()
        else:
            raise ValueError("please provide an analysis object to be plotted for.")

        # place to write results to
        if not res_folder:
            # want to write to the graph directory
            self.results_folder = self._analysis_object.graph_dir
        else:
            self.results_folder = validate.folder_path(res_folder, create=True)

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
        
        # file extension
        self._file_ext()

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
    
    def set_file_ext(self, file_ext):

        self.file_ext = validate.string(file_ext)


    def _analysis_dicts_to_df(self):

        values_dict = self._overall_dict()
        freenrg_df_dict = self._dict_to_df()


    def _overall_dict(self):

        values_dict = {}
        for eng in self.engines:
            values_dict.update({eng:{}})
        values_dict.update({"experimental":{}})

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

        self.values_dict = values_dict

        return values_dict

    def _dict_to_df(self):

        freenrg_df_dict = {}
        for eng in self.engines:
            freenrg_df_dict.update({eng:{}})

        # construct dict with experimental freenrg and error and computed
        for eng in self.engines:

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
                freenrg_df.to_csv(f"{self.results_folder}/fep_{pv}_results_table_{self.file_ext}_{eng}.csv")
        
        self.freenrg_df_dict = freenrg_df_dict

        return freenrg_df_dict

    def _set_style(self):

        self.set_colours()
        self._get_bar_spacing()


    def set_colours(self, colour_dict=None):

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

        self.colours = colour_dict
        
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

        pert_val = validate.string(pert_val)
        if pert_val not in ["pert","val"]:
            raise ValueError("pert_val must be either 'pert' or 'val'")

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
        plt.title(f"Computed vs Experimental for {self.file_ext.replace('_',',')}")
        if pert_val == "pert":
            plt.ylabel("$\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
            plt.xlabel("perturbations")
        elif pert_val == "val":
            plt.ylabel("$\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
            plt.xlabel("ligands")
        plt.xticks(x_locs, freenrg_df_plotting.index, rotation=70, ha="right")
        plt.legend()
        if len(engines) == 3:
            eng_name = "all"
        else:
            eng_name = "_".join(str(eng) for eng in engines)
        plt.savefig(f"{self.results_folder}/fep_vs_exp_barplot_{pert_val}_{self.file_ext}_{eng_name}.png", dpi=300, bbox_inches='tight')
        plt.show()


    def scatter(self, pert_val=None, engines=None):

        pert_val = validate.string(pert_val)
        if pert_val not in ["pert","val"]:
            raise ValueError("pert_val must be either 'pert' or 'val'")

        if engines:
            engines = self._plotting_engines(engines)
        # if no engines provided, use the defaults that were set based on the analysis object
        else:
            engines = self.engines

        # plot a scatter plot
        plt.rc('font', size=12)
        fig, ax = plt.subplots(figsize=(10,10))

        lines = []

        for eng in engines:

            col = self.colours[eng]

            freenrg_df_plotting = self.freenrg_df_dict[eng][pert_val].dropna()
            x = freenrg_df_plotting["freenrg_exp"]
            y = freenrg_df_plotting["freenrg_fep"]
            x_er = freenrg_df_plotting["err_exp"]
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

        # get the bounds. This can be done with min/max or simply by hand.
        all_freenrg_values_pre = []
        for eng in engines:
            freenrg_df_plotting = self.freenrg_df_dict[eng][pert_val].dropna()
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

        # for a scatterplot we want the axis ranges to be the same. 
        plt.xlim(min_lim*1.3, max_lim*1.3)
        plt.ylim(min_lim*1.3, max_lim*1.3)

        plt.axhline(color="black", zorder=1)
        plt.axvline(color="black", zorder=1)

        #plt.xlabel('ΔΔG for experimental (kcal/mol)')
        #plt.ylabel('ΔΔG for calculated (kcal/mol)')
        plt.title(f"Computed vs Experimental for {self.file_ext.replace('_',',')}")
        if pert_val == "pert":
            plt.ylabel("$\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
            plt.xlabel("Experimental $\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
        elif pert_val == "val":
            plt.ylabel("$\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
            plt.xlabel("Experimental $\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")

        if len(engines) == 3:
            eng_name = "all"
        else:
            eng_name = "_".join(str(eng) for eng in engines)

        plt.savefig(f"{self.results_folder}/fep_vs_exp_scatterplot_{pert_val}_{self.file_ext}_{eng_name}.png", dpi=300, bbox_inches='tight')
        plt.show()


    def outlier(self):
        number_outliers_to_annotate = 5
        engine = "GROMACS"

        freenrg_df_plotting_scatter = values_dict[engine]["freenrg_df_pert"].dropna()

        x = freenrg_df_plotting_scatter["freenrg_exp"]
        y = freenrg_df_plotting_scatter["freenrg_fep"]
        x_er = freenrg_df_plotting_scatter["err_exp"]
        y_er = freenrg_df_plotting_scatter["err_fep"]    

        # get an array of the MUE values comparing experimental and FEP values. Take the absolute values.
        mue_values = abs(freenrg_df_plotting_scatter["freenrg_exp"] - freenrg_df_plotting_scatter["freenrg_fep"])

        # find the n ligand names that are outliers.
        outlier_names = mue_values.nlargest(number_outliers_to_annotate).index.values.tolist()
        print(outlier_names)

        # construct a list of labels to annotate the scatterplot with.
        annot_labels = []
        colours = []
        for ligand in freenrg_df_plotting_scatter.index.values:
            # if the ligand is an outlier, append the name to the annotation labels list.
            if ligand in outlier_names:
                annot_labels.append(ligand)
                colours.append("hotpink")
            else:
                # if the ligand is not an outlier, append an empty string to the annotation labels list.
                annot_labels.append("")
                colours.append("teal")

        # Create the same scatterplot as above. Can include some more of the formatting if needed.
        plt.rc('font', size=12)
        plt.figure(figsize=(10,10))

        plt.scatter(x,y, zorder=10, c=colours)

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

        # get the bounds. This can be done with min/max or simply by hand.
        all_freenrg_values = np.concatenate([x.values,y.values])
        min_lim = min(all_freenrg_values)   
        max_lim = max(all_freenrg_values)

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

        # for a scatterplot we want the axis ranges to be the same. 
        plt.xlim(min_lim*1.3, max_lim*1.3)
        plt.ylim(min_lim*1.3, max_lim*1.3)

        plt.axhline(color="black", zorder=1)
        plt.axvline(color="black", zorder=1)

        plt.ylabel("Computed $\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
        plt.xlabel("Experimental $\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")

        # then, after generating the figure, we can annotate:
        for i, txt in enumerate(annot_labels):
            plt.annotate(txt, 
                        (freenrg_df_pert_plotting_scatter["freenrg_exp"].values.tolist()[i]+0.1,     # x coords
                        freenrg_df_pert_plotting_scatter["freenrg_fep"].values.tolist()[i]+0.1),    # y coords
                        size=15, color="hotpink")

        # plt.savefig("analysis/fep_vs_exp_outlier_plot.png", dpi=300, bbox_inches='tight')
        plt.show()    



    def mae_df_make():
        # TODO funcion for this in dictionaries?

        mae_pert_df, mae_pert_df_err = calc_mae(values_dict, "perts")

        print(mae_pert_df)
        print(mae_pert_df_err)

        mae_pert_df.to_csv(f"{res_folder}/mae_pert_{self.file_ext}.csv", sep=" ")
        mae_pert_df_err.to_csv(f"{res_folder}/mae_pert_err_{self.file_ext}.csv", sep=" ")


    # can plot diff data series cinnabar?
    # plot cycle closure things?



def plot_convergence():
    engine = ['AMBER','SOMD','GROMACS']
    # trans = ['ejm42~ejm31','ejm42~ejm55','ejm54~ejm42','ejm55~ejm54','2w~2z','67~60','60~63']
    trans = ['ejm55~ejm54']

    colour = ['orange','orchid','darkturquoise','midnightblue']
    colour_dict = {"AMBER":"orange","SOMD":"darkturquoise","GROMACS":"orchid","experimental":"midnightblue"}
    # plot the convergence w time
    for tra in trans:
        prot = "tyk2"

        for leg in [ 'free','bound']:
            plt.figure()
            lines = []
            for eng,col in zip(engine,colour):
                # with open (f'/home/anna/Documents/benchmark/{prot}/outputs/{eng}/{tra}/{leg}_pmf_{tra}_{eng}.pickle', 'rb') as handle:
                with open (f'/home/anna/Documents/benchmark/{prot}/pickles/{leg}_pmf_{tra}_{eng}.pickle', 'rb') as handle:
                    pmf_dict = pickle.load(handle)
                lines += plt.plot(0,0,c=col, label=eng)
                for repeat in pmf_dict:
                    pmf = pmf_dict[repeat]
                    x =[]
                    y=[]
                    yerr = []
                    for p in pmf:
                        x.append(p[0])
                        y.append(p[1]*(1/BSS.Units.Energy.kcal_per_mol))
                        yerr.append(p[2]*(1/BSS.Units.Energy.kcal_per_mol))
                    plt.errorbar(x,y,yerr=yerr,color=col, ecolor='black')
            plt.xlim(xmin=0,xmax=1)
            plt.ylabel("Computed $\Delta$G$_{transformation}$ / kcal$\cdot$mol$^{-1}$")
            plt.xlabel("Lambda")
            labels = [l.get_label() for l in lines]
            plt.legend(lines, labels)
            plt.title(f"Convergence, {leg} for {tra}")
            plt.savefig(f'/home/anna/Documents/benchmark/{prot}/outputs/{eng}/{tra}/convergence_{leg}.png')

        # plotting delta delta G

        plt.figure()
        lines = []
        for eng in engine:
            # with open (f'/home/anna/Documents/benchmark/{prot}/outputs/{eng}/{tra}/bound_pmf_{tra}_{eng}.pickle', 'rb') as handle:
            with open (f'/home/anna/Documents/benchmark/{prot}/pickles/bound_pmf_{tra}_{eng}.pickle', 'rb') as handle:
                bound_pmf_dict = pickle.load(handle)
            # with open (f'/home/anna/Documents/benchmark/{prot}/outputs/{eng}/{tra}/free_pmf_{tra}_{eng}.pickle', 'rb') as handle:
            with open (f'/home/anna/Documents/benchmark/{prot}/pickles/free_pmf_{tra}_{eng}.pickle', 'rb') as handle:
                free_pmf_dict = pickle.load(handle)
            lines += plt.plot(0,0,c=colour_dict[eng], label=eng)
            for repf,repb in zip(free_pmf_dict,bound_pmf_dict):
                bound_pmf = bound_pmf_dict[repb]
                free_pmf = free_pmf_dict[repf]
                x = []
                y = []
                yerr = []
                for pb,pf in zip(bound_pmf,free_pmf):
                    x.append(pb[0])
                    y.append((pb[1]*(1/BSS.Units.Energy.kcal_per_mol))-(pf[1]*(1/BSS.Units.Energy.kcal_per_mol)))
                    yerr.append((pb[2]*(1/BSS.Units.Energy.kcal_per_mol))+(pf[2]*(1/BSS.Units.Energy.kcal_per_mol)))
                plt.errorbar(x,y,yerr=yerr,color=colour_dict[eng])
        plt.xlim(xmin=0,xmax=1)
        plt.ylabel("Computed $\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
        plt.xlabel("Lambda")
        labels = [l.get_label() for l in lines]
        plt.legend(lines, labels)
        plt.title(f"Convergence for {tra}")
        plt.savefig(f'/home/anna/Documents/benchmark/{prot}/outputs/{eng}/{tra}/convergence_deltadeltaG.png')

