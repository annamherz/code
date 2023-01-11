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

class plotting():

    def __init__(self, analysis_object):

        self._analysis_object = analysis_object

        plotting.colours()

        # get info about things for plotting from analysis

        # get number of engines - if multiple

        

    def colours(self, colour_dict=None):

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

    def overall_dict():

        # good idea to check the network considering etc.
        # can use the perts obtained from the get_info_network, and then plot this.
        # can collate all perts and ligands for all engines as well

        all_perturbations = []
        all_ligands = []

        # run for all engines with selected network and populate the dictionary for plotting
        for eng in engines:
            results_all_files = glob.glob(f"{out_folder}/repeat_*_{eng}_{file_ext}.csv")

            # first get network info for each engine
            values = get_info_network(eng, results_files=results_all_files,
                                        net_file=net_file,
                                        output_folder=temp_folder)
            values_dict[eng]["results_files"] = values[2]
            values_dict[eng]["perts"] = values[0]
            values_dict[eng]["ligs"] = values[1]

            for pert in values[0]:
                if pert not in all_perturbations:
                    all_perturbations.append(pert)
            for lig in values[1]:
                if lig not in all_ligands:
                    all_ligands.append(lig)

            # make a dict of the computed results
            comp_dict = make_dict_comp_results(results_files=values_dict[eng]["results_files"], 
                                                perturbations=values_dict[eng]["perts"],
                                                file_name=f"{res_folder}/{comp_pert_file_name}_{eng}",
                                                engine=eng)
            values_dict[eng]["pert_results"] = comp_dict

        # also want to save exp results for plotting and comparison

        # first need to convert the yml file into one useable by freenergworkflows
        convert_yml_into_freenrgworkflows(exp_file, exp_file_dat)
        # TODO diff exp conversion that doesnt rely on freenrgworkflows
        # using freenergworkflows to convert
        experiments = experiments.ExperimentalData()
        experiments.compute_affinities(exp_file_dat, data_type="IC50", comments="#", delimiter=",")
        experimental_DDGs = experiments.freeEnergiesInKcal

        exp_pert_dict,exp_lig_dict = freenrgworkflows_into_dict(experimental_DDGs, all_ligands, all_perturbations)

        values_dict["experimental"]["perts"] = all_perturbations
        values_dict["experimental"]["ligs"] = all_ligands
        values_dict["experimental"]["pert_results"] = exp_pert_dict
        values_dict["experimental"]["val_results"] = exp_lig_dict
        values_dict["experimental"]["freenrgworkflows_ouput"] = experimental_DDGs




    def convergence():
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


    def dict_to_df():

        # construct dict with experimental freenrg and error and computed
        for eng in engines:
            freenrg_pert_dict = {}
            for pert in values_dict[eng]["perts"]:
                exp_ddG = values_dict["experimental"]["pert_results"][pert][0]
                exp_err = values_dict["experimental"]["pert_results"][pert][1]
                comp_ddG = values_dict[eng]["pert_results"][pert][0]
                comp_err = values_dict[eng]["pert_results"][pert][1]
                freenrg_pert_dict[pert] = [exp_ddG, exp_err, comp_ddG, comp_err]
            freenrg_df_pert = pd.DataFrame(freenrg_pert_dict, index=["freenrg_exp", "err_exp", "freenrg_fep", "err_fep"]).transpose()
            values_dict[eng]["freenrg_df_pert"] = freenrg_df_pert

            # save our results to a file that can be opened in e.g. Excel.
            freenrg_df_pert.to_csv(f"{res_folder}/fep_diff_results_table_{file_ext_out}_{eng}.csv")


    def mae_df_make():
        # TODO funcion for this in dictionaries?

        mae_pert_df, mae_pert_df_err = calc_mae(values_dict, "perts")

        print(mae_pert_df)
        print(mae_pert_df_err)

        mae_pert_df.to_csv(f"{res_folder}/mae_pert_{file_ext_out}.csv", sep=" ")
        mae_pert_df_err.to_csv(f"{res_folder}/mae_pert_err_{file_ext_out}.csv", sep=" ")


    def one_scatter():

        for engine in engines:
            print(engine)
            # check the outputed data frame - if want to get rid of the NaN values entirely can do
            freenrg_df_pert_plotting_scatter = values_dict[engine]["freenrg_df_pert"].dropna()

            # plot a scatter plot
            plt.rc('font', size=12)
            fig, ax = plt.subplots(figsize=(10,10))

            # get these based on which column the data is in.
            x = freenrg_df_pert_plotting_scatter["freenrg_exp"]
            y = freenrg_df_pert_plotting_scatter["freenrg_fep"]
            x_er = freenrg_df_pert_plotting_scatter["err_exp"]
            y_er = freenrg_df_pert_plotting_scatter["err_fep"]

            # plotting the scatterplot
            scatterplot = [plt.scatter(x, y, zorder=10)]

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
            all_freenrg_values = np.concatenate([freenrg_df_pert_plotting_scatter["freenrg_exp"].values,freenrg_df_pert_plotting_scatter["freenrg_fep"].values])
            min_lim = min(all_freenrg_values)   
            max_lim = max(all_freenrg_values)

            # for a scatterplot we want the axis ranges to be the same. 
            plt.xlim(min_lim*1.3, max_lim*1.3)
            plt.ylim(min_lim*1.3, max_lim*1.3)

            plt.axhline(color="black", zorder=1)
            plt.axvline(color="black", zorder=1)

            plt.title(f"Computed vs Experimental with {engine} and {file_ext_out.replace('_',',')}")
            plt.ylabel("Computed $\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
            plt.xlabel("Experimental $\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")

            plt.savefig(f"{res_folder}/fep_vs_exp_scatterplot_pert_{file_ext_out}_{engine}.png", dpi=300, bbox_inches='tight')
            plt.show()


    def one_bar():

        for engine in engines:
            # else, if want to substitute w 0 if it is the experimental values, eg for bar charts
            # THIS DOES NOT WORK FOR SCATTER PLOTS - they will all be at 0, this is only useable for bar plots
            freenrg_df_pert_plotting_bar = values_dict[engine]["freenrg_df_pert"].fillna(0)

            # initiate an empty figure with fixed dimensions.
            fig, ax = plt.subplots(figsize=(15,8))

            # determine positions for X axis labels.
            x_locs = np.arange(len(freenrg_df_pert_plotting_bar))

            # set bar width
            width = 0.35  

            # plot both our experimental and FEP free energies using an offset on the x position so bars don't overlap.
            ax.bar(x_locs - width/2, height=freenrg_df_pert_plotting_bar["freenrg_exp"], width=width, yerr=freenrg_df_pert_plotting_bar["err_exp"],
                            label='Experimental')
            ax.bar(x_locs + width/2, height=freenrg_df_pert_plotting_bar["freenrg_fep"], width=width, yerr=freenrg_df_pert_plotting_bar["err_fep"],
                            label='FEP')
            
            # format the plot further.
            plt.axhline(color="black")
            plt.title(f"Computed vs Experimental with {engine} and {file_ext_out.replace('_',',')}")
            plt.ylabel("$\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
            plt.xticks(x_locs, freenrg_df_pert_plotting_bar.index, rotation=70, ha="right")
            plt.legend()

            plt.savefig(f"{res_folder}/fep_vs_exp_barplot_pert_{file_ext_out}_{engine}.png", dpi=300, bbox_inches='tight')
            plt.show()

# TODO one bar and one scatter can also go into all, just using one engine


    def all_scatter():

        plt.rc('font', size=12)
        fig, ax = plt.subplots(figsize=(8,8))

        colour = ['orange','orchid','darkturquoise','midnightblue']
        lines = []

        for eng,col in zip(engines,colour):

            freenrg_df_pert_plotting_scatter = values_dict[eng]["freenrg_df_pert"].dropna()
            x = freenrg_df_pert_plotting_scatter["freenrg_exp"]
            y = freenrg_df_pert_plotting_scatter["freenrg_fep"]
            x_er = freenrg_df_pert_plotting_scatter["err_exp"]
            y_er = freenrg_df_pert_plotting_scatter["err_fep"]    

            scatterplot = [plt.scatter(x, y, zorder=10, c=col)]    

            # diff shapes for diff proteins
            # scatterplot = [plt.scatter(x[:4], y[:4], zorder=10, c=col, label="TYK2"),
            #                plt.scatter(x[4:5], y[4:5], zorder=10, c=col, marker="D", label="p38"),
            #                plt.scatter(x[5:], y[5:], zorder=10, c=col, marker="s",label="MCL1")]
            lines += plt.plot(0,0,c=col, label=eng)

            # x_er = np.array((eng_dict_plot_exp_for_eng[eng]).iloc[1:len(eng_dict_plot_exp_for_eng[eng]), 2])
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

            #plotting lines - need to change this part so incl from the correct written equations if want a linear fit line
            # x_line = np.linspace(-2,2,20)
            # y_line = (slope)*(x_line) + (intercept)
            # ax.plot(x_line, y_line, label=eng)

        labels = [l.get_label() for l in lines]
        plt.legend(lines, labels, loc='upper left')
        # plt.legend(scatterplot, ["TYK2","p38","MCL1"])

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
            x = np.array(freenrg_df_pert_plotting_scatter["freenrg_exp"]).tolist()
            y = np.array(freenrg_df_pert_plotting_scatter["freenrg_fep"]).tolist()
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
        # plt.title(f"Computed vs Experimental FEP for {protein.upper()}, {file_ext_out.replace('_', ',')}")
        plt.ylabel("Computed $\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
        plt.xlabel("Experimental $\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")

        # plt.legend()

        plt.savefig(f"{res_folder}/fep_vs_exp_scatterplot_pert_{file_ext_out}_all.png", dpi=300, bbox_inches='tight')
        plt.show()


    def all_bar():

        plt.rc('font', size=12)
        fig, ax = plt.subplots(figsize=(15,8))

        colour = ['darkturquoise','orange','orchid','midnightblue']
        # set bar width
        width = 0.15  
        placement = [-width*(3/2), -width*(1/2), width*(1/2), width*(3/2)]

        for eng,col,place in zip(engines,colour,placement):

            freenrg_df_pert_plotting_bar = values_dict[eng]["freenrg_df_pert"].fillna(0)
            x = freenrg_df_pert_plotting_bar["freenrg_exp"]
            y = freenrg_df_pert_plotting_bar["freenrg_fep"]
            x_er = freenrg_df_pert_plotting_bar["err_exp"]
            y_er = freenrg_df_pert_plotting_bar["err_fep"]    

            # determine positions for X axis labels.
            x_locs = np.arange(len(freenrg_df_pert_plotting_bar))

            # plot both our experimental and FEP free energies using an offset on the x position so bars don't overlap.

            ax.bar(x_locs + place, height=freenrg_df_pert_plotting_bar["freenrg_fep"], width=width, yerr=freenrg_df_pert_plotting_bar["err_fep"],
                            label=eng, color=col)

        # plot experimental
        ax.bar(x_locs + placement[-1], height=freenrg_df_pert_plotting_bar["freenrg_exp"], width=width, yerr=freenrg_df_pert_plotting_bar["err_exp"],
                        label='Experimental', color=colour[-1]) 

        #plt.xlabel('ΔΔG for experimental (kcal/mol)')
        #plt.ylabel('ΔΔG for calculated (kcal/mol)')
        # format the plot further.
        plt.axhline(color="black")
        plt.title(f"Computed vs Experimental for {protein.upper()} and {file_ext_out.replace('_',',')}")
        plt.ylabel("$\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
        plt.xticks(x_locs, freenrg_df_pert_plotting_bar.index, rotation=70, ha="right")
        plt.legend()

        plt.savefig(f"{res_folder}/fep_vs_exp_barplot_pert_{file_ext_out}_all.png", dpi=300, bbox_inches='tight')
        plt.show()

    def _bar_spacing():
        pass

    def outlier():
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


    # TODO plot cinnabar, someway to get the stats things from cinnabar
    # can plot diff data series cinnabar?
    # plot cycle closure things?
    # remove all freenrg workflow things
    # rewrite exp calculation