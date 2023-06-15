#!/usr/bin/python3

import BioSimSpace as BSS
from BioSimSpace import Units as _Units
import os
import sys
from argparse import ArgumentParser

print("adding code to the pythonpath...")
code = '/home/anna/Documents/code/python'
if code not in sys.path:
    sys.path.insert(1, code)
import pipeline

from pipeline.prep import *
from pipeline.utils import *
from pipeline.analysis import *

import pandas as pd
import matplotlib.pyplot as plt

def truncation_check(path_to_dir, estimator, start_end="start", statsineff=False, eq=False):

    start_end = validate.truncate_keep(start_end)

    # analyse the work dir
    analysed_pert = analyse(path_to_dir)

    truncate_percentage = [5, 10, 15, 20, 25, 30, 40, 50, 60, 80, 100]
    
    results_dict = {}
    bound_dict = {}
    free_dict = {}

    for trunc_per in truncate_percentage:

        analysed_pert.set_options({ "estimator" : estimator,
                                    "truncate percentage": trunc_per,
                                    "try pickle": True,
                                    "truncate keep":start_end,
                                    "statistical inefficiency":statsineff,
                                    "auto equilibration": eq
                                    })

        analysed_pert.analyse_all_repeats()

        # save output in the pickle dir 
        write_analysis_file(analysed_pert, analysed_pert._pickle_dir, method=trunc_per)

        results_dict[trunc_per] = (analysed_pert.freenrg_val, analysed_pert.freenrg_err)
        bound_dict[trunc_per] = (analysed_pert.bound_val, analysed_pert.bound_err)
        free_dict[trunc_per] = (analysed_pert.free_val, analysed_pert.free_err)

    return results_dict, bound_dict, free_dict


def pert_dict_into_df(pert_dict, plot_error=False, plot_difference=True):

    df = pd.DataFrame.from_dict(pert_dict)
    perts = list(df.columns)
    df = df.reset_index().dropna()

    index_dict = {}
    for x in df['index']:
        index_dict[x] = []

    for pert in perts:
        x_vals = []
        y_vals = []
        for a,x in zip(df[pert], df['index']):
            if plot_error:
                ind = 1
            else:
                ind = 0
            if plot_difference:
                y_val = (df.iloc[-1][pert][ind] - a[ind])
            else:
                y_val = a[ind]
            if not index_dict[x]:
                index_dict[x] = [y_val]
            else:
                index_dict[x].append(y_val)
            y_vals.append(y_val)
            x_vals.append(x)

    for x in df['index']:
        try:
            val_list = [x for x in index_dict[x] if pd.notna(x)]
            avg = np.mean(val_list)
            min_val = min(val_list)
            max_val = max(val_list)
        except:
            avg = None
            min_val = None
            max_val = None

        index_dict[x] = (avg, min_val, max_val)

    df = pd.DataFrame.from_dict(index_dict, orient="index", columns=["avg","min","max"])
    df = df.dropna()
    
    return df

def single_pert_dict_into_df(pert_dict):

    df = pd.DataFrame.from_dict(pert_dict)
    df = df.transpose()
    df.columns = ["avg","max"]
    df["min"] = df.loc[:, "max"]*-1
    
    return df

def plot_truncated(sdf, edf, file_path=None, plot_error=False, plot_difference=True):

    include_key = True

    plt.rc('font', size=12)
    plt.rcParams['axes.xmargin'] = 0 # plt.margins(x=0)
    fig, ax = plt.subplots(figsize=(10,10))
    lines = []

    # fill in final value and its error
    plt.axhline(y=sdf['avg'].iloc[-1], color="c")
    plt.fill_between(sdf.index, sdf['avg'].iloc[-1]+sdf['min'].iloc[-1], sdf['avg'].iloc[-1]+sdf['max'].iloc[-1],  color="paleturquoise", alpha=1)
    lines += plt.plot(0,0,c="c", label="final estimate")
    scatterplot = [plt.plot(sdf.index, sdf['avg'], c="lightcoral")]
    plt.fill_between(sdf.index, sdf['avg']+sdf['min'], sdf['avg']+sdf['max'], color="mistyrose", alpha=.4)
    lines += plt.plot(0,0,c="lightcoral", label="forward")
    scatterplot = [plt.plot(edf.index, edf['avg'], c="cornflowerblue")]
    plt.fill_between(edf.index, edf['avg']+edf['min'], edf['avg']+edf['max'], color="lightskyblue", alpha=.3)
    lines += plt.plot(0,0,c="cornflowerblue", label="reverse")

    if include_key:
        labels = [l.get_label() for l in lines]
        plt.legend(lines, labels, loc='upper left')
    
    plt.xlabel('Percentage of run used')
    if file_path:
        plt.title(f"{file_path.split('/')[-1].split('.')[0].replace('_',' ')}")
    else:
        pass

    if plot_error:
        if plot_difference:
            plt.ylabel("Computed Error difference to final / kcal$\cdot$mol$^{-1}$")
        else:
            plt.ylabel("Computed Error / kcal$\cdot$mol$^{-1}$")
    else:
        if plot_difference:
            plt.ylabel("difference to final result for computed $\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
        else:
            plt.ylabel("Computed $\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")

    if file_path:
        plt.savefig(file_path)

def plot_pert():
    pass

def check_arguments(args):

    # pass the checks to the other check functions
    if args.main_folder:
        main_folder = args.main_folder
    else:
        main_folder = str(input("what is the main folder where all the files should go?: ")).strip()

    if args.network_file:
        network_file = args.network_file
    else:
        network_file = str(input("what is the path to the network file?: ").strip())

    if args.estimator:
        estimator = validate.estimator(args.estimator)
    else:
        estimator = validate.estimator(str(input("what is the analysis method (MBAR/TI) ? : ").strip()))

    if args.stats:
        statsineff = validate.boolean(args.stats)
    else:
        statsineff = False

    if args.autoeq:
        eq = validate.boolean(args.autoeq)
    else:
        eq = False

    return main_folder, network_file, estimator, statsineff, eq

def main():

    # accept all options as arguments
    parser = ArgumentParser(description="run the ligprep")
    parser.add_argument("-mf", "--main_folder", type=str, default=None, help="main folder path to create for all the runs")
    parser.add_argument("-n", "--network_file", type=str, default=None, help="path to network file")
    parser.add_argument("-e", "--estimator", type=str, default="MBAR", help="analysis method (MBAR or TI)")
    parser.add_argument("-s", "--stats", type=str, default=False, help="if using statistical inefficiency")
    parser.add_argument("-eq", "--autoeq", type=str, default=False, help="if using auto equilibration")
    args = parser.parse_args()

    # check arguments
    print("checking the provided command line arguments...")
    main_dir, network_file, estimator, statsineff, eq = check_arguments(args)

    # check all the perturbations in the network
    perts, ligs = get_info_network(net_file=network_file)

    engines = validate.engines()

    file_ext = f"{estimator}_stats{statsineff}_eq{eq}"

    for engine in engines:

        spert_results_dict = {}
        spert_bound_dict = {}
        spert_free_dict = {}
        epert_results_dict = {}
        epert_bound_dict = {}
        epert_free_dict = {}
    
        for pert in perts:
            # find the extract folder
            # find correct path, use extracted if it exists
            if os.path.exists(f"{main_dir}/outputs_extracted/{engine}/{pert}"):
                path_to_dir = f"{main_dir}/outputs_extracted/{engine}/{pert}"
                final_results_folder = f"{main_dir}/outputs_extracted/results"
            elif os.path.exists(f"{main_dir}/outputs/{engine}_extracted/{pert}"):
                path_to_dir = f"{main_dir}/outputs/{engine}_extracted/{pert}"
                final_results_folder = f"{main_dir}/outputs/results"
            else:
                path_to_dir = f"{main_dir}/outputs/{engine}/{pert}"
                final_results_folder = f"{main_dir}/outputs/results"

            try:
                # calculate all the truncated data
                sresults_dict, sbound_dict, sfree_dict = truncation_check(path_to_dir, estimator, "start", statsineff, eq)
                eresults_dict, ebound_dict, efree_dict = truncation_check(path_to_dir, estimator, "end", statsineff, eq)

                spert_results_dict[pert] = sresults_dict
                spert_bound_dict[pert] = sbound_dict
                spert_free_dict[pert] = sfree_dict
                epert_results_dict[pert] = eresults_dict
                epert_bound_dict[pert] = ebound_dict
                epert_free_dict[pert] = efree_dict

                sdf = single_pert_dict_into_df(spert_results_dict[pert])
                edf = single_pert_dict_into_df(epert_results_dict[pert])
                # plot individually for perts
                print(f"plotting for {pert}, {engine} in {path_to_dir}/forward_reverse_{pert}.png...")
                plot_truncated(sdf, edf, file_path=f"{path_to_dir}/forward_reverse_{pert}_{file_ext}.png", plot_difference=False)

            except:
                print(f"could not calculate and plot for {pert} in {engine}.")

        print(f"plotting diff to final result for all perturbations in {engine}...")
        sdf = pert_dict_into_df(spert_results_dict, plot_error=False, plot_difference=True)
        edf = pert_dict_into_df(epert_results_dict, plot_error=False, plot_difference=True)
        plot_truncated(sdf, edf, f"{final_results_folder}/plt_truncated_{engine}_difference_forward_reverse_{file_ext}.png", plot_difference=True)
        print(f"saved images in {final_results_folder}.")

if __name__ == "__main__":
    main()