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

def truncation_check(path_to_dir, estimator):

    # analyse the work dir
    analysed_pert = analyse(path_to_dir)

    truncate_percentage = [5, 10, 15, 20, 25, 30, 40, 50, 60, 80, 100]
    
    results_dict = {}
    bound_dict = {}
    free_dict = {}

    for trunc_per in truncate_percentage:

        analysed_pert.set_options({ "estimator" : estimator,
                                    "truncate percentage": trunc_per,
                                    "truncate keep":"start"})

        analysed_pert.analyse_all_repeats()

        # save output in the pickle dir 
        write_analysis_file(analysed_pert, analysed_pert._pickle_dir, method=trunc_per)

        results_dict[trunc_per] = (analysed_pert.freenrg_val, analysed_pert.freenrg_err)
        bound_dict[trunc_per] = (analysed_pert.bound_val, analysed_pert.bound_err)
        free_dict[trunc_per] = (analysed_pert.free_val, analysed_pert.free_err)

    return results_dict, bound_dict, free_dict

def plot_truncated(pert_dict, file_path=None, plot_error=False, plot_difference=False, plot_diff="end"):

    df = pd.DataFrame.from_dict(pert_dict)
    perts = list(df.columns)
    df = df.reset_index().dropna()
    if file_path:
        df.to_csv(f"{file_path.replace('.png','.csv')}")

    include_key = False

    plt.rc('font', size=12)
    fig, ax = plt.subplots(figsize=(10,10))

    for pert in perts:
        col = next(ax._get_lines.prop_cycler)['color']
        x_vals = []
        y_vals = []
        lines = []
        for a,x in zip(df[pert], df['index']):
            if plot_error:
                if plot_difference:
                    if plot_diff == "end":
                        y_vals.append(abs(df.iloc[-1][pert][1] - a[1]))
                    elif plot_diff == "average":
                        y_vals.append(df.iloc[-1][pert][1] - a[1])
                else:
                    y_vals.append(a[1])
            else:
                if plot_difference:
                     y_vals.append(df.iloc[-1][pert][0] - a[0])
                else:
                    y_vals.append(a[0])
            x_vals.append(x)
        scatterplot = [plt.scatter(x_vals, y_vals, c=col)]   
        lines += plt.plot(0,0,label=pert)
        plt.errorbar(x_vals, y_vals, ecolor='none', color=col)

        if include_key:
            labels = [l.get_label() for l in lines]
            plt.legend(lines, labels, loc='upper left')
        
        plt.xlabel('Percentage of run used (from start)')
        plt.title(f"{file_path.split('/')[-1].split('.')[0].replace('_',' ')}")

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

    return main_folder, network_file, estimator

def main():

    # accept all options as arguments
    parser = ArgumentParser(description="run the ligprep")
    parser.add_argument("-mf", "--main_folder", type=str, default=None, help="main folder path to create for all the runs")
    parser.add_argument("-n", "--network_file", type=str, default=None, help="path to network file")
    parser.add_argument("-e", "--estimator", type=str, default="MBAR", help="analysis method (MBAR or TI)")
    args = parser.parse_args()

    # check arguments
    print("checking the provided command line arguments...")
    main_dir, network_file, estimator = check_arguments(args)

    # check all the perturbations in the network
    perts, ligs = get_info_network(net_file=network_file)

    engines = validate.engines()

    for engine in engines:

        pert_results_dict = {}
        pert_bound_dict = {}
        pert_free_dict = {}

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

            # calculate all the truncated data
            results_dict, bound_dict, free_dict = truncation_check(path_to_dir, estimator)

            pert_results_dict[pert] = results_dict
            pert_bound_dict[pert] = bound_dict
            pert_free_dict[pert] = free_dict


        # plot all the truncated data 
        print("plotting...")
        plot_truncated(pert_results_dict, f"{final_results_folder}/plt_truncated_{engine}.png")
        plot_truncated(pert_results_dict, f"{final_results_folder}/plt_truncated_error_{engine}.png", plot_error=True)
        plot_truncated(pert_bound_dict, f"{final_results_folder}/plt_truncated_bound_{engine}.png")
        plot_truncated(pert_bound_dict, f"{final_results_folder}/plt_truncated_bound_error_{engine}.png", plot_error=True)
        plot_truncated(pert_free_dict, f"{final_results_folder}/plt_truncated_free_{engine}.png")
        plot_truncated(pert_free_dict, f"{final_results_folder}/plt_truncated_free_error_{engine}.png", plot_error=True)
        # difference to final result
        plot_truncated(pert_results_dict, f"{final_results_folder}/plt_truncated_{engine}_difference.png", plot_difference=True)
        plot_truncated(pert_results_dict, f"{final_results_folder}/plt_truncated_error_{engine}_difference.png", plot_error=True, plot_difference=True)
        plot_truncated(pert_bound_dict, f"{final_results_folder}/plt_truncated_bound_{engine}_difference.png", plot_difference=True)
        plot_truncated(pert_bound_dict, f"{final_results_folder}/plt_truncated_bound_error_{engine}_difference.png", plot_error=True, plot_difference=True)
        plot_truncated(pert_free_dict, f"{final_results_folder}/plt_truncated_free_{engine}_difference.png", plot_difference=True)
        plot_truncated(pert_free_dict, f"{final_results_folder}/plt_truncated_free_error_{engine}_difference.png", plot_error=True, plot_difference=True)
        # difference to average
        print(f"saved images in {final_results_folder}.")

if __name__ == "__main__":
    main()