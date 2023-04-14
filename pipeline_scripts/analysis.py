#!/usr/bin/python3

# analysis script for a single pert
# anna

import BioSimSpace as BSS
import sys
import os as _os
from argparse import ArgumentParser

BSS.setVerbose = True

from pipeline.analysis import *
from pipeline.utils import write_analysis_file
from pipeline.prep import *


def analysis(pert, engine, ana_file, main_dir):

    # options
    analysis_options = analysis_protocol(ana_file, auto_validate=True)
    analysis_options.rewrite_protocol()

    # TODO final results folder has a method for the analysis option used?

    # find correct path, use extracted if it exists
    if _os.path.exists(f"{main_dir}/outputs_extracted/{engine}/{pert}"):
        path_to_dir = f"{main_dir}/outputs_extracted/{engine}/{pert}"
        final_results_folder = f"{main_dir}/outputs_extracted/results"
    elif _os.path.exists(f"{main_dir}/outputs/{engine}_extracted/{pert}"):
        path_to_dir = f"{main_dir}/outputs/{engine}_extracted/{pert}"
        final_results_folder = f"{main_dir}/outputs/results"
    else:
        path_to_dir = f"{main_dir}/outputs/{engine}/{pert}"
        final_results_folder = f"{main_dir}/outputs/results"

    if not _os.path.exists(path_to_dir):
        raise OSError(f"{path_to_dir} does not exist.")

    print(f'analysing results for {path_to_dir}')
    print(f"using {analysis_options} for analysis")

    # using the pipeline module for analysis
    analysed_pert = analyse(path_to_dir)
    analysed_pert.set_options(analysis_options)
    avg, error, repeats_tuple_list = analysed_pert.analyse_all_repeats()
    analysed_pert.plot_graphs()

    # write the final result
    write_analysis_file(analysed_pert, final_results_folder)

def check_arguments(args):

    # pass the checks to the other check functions
    if args.perturbation:
        perturbation = args.perturbation
    else:
        perturbation = str(input("what is the perturbation?: ")).strip()

    if args.engine:
        engine = args.engine
    else:
        engine = str(input("what is the engine?: ").strip())

    if args.main_folder:
        main_folder = args.main_folder
    else:
        main_folder = str(input("what is the main folder of the run?: ")).strip()

    if args.analysis_file:
        analysis_file = args.analysis_file
    else:
        analysis_file = str(input("what is the path to the analysis protocol file?: ").strip())

    return perturbation, engine, analysis_file, main_folder

def main():

    # accept all options as arguments
    parser = ArgumentParser(description="run the fepprep")
    parser.add_argument("-pert", "--perturbation", type=str, default=None, help="name of perturbation")
    parser.add_argument("-eng", "--engine", type=str, default=None, help="engine of the run")
    parser.add_argument("-mf", "--main_folder", type=str, default=None, help="main folder path for all the runs")
    parser.add_argument("-a", "--analysis_file", type=str, default=None, help="path to analysis protocol file")
    args = parser.parse_args()

    # check arguments
    print("checking the provided command line arguments...")
    pert, engine, ana_file, main_dir = check_arguments(args)

    print(f"analysis for {pert, engine, main_dir}")

    analysis(pert, engine, ana_file, main_dir)

if __name__ == "__main__":
    main()
