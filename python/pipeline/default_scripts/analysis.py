#!/usr/bin/python3

import BioSimSpace as BSS
import sys
import os as _os
from argparse import ArgumentParser
import logging

BSS.setVerbose = True

from pipeline.analysis import *
from pipeline.utils import write_analysis_file
from pipeline.prep import *


def analysis(pert, engine, ana_file, main_dir, prot_file=None):
    # options
    analysis_options = analysis_protocol(ana_file, auto_validate=True)
    pert_name = None

    if prot_file:
        protocol = pipeline_protocol(prot_file)  # instantiate the protocol as an object
        protocol.validate()  # validate all the input

        if protocol.name():
            pert_name = f"{pert}_{protocol.name()}"
            analysis_options.name(protocol.name())

    if not pert_name:
        pert_name = pert

    # find correct path, use extracted if it exists
    if _os.path.exists(f"{main_dir}/outputs_extracted/{engine}/{pert_name}"):
        path_to_dir = f"{main_dir}/outputs_extracted/{engine}/{pert_name}"
        final_results_folder = f"{main_dir}/outputs_extracted/results"
    else:
        path_to_dir = f"{main_dir}/outputs/{engine}/{pert_name}"
        final_results_folder = f"{main_dir}/outputs/results"

    if not _os.path.exists(path_to_dir):
        raise OSError(f"{path_to_dir} does not exist.")

    print(f"analysing results for {path_to_dir}")
    print(f"using {analysis_options.print_protocol()} for analysis")

    # using the pipeline module for analysis
    analysed_pert = analyse(path_to_dir, pert=pert, analysis_prot=analysis_options)
    avg, error, repeats_tuple_list = analysed_pert.analyse_all_repeats()
    analysed_pert.check_convergence()
    analysed_pert.plot_graphs()

    # write the final result
    write_analysis_file(analysed_pert, final_results_folder)

    # plot the convergence
    # analysed_pert.calculate_convergence()
    # analysed_pert.plot_convergence()
    # analysed_pert.plot_across_lambda()

    # write for edgembar
    # analysed_pert.format_for_edgembar()


def analysis_work_dir(work_dir, pert, engine, ana_file):
    if ana_file:
        analysis_options = analysis_protocol(ana_file, auto_validate=True)
    else:
        analysis_options = None

    analysed_pert = analyse(work_dir, pert, engine, analysis_options)
    avg, error, repeats_tuple_list = analysed_pert.analyse_all_repeats()
    analysed_pert.plot_graphs()

    print(avg, error, repeats_tuple_list)


def check_arguments(args):
    # pass the checks to the other check functions

    if args.work_dir:
        work_dir = args.work_dir
        perturbation = None
        engine = None
        main_folder = None
        analysis_file = None
        prot_file = None
    else:
        work_dir = None

        if args.main_folder:
            main_folder = args.main_folder
        else:
            main_folder = str(input("what is the main folder of the run?: ")).strip()

        if args.protocol_file:
            prot_file = args.protocol_file
        else:
            prot_file = None

    if args.perturbation:
        perturbation = args.perturbation
    else:
        perturbation = str(input("what is the perturbation?: ")).strip()

    if args.engine:
        engine = args.engine
    else:
        engine = str(input("what is the engine?: ").strip())

    if args.analysis_file:
        analysis_file = args.analysis_file
    else:
        if work_dir:
            analysis_file = None
        else:
            analysis_file = str(
                input("what is the path to the analysis protocol file?: ").strip()
            )

    return perturbation, engine, analysis_file, main_folder, prot_file, work_dir


def main():
    # accept all options as arguments
    parser = ArgumentParser(description="run the fepprep")
    parser.add_argument(
        "-pert", "--perturbation", type=str, default=None, help="name of perturbation"
    )
    parser.add_argument(
        "-eng", "--engine", type=str, default=None, help="engine of the run"
    )
    parser.add_argument(
        "-mf",
        "--main_folder",
        type=str,
        default=None,
        help="main folder path for all the runs",
    )
    parser.add_argument(
        "-a",
        "--analysis_file",
        type=str,
        default=None,
        help="path to analysis protocol file",
    )
    parser.add_argument(
        "-p", "--protocol_file", type=str, default=None, help="path to protocol file"
    )
    parser.add_argument(
        "-wd",
        "--work_dir",
        type=str,
        default=None,
        help="work dir of run, will ignore mf and protocol and ana args.",
    )
    args = parser.parse_args()

    # check arguments
    print("checking the provided command line arguments...")
    pert, engine, ana_file, main_dir, prot_file, work_dir = check_arguments(args)

    if work_dir:
        analysis_work_dir(work_dir, pert, engine, ana_file)
    else:
        print(f"analysis for {pert, engine, main_dir}")
        analysis(pert, engine, ana_file, main_dir, prot_file)


if __name__ == "__main__":
    main()
