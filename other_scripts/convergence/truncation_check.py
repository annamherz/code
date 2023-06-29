#!/usr/bin/python3

import BioSimSpace as BSS
from BioSimSpace import Units as _Units
import os
import sys
from argparse import ArgumentParser

print("adding code to the pythonpath...")
code = "/home/anna/Documents/code/python"
if code not in sys.path:
    sys.path.insert(1, code)
import pipeline

from pipeline.prep import *
from pipeline.utils import *
from pipeline.analysis import *

import pandas as pd
import matplotlib.pyplot as plt


def check_arguments(args):
    # pass the checks to the other check functions
    if args.main_folder:
        main_folder = args.main_folder
    else:
        main_folder = str(
            input("what is the main folder where all the files should go?: ")
        ).strip()

    if args.network_file:
        network_file = args.network_file
    else:
        network_file = str(input("what is the path to the network file?: ").strip())

    if args.estimator:
        estimator = validate.estimator(args.estimator)
    else:
        estimator = validate.estimator(
            str(input("what is the analysis method (MBAR/TI) ? : ").strip())
        )

    if args.stats:
        statsineff = validate.boolean(args.stats)
    else:
        statsineff = False

    if args.autoeq:
        eq = validate.boolean(args.autoeq)
    else:
        eq = False

    if args.protocol_file:
        protocol_file = args.protocol_file
    else:
        protocol_file = None

    return main_folder, network_file, estimator, statsineff, eq, protocol_file


def main():
    # accept all options as arguments
    parser = ArgumentParser(description="run the convergence analysis")
    parser.add_argument(
        "-mf",
        "--main_folder",
        type=str,
        default=None,
        help="main folder path to create for all the runs",
    )
    parser.add_argument(
        "-n", "--network_file", type=str, default=None, help="path to network file"
    )
    parser.add_argument(
        "-e",
        "--estimator",
        type=str,
        default="MBAR",
        help="analysis method (MBAR or TI)",
    )
    parser.add_argument(
        "-s",
        "--stats",
        type=str,
        default=False,
        help="if using statistical inefficiency",
    )
    parser.add_argument(
        "-eq", "--autoeq", type=str, default=False, help="if using auto equilibration"
    )
    parser.add_argument(
        "-p", "--protocol_file", type=str, default=None, help="path to protocol file"
    )
    args = parser.parse_args()

    # check arguments
    print("checking the provided command line arguments...")
    main_dir, network_file, estimator, statsineff, eq, prot_file = check_arguments(args)

    # check all the perturbations in the network
    perts, ligs = get_info_network(net_file=network_file)

    engines = validate.engines()

    file_ext = f"{estimator}_stats{statsineff}_eq{eq}"

    for engine in engines:
        for pert in perts:
            # options
            analysis_options = analysis_protocol(
                {
                    "estimator": estimator,
                    "statistical inefficiency": statsineff,
                    "auto equilibration": eq,
                },
                auto_validate=True,
            )

            pert_name = None

            if prot_file:
                protocol = pipeline_protocol(
                    prot_file
                )  # instantiate the protocol as an object
                protocol.validate()  # validate all the input

                if protocol.name():
                    pert_name = f"{pert}_{protocol.name()}"
                    analysis_options.name(protocol.name())

            if not pert_name:
                pert_name = pert

            # find correct path, use extracted if it exists
            if os.path.exists(f"{main_dir}/outputs_extracted/{engine}/{pert_name}"):
                path_to_dir = f"{main_dir}/outputs_extracted/{engine}/{pert_name}"
                final_results_folder = f"{main_dir}/outputs_extracted/results"
            elif os.path.exists(f"{main_dir}/outputs/{engine}_extracted/{pert_name}"):
                path_to_dir = f"{main_dir}/outputs/{engine}_extracted/{pert_name}"
                final_results_folder = f"{main_dir}/outputs/results"
            else:
                path_to_dir = f"{main_dir}/outputs/{engine}/{pert_name}"
                final_results_folder = f"{main_dir}/outputs/results"

            if not os.path.exists(path_to_dir):
                raise OSError(f"{path_to_dir} does not exist.")

            print(f"analysing results for {path_to_dir}")

            # using the pipeline module for analysis
            analysed_pert = analyse(
                path_to_dir, pert=pert, analysis_prot=analysis_options
            )

            # plot the convergence
            analysed_pert.calculate_convergence()
            analysed_pert.plot_convergence()


if __name__ == "__main__":
    main()
