#!/usr/bin/python3

from argparse import ArgumentParser
import os
from pipeline.setup import *


def ask_things():
    lig_ff = str(
        input("what is the Ligand forcefield? (allowed: GAFF2, Parsely, Sage): ")
    ).strip()
    sampling_time = str(input("what is the sampling time? (in ns): ")).strip()
    engine = str(input("what is the engine? (allowed: SOMD, AMBER, GROMACS): ")).strip()
    hmr = str(input("should HMR be applied? (allowed: True, False): ")).strip()
    repeats = str(input("how many repeats? (recommended: 3): ")).strip()
    trajectories = str(
        input("save the trajectories? (allowed: 'None', '0,0.5,1', '0,1', 'All'): ")
    ).strip()

    # make these into a dictionary for the setup
    protocol_dict = {
        "ligand forcefield": lig_ff,
        "sampling": sampling_time,
        "hmr": hmr,
        "repeats": repeats,
        "trajectories": trajectories,
        "engine": engine,
    }

    return protocol_dict


def check_arguments(pl, args):
    # pass the checks to the other check functions
    if args.ligands_folder:
        pl.ligands_folder(args.ligands_folder)
    else:
        pl.ligands_folder(str(input("what is the ligands folder?: ")).strip())

    if args.protein_path:
        pl.protein_path(args.protein_path)
    else:
        pl.protein_path(str(input("what is the parameterised protein path?: ")).strip())

    if args.main_folder:
        pl.main_folder(args.main_folder)
    else:
        pl.main_folder(
            str(
                input("what is the main folder where all the files should go?: ")
            ).strip()
        )

    return pl


def main():
    print("setup the pipeline! first choose options")

    # accept all options as arguments
    parser = ArgumentParser(description="set up simulations")
    parser.add_argument(
        "-lf",
        "--ligands_folder",
        type=str,
        default=None,
        help="folder path to the ligand files",
    )
    parser.add_argument(
        "-pf",
        "--protein_path",
        type=str,
        default=None,
        help="path to parameterised protein *.prm7 and *.rst7",
    )
    parser.add_argument(
        "-mf",
        "--main_folder",
        type=str,
        default=None,
        help="main folder path to create for all the runs",
    )
    parser.add_argument(
        "-m",
        "--method",
        type=str,
        default="benchmark",
        help="descriptor of what this run is for.",
    )
    args = parser.parse_args()

    # intialise setup class
    pl = initialise_pipeline()

    # check arguments
    print("checking the provided command line arguments...")
    pl = check_arguments(pl, args)

    print("please decide some basic protocol settings:")
    protocol_dict = ask_things()
    # fill in rest with default and save the files
    pl.setup_protocols(protocol_dict)

    # make the files needed
    pl.setup_ligands()
    pl.setup_network()

    # make run_all_slurm to the main folder that was made.
    pl.write_run_all()

    print(
        f"made all files in {pl.main_folder()}. Run run_all_slurm.sh to run all pipeline components."
    )
    print(
        "carefully check the run all script first to make sure all paths and executables are correct and can be sourced correctly."
    )
    print("modify any of the input files as required.")


if __name__ == "__main__":
    main()
