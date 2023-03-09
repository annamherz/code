#!/usr/bin/python3

from argparse import ArgumentParser
import os
from pipeline.setup import *
from pipeline.utils import *

# ligand ff
# samplign time
# repeats
# hmr 
# trajectories
# FEP engine

def ask_things():

    is_change = True
    change_dict = {}

    while is_change:
        change = str(input("what is the parameter to change? : ")).strip()
        change_value = str(input("what is the change? : ")).strip()
        change_dict[change] = change_value
        is_change = validate.boolean(input("change something else? (True / False) : "))

    return change_dict

def check_arguments(pl, args):

    # current main folder

    # read current network
    # choose if same network or not
    # if not, choose which perturbations are keeping instead
    # or if there are any others
    # take the changed protocol - read the current, apply changes
    # write new file with the method name for this

    # save the pipeline object as a pickle as well that can be opened?

    if args.main_folder:
        pl.main_folder(args.main_folder)
    else:
        pl.main_folder(str(input("what is the main folder where all the files should go?: ")).strip())

    return pl

def main():

    print("setup the pipeline! first choose options")

    # accept all options as arguments
    parser = ArgumentParser(description="set up simulations")
    parser.add_argument("-lf", "--ligands_folder", type=str, default=None, help="folder path to the ligand files")
    parser.add_argument("-pf", "--protein_path", type=str, default=None, help="path to parameterised protein *.prm7 and *.rst7")
    parser.add_argument("-mf", "--main_folder", type=str, default=None, help="main folder path to create for all the runs")
    parser.add_argument("-m", "--method", type=str, default="benchmark", help="descriptor of what this run is for.")
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

    print(f"made all files in {pl.main_folder()}. Run run_all_slurm.sh to run all pipeline components.")
    print("carefully check the run all script first to make sure all paths are correct.")
    print("modify any of the input files as required.")

# TODO add setting of environment variables eg amber etc

if __name__ == "__main__":
    main()
