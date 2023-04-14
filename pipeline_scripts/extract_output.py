#!/usr/bin/python3

# import BioSimSpace as BSS
from argparse import ArgumentParser

from pipeline.utils import *
from pipeline.prep import *

def extract_output(folder, prot_file):
    # read in protocol
    protocol = pipeline_protocol(prot_file, auto_validate=True)
    if protocol.trajectories() == "None":
        traj_lambdas = []
    if protocol.trajectories() == "0,0.5,1":
        traj_lambdas = ["0.0000","0.5000","1.0000"]
    if protocol.trajectories() == "0,1":
        traj_lambdas = ["0.0000","1.0000"]
    if protocol.trajectories() == "All":
        # TODO get lambda windows from network file
        traj_lambdas = ["0.0000","0.5000","1.0000"]

    # simfile header
    if "SOMD" in folder:
        try:
            add_header_simfile(folder)
        except:
            print(f"could not add the header to simfile in {folder}")

    # extract to output folder
    # initialise
    extraction = extract(folder)

    # get the output from the folder to new folder
    extraction.extract_output()
    # get trajectory, will get rmsd by default
    extraction.extract_frames(traj_lambdas=traj_lambdas, overwrite=True)

def check_arguments(args):

    # pass the checks to the other check functions
    if args.folder:
        folder = args.folder
    else:
        folder = str(input("what is the folder to extract?: ")).strip()

    if args.protocol_file:
        protocol_file = args.protocol_file
    else:
        protocol_file = str(input("what is the path to the protocol file?: ").strip())

    return folder, protocol_file

def main():

    # accept all options as arguments
    parser = ArgumentParser(description="extract the output")
    parser.add_argument("-f", "--folder", type=str, default=None, help="folder to extract")
    parser.add_argument("-p", "--protocol_file", type=str, default=None, help="path to protocol file")
    args = parser.parse_args()

    # check arguments
    print("checking the provided command line arguments...")
    folder, protocol_file = check_arguments(args)

    extract_output(folder, protocol_file)

if __name__ == "__main__":
    main()
