# import BioSimSpace as BSS
import sys

from pipeline.utils import *

# need the transdir, sys argv 1
folder = sys.argv[1]
# folder format is $MAINDIRECTORY/outputs/$2/$1

prot_file = os.environ["prot_file"] # protocol file

# read in protocol
protocol = pipeline_protocol(prot_file, auto_validate=True)
if protocol.trajectories == "None":
    traj_lambdas = []
if protocol.trajectories == "0,0.5,1":
    traj_lambdas = ["0.0000","0.5000","1.0000"]
if protocol.trajectories == "0,1":
    traj_lambdas = ["0.0000","1.0000"]
if protocol.trajectories == "All":
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
