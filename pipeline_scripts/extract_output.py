# import BioSimSpace as BSS
import sys

try:
    import pipeline
except:
    print("adding code to the pythonpath...")
    code = '/home/anna/Documents/code/python'
    if code not in sys.path:
        sys.path.insert(1, code)
    import pipeline
    
from pipeline.utils import *

# need the transdir, sys argv 1
folder = sys.argv[1]
# folder format is $MAINDIRECTORY/outputs/$2/$1

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
