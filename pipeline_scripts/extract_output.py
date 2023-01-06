# import BioSimSpace as BSS
import sys

from pipeline.utils import *

# need the transdir, sys argv 1
folder = sys.argv[1]
# folder format is $MAINDIRECTORY/outputs/$2/$1

# simfile header
if "SOMD" in folder:
    try:
        add_header_simfile(folder)
    except:
        print(f"could not add the header to simfile in {folder}")

# extract to output folder
extract_output_single(folder)

# get trajectory