#!/usr/bin/python3

# FEP prep
# get all inputs from the overall pipeline
# create the free and bound merged systems
# validate the protocol
# create fepprep objects with the systems and protocol
# generate the folders

# Import libraries
import BioSimSpace as BSS
import os
import sys
import csv
from distutils.dir_util import remove_tree

from pipeline.prep import *
from pipeline.utils import *

# Values from the corresponding bash arrays.
# The ligands for which the transformation is to be carried out.
trans = sys.argv[1]
lig_1 = trans.split('~')[0]
lig_2 = trans.split('~')[1]
engine_query = str(sys.argv[2]).upper()
num_lambda_query = int(sys.argv[3])

# files that were set in the run_all script
pmemd_path = os.environ["amber"] + "/bin/pmemd.cuda"
main_dir = os.environ["MAINDIRECTORY"]
net_file = os.environ["net_file"] # network file
prot_file = os.environ["prot_file"] # protocol file
prep_dir = f"{main_dir}/prep"  # define lig prep location
workdir = f"{main_dir}/outputs/{engine_query}/{lig_1}~{lig_2}" # pert dir

# check if the trans, eng and window are in the network file
found = False
with open(net_file, "r") as lambdas_file:
    reader = csv.reader(lambdas_file, delimiter=" ")
    for row in reader:
        if (row[0] == lig_1 and row[1] == lig_2) or (row[1] == lig_1 and row[0] == lig_2):
            if int(row[2]) == num_lambda_query:
                if str(row[-1]).upper() == engine_query:
                    found = True

if not found:
    raise NameError(
        f"The perturbation {trans} (or the reverse) with {num_lambda_query} windows using {engine_query} was not found in {net_file}.")

# parse protocol file
print("reading in the protocol file...")
protocol = pipeline_protocol(prot_file) # instantiate the protocol as an object
print("validating the protocol file...")
protocol.validate() # validate all the input
# print("rewriting the protocol file...")
# protocol.rewrite_protocol() # rewrite protocol file
# print("the protocol is now:")
# protocol.print_protocol()
# add the number of lambdas and engine to the protocol
protocol.num_lambda(num_lambda_query)
protocol.engine(engine_query)

# instantiate each system as a fepprep class with the protocol
fepprep_obj = fepprep(protocol=protocol)

for name, leg in zip(["lig", "sys"], ["free", "bound"]):

    # Load equilibrated inputs for both ligands
    system_1 = BSS.IO.readMolecules(
        [f"{prep_dir}/{lig_1}_{name}_equil_solv.rst7", f"{prep_dir}/{lig_1}_{name}_equil_solv.prm7"])
    system_2 = BSS.IO.readMolecules(
        [f"{prep_dir}/{lig_2}_{name}_equil_solv.rst7", f"{prep_dir}/{lig_2}_{name}_equil_solv.prm7"])

    fepprep_obj.add_system(system_1, free_bound=leg, start_end="start")
    fepprep_obj.add_system(system_2, free_bound=leg, start_end="end")

# remove any existing files in the workdir
try:    
    remove_tree(workdir)
except:
    pass
# generate folder based on fepprep protocol (both or start)
fepprep_obj.generate_folders(workdir)
