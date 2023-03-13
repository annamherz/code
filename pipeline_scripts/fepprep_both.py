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
import shutil

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
pmemd_path = os.environ["AMBERHOME"] + "/bin/pmemd.cuda"
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

# create the system for each the free and the bound leg.
system_free = None
system_bound = None

for name, leg in zip(["lig", "sys"], ["free", "bound"]):
    # Load equilibrated inputs for both ligands
    system_1 = BSS.IO.readMolecules(
        [f"{prep_dir}/{lig_1}_{name}_equil_solv.rst7", f"{prep_dir}/{lig_1}_{name}_equil_solv.prm7"])
    system_2 = BSS.IO.readMolecules(
        [f"{prep_dir}/{lig_2}_{name}_equil_solv.rst7", f"{prep_dir}/{lig_2}_{name}_equil_solv.prm7"])

    print(f"Preparing the {leg} leg...")
    if leg == "free":
        system_free = merge.merge_system(system_1, system_2)
    if leg == "bound":
        system_bound = merge.merge_system(system_1, system_2)

# instantiate each system as a fepprep class with the protocol
fepprep_obj = fepprep(system_free, system_bound, protocol)
# then generate all folders starting from 0.5
fepprep_obj.generate_folders(workdir)

workdir2 = f"{main_dir}/outputs/{engine_query}/{lig_1}~{lig_2}_temp"
# create the system for each the free and the bound leg.
system_free = None
system_bound = None

for name, leg in zip(["lig", "sys"], ["free", "bound"]):
    # Load equilibrated inputs for both ligands
    system_1 = BSS.IO.readMolecules(
        [f"{prep_dir}/{lig_1}_{name}_equil_solv.rst7", f"{prep_dir}/{lig_1}_{name}_equil_solv.prm7"])
    system_2 = BSS.IO.readMolecules(
        [f"{prep_dir}/{lig_2}_{name}_equil_solv.rst7", f"{prep_dir}/{lig_2}_{name}_equil_solv.prm7"])

    print(f"Preparing the {leg} leg...")
    # swap order of mapping systems so maps to system 2
    if leg == "free":
        system_free = merge.merge_system(system_2, system_1)
    if leg == "bound":
        system_bound = merge.merge_system(system_2, system_1)

# instantiate each system as a fepprep class with the protocol
fepprep_obj = fepprep(system_free, system_bound, protocol)
# then generate all folders starting from 0.5
fepprep_obj.generate_folders(workdir2)

end_lambdas = ["0.6000", "0.7000", "0.8000", "0.9000", "1.0000"]
start_lambdas = ["0.0000", "0.1000", "0.2000", "0.3000", "0.4000", "0.5000"]

for lam in end_lambdas:
    os.rmdir(f"{workdir}/*{lam}")

for lam in end_lambdas:
    print("moving end lambdas to main run folder")
    shutil.move(f"{workdir2}/lambda_{lam}", f"{workdir}/lambda_{lam}")
    shutil.move(f"{workdir2}/min/lambda_{lam}", f"{workdir}/min/lambda_{lam}")
    shutil.move(f"{workdir2}/eq/lambda_{lam}", f"{workdir}/eq/lambda_{lam}")
    shutil.move(f"{workdir2}/heat/lambda_{lam}", f"{workdir}/heat/lambda_{lam}")

os.rmdir(f"{workdir2}")