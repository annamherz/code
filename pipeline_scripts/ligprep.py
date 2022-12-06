# ligand prep

import BioSimSpace as BSS
from BioSimSpace import _Exceptions
from BioSimSpace.Units.Length import angstrom as _angstrom
import sys
import os

BSS.setVerbose(True)

try:
    import pipeline
except:
    print("adding code to the pythonpath...")
    code = '/home/anna/Documents/code/python'
    if code not in sys.path:
        sys.path.insert(1, code)
    import pipeline

from pipeline.prep import *

# decide engine used for equilibration runs ("AMBER" or "GROMACS")
engine = "AMBER"

# set environs
# pmemd_path = os.environ["AMBERHOME"] + "/bin/pmemd.cuda"
pmemd_path = os.environ["amber"] + "/bin/pmemd.cuda"
main_dir = os.environ["MAINDIRECTORY"]
# main_dir = "/home/anna/Documents/benchmark/test_tyk2_benchmark_sage"
protein_file = os.environ["protein_file"]
ligands_folder = os.environ["ligands_folder"]
# ref the protocol.dat which was set in overall bash script
prot_file = os.environ["prot_file"]
lig_name = sys.argv[1]
print(f'prep for {lig_name}...')

# Exit if prep files are already available for ligand and protein
if os.path.exists(f"{main_dir}/prep/{lig_name}_sys_equil_solv.prm7"):
    if os.path.exists(f"{main_dir}/prep/{lig_name}_lig_equil_solv.prm7"):
        print(f"Prep files already generated for {lig_name}. Done.")
        sys.exit(0)

# make other folders that need to be made
# for the prepped files
if not os.path.exists(f"{main_dir}/prep"):
    os.mkdir(f"{main_dir}/prep")
else:
    pass

# for the system pre prep
if not os.path.exists(f"{main_dir}/prep/pre_run"):
    os.mkdir(f"{main_dir}/prep/pre_run")
else:
    pass

################
# Defining all the parameters needed

# Read the protocol.dat to get all parameter infos.
search_dict = {"ligand forcefield":None, "protein forcefield":None,
                "solvent":None, "box edges":None, "box type": None,
                "protocol":None, "sampling":None, "HMR":None,
                "repeats":None, "keep trajectories":None}

# get value regardless of the order of the protocol.dat
for search in search_dict.keys():
    with open(f"{prot_file}", "r") as file:
        for line in file:
            if search.casefold() in line.casefold():
                search_dict[search] = line.strip()
                break

# forcefields
ligff_query = (search_dict["ligand forcefield"].rstrip().replace(" ", "").split("=")[-1]).lower() # ligand ff
protff_query = search_dict["protein forcefield"].rstrip().replace(" ", "").split("=")[-1] # protein ff
solvent_query = search_dict["solvent"].rstrip().replace(" ", "").split("=")[-1] # solvent ff

# box size and type
boxsize_query = search_dict["box edges"].rstrip().replace(" ", "").split("=")[-1]
box_axis_length = boxsize_query.split("*")[0]
box_axis_unit_query = boxsize_query.split("*")[1]
boxtype_query = (search_dict["box type"].rstrip().replace(" ", "").split("=")[-1]).lower()

# load, solvate and run the systems.
# load the protein, this was paramaterised during the setup stage.
prot_wat = BSS.IO.readMolecules(
    [f"{protein_file}.rst7", f"{protein_file}.prm7"])

# load ligand
# These should already be in the correct position
# try sdf first, if not available try mol2
try:
    ligand = BSS.IO.readMolecules(f"{ligands_folder}/{lig_name}.sdf")[0]
except:
    ligand = BSS.IO.readMolecules(f"{ligands_folder}/{lig_name}.mol2")[0]

# paramaterise the ligand
print(f"Parameterising {lig_name}...")
lig_p = lig_paramaterise(ligand, ligff_query).getMolecule()

# Combine protein, ligand and crystallographic waters.
system = lig_p + prot_wat

# Solvate and run each the bound and the free system.
legs_mols = [lig_p, system]
legs = ["lig", "sys"]

# zip together the molecules in that leg with the name for that leg
for leg, leg_mol in zip(legs, legs_mols):

    # solvate
    print(f"Solvating {leg} for {lig_name}...")
    leg_mol_solvated = min_solv(leg_mol, solvent_query, boxtype_query, box_axis_length, box_axis_unit_query)

    # saving pre runs
    print(f"Saving solvated for {leg} and {lig_name}")
    BSS.IO.saveMolecules(f"{main_dir}/prep/pre_run/{lig_name}_{leg}_s",
                         leg_mol_solvated, ["PRM7", "RST7", "PDB"])

    # minimise and eqiulibrate these
    print(f"minimising and equilibrating for {leg} and {lig_name}")
    leg_equil_final = lig_paramateriserep(leg_mol_solvated, leg)

    # finally, save last snapshot
    print(f"Saving solvated/equilibrated for {leg} and {lig_name}")
    BSS.IO.saveMolecules(
        f"{main_dir}/prep/{lig_name}_{leg}_equil_solv", leg_equil_final, ["PRM7", "RST7", "PDB"])
    print(f"Done for {lig_name}.")
