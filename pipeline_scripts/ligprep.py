# ligand prep

import BioSimSpace as BSS
from BioSimSpace import _Exceptions
from BioSimSpace.Units.Length import angstrom as _angstrom
import sys
import os

from pipeline.utils import *
from pipeline.prep import *

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
pmemd_path = os.environ["amber"] + "/bin/pmemd.cuda"
main_dir = os.environ["MAINDIRECTORY"]
protein_file = os.environ["protein_file"]
ligands_folder = os.environ["ligands_folder"]
prot_file = os.environ["prot_file"] # protocol.dat which was set in overall bash script
lig_name = sys.argv[1]
print(f'prep for {lig_name}...')

# Exit if prep files are already available for ligand and protein
if os.path.exists(f"{main_dir}/prep/{lig_name}_sys_equil_solv.prm7"):
    if os.path.exists(f"{main_dir}/prep/{lig_name}_lig_equil_solv.prm7"):
        print(f"Prep files already generated for {lig_name}. Done.")
        sys.exit(0)

# make other folders that need to be made for the prepped files
folders = [f"{main_dir}/prep", f"{main_dir}/prep/pre_run"]
for folder in folders:
    if not os.path.exists(folder):
        os.mkdir(folder)
    else:
        pass

# parse protocol file
protocol = check_protocol(prot_file) # instantiate the protocol as an object
protocol.validate() # validate all the input into needed format.

# load, solvate and run the systems.
# load the protein, this was paramaterised during the setup stage.
prot_wat = BSS.IO.readMolecules([f"{protein_file}.rst7", f"{protein_file}.prm7"])

# load ligand, these should already be in the correct position
try: # sdf first
    ligand = BSS.IO.readMolecules(f"{ligands_folder}/{lig_name}.sdf")[0]
except: # mol2 if sdf is not available
    ligand = BSS.IO.readMolecules(f"{ligands_folder}/{lig_name}.mol2")[0]

# paramaterise the ligand
print(f"Parameterising {lig_name}...")
lig_p = lig_paramaterise(ligand, protocol.ligand_forcefield).getMolecule()

# Combine protein, ligand and crystallographic waters.
system = lig_p + prot_wat

# Solvate and run each the bound and the free system.
legs_mols, legs = [lig_p, system], ["lig", "sys"]

# zip together the molecules in that leg with the name for that leg
for leg, leg_mol in zip(legs, legs_mols):

    # solvate
    print(f"Solvating {leg} for {lig_name}...")
    leg_mol_solvated = minimum_solvation(leg_mol,
                                protocol.solvent,
                                protocol.box_type,
                                protocol.box_edges,
                                protocol.box_edges_unit)

    # saving pre runs
    print(f"Saving solvated for {leg} and {lig_name}")
    BSS.IO.saveMolecules(f"{main_dir}/prep/pre_run/{lig_name}_{leg}_s",
                         leg_mol_solvated, ["PRM7", "RST7", "PDB"])

    # minimise and eqiulibrate these
    print(f"minimising and equilibrating for {leg} and {lig_name}")
    leg_equil_final = minimise_equilibrate_leg(leg_mol_solvated, leg, engine, pmemd_path)

    # finally, save last snapshot
    print(f"Saving solvated/equilibrated for {leg} and {lig_name}")
    BSS.IO.saveMolecules(
        f"{main_dir}/prep/{lig_name}_{leg}_equil_solv", leg_equil_final, ["PRM7", "RST7", "PDB"])
    print(f"Done for {lig_name}.")
