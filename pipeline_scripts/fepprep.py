# FEP prep

# Import libraries
import BioSimSpace as BSS
from BioSimSpace import _Exceptions
import os
import sys
import csv
from distutils.dir_util import copy_tree

try:
    import pipeline
except:
    print("adding code to the pythonpath...")
    code = '/home/anna/Documents/code/python'
    if code not in sys.path:
        sys.path.insert(1, code)
    import pipeline

from pipeline.prep import *
from pipeline.utils import *

# Values from the corresponding bash arrays.
# The ligands for which the transformation is to be carried out.
trans = sys.argv[1]
lig_1 = trans.split('~')[0]
lig_2 = trans.split('~')[1]
engine_query = str(sys.argv[2]).upper()
num_lambda = int(sys.argv[3])

main_dir = os.environ["MAINDIRECTORY"]
net_file = os.environ["net_file"] # network file
prot_file = os.environ["prot_file"] # protocol file
prep_dir = f"{main_dir}/prep"  # define lig prep location
# define trans dir
workdir = f"{main_dir}/outputs/{engine_query}/{lig_1}~{lig_2}"

# check if the trans, eng and window are in the network file
found = False
with open(net_file, "r") as lambdas_file:
    reader = csv.reader(lambdas_file, delimiter=" ")
    for row in reader:
        if (row[0] == lig_1 and row[1] == lig_2) or (row[1] == lig_1 and row[0] == lig_2):
            if int(row[2]) == num_lambda:
                if str(row[-1]).upper() == engine_query:
                    found = True

if not found:
    raise NameError(
        f"The perturbation {trans} (or the reverse) with {num_lambda} windows using {engine_query} was not found in {net_file}.")

# parse protocol file
protocol = pipeline_protocol(prot_file) # instantiate the protocol as an object
protocol.validate() # validate all the input
protocol.num_lambda = validate.num_lambda(num_lambda)
protocol.engine = validate.engine(engine_query)

# Load equilibrated free inputs for both ligands.
# create the system for each the free and the bound leg.
system_free = None
system_bound = None

for name, leg in zip(["lig", "sys"], ["free", "bound"]):
    # Load equilibrated free inputs
    system_1 = BSS.IO.readMolecules(
        [f"{prep_dir}/{lig_1}_{name}_equil_solv.rst7", f"{prep_dir}/{lig_1}_{name}_equil_solv.prm7"])
    system_2 = BSS.IO.readMolecules(
        [f"{prep_dir}/{lig_2}_{name}_equil_solv.rst7", f"{prep_dir}/{lig_2}_{name}_equil_solv.prm7"])

    print(f"Preparing the {leg} leg...")

    if leg == "free":
        # Extract ligands.
        ligand_1 = system_1.getMolecule(0)
        ligand_2 = system_2.getMolecule(0)

    elif leg == "bound":
        # Extract ligands and protein. Do this based on nAtoms and nResidues, as sometimes
        # the order of molecules is switched, so we can't use index alone.
        ligand_1 = None
        protein = None
        n_residues = [mol.nResidues() for mol in system_1]
        n_atoms = [mol.nAtoms() for mol in system_1]
        for i, (n_resi, n_at) in enumerate(zip(n_residues[:20], n_atoms[:20])):
            if n_resi == 1 and n_at > 5:
                ligand_1 = system_1.getMolecule(i)
            elif n_resi > 1:
                protein = system_1.getMolecule(i)
            else:
                pass

        # loop over molecules in system to extract the ligand
        ligand_2 = None
        n_residues = [mol.nResidues() for mol in system_2]
        n_atoms = [mol.nAtoms() for mol in system_2]
        for i, (n_resi, n_at) in enumerate(zip(n_residues, n_atoms)):
            # grab the system's ligand and the protein. ignore the waters.
            if n_resi == 1 and n_at > 5:
                ligand_2 = system_2.getMolecule(i)
            else:
                pass

        if protein:
            pass
        else:
            raise _Exceptions.AlignmentError(
                "Could not extract protein from input systems. Check that your ligands/proteins are properly prepared!")

    if ligand_1 and ligand_2:
        pass
    else:
        raise _Exceptions.AlignmentError(
            "Could not extract ligands from input systems. Check that your ligands/proteins are properly prepared!")

    # merge the ligands based on the engine.
    print("mapping, aligning and merging the ligands...")
    merged_trans = merge_ligands(ligand_1, ligand_2, engine_query)

    # Get equilibrated waters and waterbox information for both bound and free. Get all information from lambda==0
    # Following is work around because setBox() doesn't validate correctly boxes with lengths and angles
    system_1.removeMolecules(ligand_1)
    system_final = merged_trans + system_1

    if leg == "free":
        system_free = system_final
    if leg == "bound":
        system_bound = system_final

# repartition the hydrogen masses
if protocol.hmr == True:
    print("repartitioning hydrogen masses for 4fs timestep...")
    if engine_query == "AMBER":
        system_free.repartitionHydrogenMass(factor=3)
        system_bound.repartitionHydrogenMass(factor=3)
    elif engine_query == "GROMACS":
        system_free.repartitionHydrogenMass(factor=4)
        system_bound.repartitionHydrogenMass(factor=4)
    elif engine_query == "SOMD":
        pass
elif protocol.hmr == False:
    pass

# define the free energy protocol with all this information.
if engine_query == 'AMBER' or engine_query == 'GROMACS':
    min_protocol = BSS.Protocol.FreeEnergyMinimisation(
        num_lam=num_lambda, steps=min_steps)
    heat_protocol = BSS.Protocol.FreeEnergyEquilibration(timestep=protocol.timestep*protocol.timestep_unit,
                                                         num_lam=protocol.num_lambda,
                                                         runtime=protocol.eq_runtime*protocol.eq_runtime_unit,
                                                         pressure=None,
                                                         temperature_start=protocol.start_temperature*protocol.temperature_unit,
                                                         temperature_end=protocol.end_temperature*protocol.temperature_unit
                                                         )
    eq_protocol = BSS.Protocol.FreeEnergyEquilibration(timestep=protocol.timestep*protocol.timestep_unit,
                                                       num_lam=protocol.num_lambda,
                                                       runtime=protocol.eq_runtime*protocol.eq_runtime_unit,
                                                       temperature=protocol.temperature*protocol.temperature_unit,
                                                       pressure=protocol.pressure*protocol.pressure_unit,
                                                       restart=True
                                                       )
    freenrg_protocol = BSS.Protocol.FreeEnergy(timestep=protocol.timestep*protocol.timestep_unit,
                                                num_lam=protocol.num_lambda,
                                                runtime=protocol.sampling*protocol.sampling_unit,
                                                temperature=protocol.temperature*protocol.temperature_unit,
                                                pressure=protocol.pressure*protocol.pressure_unit,
                                                restart=True
                                               )

elif engine_query == 'SOMD':
    eq_protocol = BSS.Protocol.FreeEnergy(timestep=protocol.timestep*protocol.timestep_unit,
                                          num_lam=num_lambda,
                                          temperature=protocol.temperature*protocol.temperature_unit,
                                          runtime=(protocol.eq_runtime*2)*protocol.eq_runtime_unit,
                                          pressure=protocol.pressure*protocol.pressure_unit,
                                          )
    freenrg_protocol = BSS.Protocol.FreeEnergy(timestep=protocol.timestep*protocol.timestep_unit,
                                                num_lam=protocol.num_lambda,
                                                runtime=protocol.sampling*protocol.sampling_unit,
                                                temperature=protocol.temperature*protocol.temperature_unit,
                                                pressure=protocol.pressure*protocol.pressure_unit,
                                               )

# now set up the MD directories.
print(f"setting up FEP run in {workdir}...")

if engine_query == 'AMBER' or engine_query == 'GROMACS':

    # set up for each the bound and the free leg
    for leg, syst in zip(["bound", "free"], [system_bound, system_free]):

        BSS.FreeEnergy.Relative(
            syst,
            min_protocol,
            engine=f"{engine_query}",
            work_dir=f"{workdir}/{leg}_0/min"
        )

        BSS.FreeEnergy.Relative(
            syst,
            heat_protocol,
            engine=f"{engine_query}",
            work_dir=f"{workdir}/{leg}_0/heat",
            extra_options={}
        )

        BSS.FreeEnergy.Relative(
            syst,
            eq_protocol,
            engine=f"{engine_query}",
            work_dir=f"{workdir}/{leg}_0/eq",
            extra_options={}
        )

        BSS.FreeEnergy.Relative(
            syst,
            freenrg_protocol,
            engine=f"{engine_query}",
            work_dir=f"{workdir}/{leg}_0",
            extra_options={}
        )

if engine_query == "SOMD":
    for leg, syst in zip(["bound", "free"], [system_bound, system_free]):

        BSS.FreeEnergy.Relative(
            syst,
            eq_protocol,
            engine=f"{engine_query}",
            work_dir=f"{workdir}/{leg}_0/eq",
            extra_options={'minimise': True, 'minimise maximum iterations': protocol.min_steps, 'equilibrate': False}
        )

        BSS.FreeEnergy.Relative(
            syst,
            freenrg_protocol,
            engine=f"{engine_query}",
            work_dir=f"{workdir}/{leg}_0",
            extra_options={'minimise': False, 'equilibrate': False}
        )

# default folder is with no integer.
# for the sake of analysis , doesnt matter as finds folders w names of leg
more_repeats = list(range(1, protocol.repeats))

print(f"there are {protocol.repeats} folder(s) being made for each leg...")
for r in more_repeats:
    for leg in ["bound", "free"]:
        copy_tree(f"{workdir}/{leg}_0", f"{workdir}/{leg}_{r}")

print("done.")
