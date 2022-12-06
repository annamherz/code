# FEP prep

# Import libraries
import BioSimSpace as BSS
from BioSimSpace import _Exceptions
import os
import sys
from distutils.dir_util import copy_tree

try:
    import pipeline
except:
    print("adding code to the pythonpath...")
    code = '/home/anna/Documents/code/python'
    if code not in sys.path:
        sys.path.insert(1, code)
    import pipeline

from pipeline.fepprep import *

# Values from the corresponding bash arrays.
# The ligands for which the transformation is to be carried out.
trans = sys.argv[1]
lig_1 = trans.split('~')[0]
lig_2 = trans.split('~')[1]
# engine
engine_query = str(sys.argv[2]).upper()
# number of lambda windows
num_lambda = int(sys.argv[3])

main_dir = os.environ["MAINDIRECTORY"]
# ref the network.dat which was set in overall bash script
net_file = os.environ["net_file"]
# ref the protocol.dat which was set in overall bash script
prot_file = os.environ["prot_file"]
prep_dir = f"{main_dir}/prep"  # define lig prep location
# define trans dir
workdir = f"{main_dir}/outputs/{engine_query}/{lig_1}~{lig_2}"

# check if the trans, eng and window are in the network file
# found = False
# with open(net_file, "r") as lambdas_file:
#     reader = csv.reader(lambdas_file, delimiter=" ")
#     for row in reader:
#         if (row[0] == lig_1 and row[1] == lig_2) or (row[1] == lig_1 and row[0] == lig_2):
#             if int(row[2]) == num_lambda:
#                 if str(row[-1]).upper() == engine_query:
#                     found = True

# if not found:
#     raise NameError(
#         f"The perturbation {trans} (or the reverse) with {num_lambda} windows using {engine_query} was not found in network.dat.")

if not isinstance(num_lambda, int):
    raise NameError(
        f"Number of lambda windows must be an integer.")

if engine_query not in ["SOMD", "GROMACS", "AMBER"]:
    raise NameError("Input MD engine not found/recognised. Please use any of ['SOMD', 'GROMACS', 'AMBER']"
                    + "on the end of the line of the transformation in network.dat ")
# read the protocol file to get all the requested parameters
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

# get the requested runtime.
runtime_query = search_dict["sampling"].rstrip().replace(" ", "").split("=")[-1].split("*")[0]
try:
    runtime_query = int(runtime_query)
except ValueError:
    raise NameError("Input runtime value not supported. Please use an integer"
                    + " on the seventh line of protocol.dat in the shape of (e.g.):\nsampling = 2*ns")

# make sure user has set ns or ps.
runtime_unit_query = search_dict["sampling"].rstrip().replace(" ", "").split("=")[-1].split("*")[1]
if runtime_unit_query not in ["ns", "ps"]:
    raise NameError("Input runtime unit not supported. Please use 'ns' or 'ps'"
                    + " on the seventh line of protocol.dat in the shape of (e.g.):\nsampling = 2*ns")
if runtime_unit_query == "ns":
    runtime_unit = BSS.Units.Time.nanosecond
elif runtime_unit_query == "ps":
    runtime_unit = BSS.Units.Time.picosecond

# Check if HMR
hmr = search_dict["HMR"].rstrip().replace(" ", "").split("=")[-1]
if hmr == "True":
    hmr = True
elif hmr == "False":
    hmr = False
if not isinstance(hmr, bool):
    raise NameError("Hydrogen Mass Repartitioning must be of type bool."
                    + "on the eighth line of protocol.dat in the shape of (e.g.):\nHMR = True")

# get the number of repeats
repeats = int(search_dict["repeats"].rstrip().replace(" ", "").split("=")[-1])
if not isinstance(repeats, int):
    raise NameError("Repeats must be of type int."
                    + "on the ninth line of protocol.dat in the shape of (e.g.):\nrepeats = 3")

# Set temperature and pressure values
start_temp = BSS.Types.Temperature(0, "kelvin")
end_temp = BSS.Types.Temperature(300, "kelvin")
pressure_query = BSS.Types.Pressure(1, "bar")

# define minimisation steps
min_steps = 10000 

# Load equilibrated free inputs for both ligands. Complain if input not found. These systems already contain equil. waters.
# these will have been prepped previously using the ligprep or by other means
# create the system for each the free and the bound leg.
system_free = None
system_bound = None
# zip together the molecules in that leg with the name for that leg
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
    merged_trans = mergeLigands(ligand_1, ligand_2, engine_query)

    # Get equilibrated waters and waterbox information for both bound and free. Get all information from lambda==0
    # Following is work around because setBox() doesn't validate correctly boxes with lengths and angles
    system_1.removeMolecules(ligand_1)
    system_final = merged_trans + system_1

    if leg == "free":
        system_free = system_final
    if leg == "bound":
        system_bound = system_final

# TODO
# repartition the hydrogen masses
if hmr == True:
    print("repartitioning hydrogen masses for 4fs timestep...")
    timestep_query = int(4)
    if engine_query == "AMBER":
        system_free.repartitionHydrogenMass(factor=3)
        system_bound.repartitionHydrogenMass(factor=3)
    elif engine_query == "GROMACS":
        system_free.repartitionHydrogenMass(factor=4)
        system_bound.repartitionHydrogenMass(factor=4)
    elif engine_query == "SOMD":
        pass
elif hmr == False:
    timestep_query = int(2)
else:
    raise TypeError(
        "HMR (Hydrogen Mass Repartitioning) must be either True or False.")

timestep_unit = BSS.Units.Time.femtosecond

# incl eq time in setup? or just have at this default
eq_runtime_query = int(100)  # 100 ps NVT and NPT eq
eq_runtime_unit = BSS.Units.Time.picosecond
# someway to set eq so it also is like this in somd?

# define the free energy protocol with all this information.
if engine_query == 'AMBER' or engine_query == 'GROMACS':
    min_protocol = BSS.Protocol.FreeEnergyMinimisation(
        num_lam=num_lambda, steps=min_steps)
    heat_protocol = BSS.Protocol.FreeEnergyEquilibration(timestep=timestep_query*timestep_unit, num_lam=num_lambda,
                                                         runtime=eq_runtime_query*eq_runtime_unit, pressure=None,
                                                         temperature_start=start_temp, temperature_end=end_temp
                                                         )
    eq_protocol = BSS.Protocol.FreeEnergyEquilibration(timestep=timestep_query*timestep_unit, num_lam=num_lambda,
                                                       runtime=eq_runtime_query*eq_runtime_unit, temperature=end_temp,
                                                       pressure=pressure_query, restart=True
                                                       )
    freenrg_protocol = BSS.Protocol.FreeEnergy(timestep=timestep_query*timestep_unit, num_lam=num_lambda, temperature=end_temp,
                                               runtime=runtime_query*runtime_unit, pressure=pressure_query,
                                               restart=True
                                               )

elif engine_query == 'SOMD':
    eq_protocol = BSS.Protocol.FreeEnergy(timestep=timestep_query*timestep_unit, num_lam=num_lambda, temperature=end_temp,
                                          runtime=(eq_runtime_query*2)*eq_runtime_unit, pressure=pressure_query,
                                          )
    freenrg_protocol = BSS.Protocol.FreeEnergy(timestep=timestep_query*timestep_unit, num_lam=num_lambda, temperature=end_temp,
                                               runtime=runtime_query*runtime_unit, pressure=pressure_query,
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
            extra_options={'minimise': True, 'minimise maximum iterations': min_steps, 'equilibrate': False}
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
more_repeats = list(range(1, repeats))

print(f"there are {repeats} folder(s) being made for each leg...")
for r in more_repeats:
    for leg in ["bound", "free"]:
        copy_tree(f"{workdir}/{leg}_0", f"{workdir}/{leg}_{r}")

print("done.")
