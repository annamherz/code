# ligand prep

import BioSimSpace as BSS
from BioSimSpace import _Exceptions
from BioSimSpace.Units.Length import angstrom as _angstrom
import sys
import os

BSS.setVerbose(True)

# decide engine used for equilibration runs ("AMBER" or "GROMACS")
engine = "AMBER"

# set environs
# pmemd_path = os.environ["AMBERHOME"] + "/bin/pmemd.cuda"
pmemd_path = os.environ["amber"] + "/bin/pmemd.cuda"
main_dir = os.environ["MAINDIRECTORY"]
# main_dir = "/home/anna/Documents/benchmark/test_tyk2_benchmark_sage"
protein_file = os.environ["protein_file"]
ligands_folder = os.environ["ligands_folder"]
# protein_file = "/home/anna/Documents/benchmark/final_system_inputs/tyk2_parameterised"
# ligands_folder = "/home/anna/Documents/benchmark/final_system_inputs/tyk2_ligands"

################
# Ligand used
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
stream = open(f"{main_dir}/execution_model/protocol.dat", "r")
lines = stream.readlines()

# get value regardless of the order of the protocol.dat

# ligand ff
ligff_query = (lines[0].rstrip().replace(" ", "").split("=")[-1]).lower()

# protein ff
protff_query = lines[1].rstrip().replace(" ", "").split("=")[-1]

# solvent ff
solvent_query = lines[2].rstrip().replace(" ", "").split("=")[-1]

# box size
boxsize_query = lines[3].rstrip().replace(" ", "").split("=")[-1]
box_axis_length = boxsize_query.split("*")[0]
box_axis_unit_query = boxsize_query.split("*")[1]
        
# get the box type settings.
boxtype_query = (lines[4].rstrip().replace(" ", "").split("=")[-1]).lower()

stream.close()

#################
# FUNCTIONS

def lig_paramaterise(molecule, ligff_query):
    # dicitonary of functions available
    ligff_dict = {"gaff": BSS.Parameters.gaff,
                  "gaff2": BSS.Parameters.gaff2,
                  "sage": BSS.Parameters.parameterise,
                  "parsely": BSS.Parameters.parameterise}
    func = ligff_dict[ligff_query]
    if ligff_query == "gaff" or ligff_query == "gaff2":
        return func(molecule)
    if ligff_query == "sage":
        return func(molecule, 'openff_unconstrained-2.0.0')
    if ligff_query == "parsely":
        return func(molecule, 'openff_unconstrained-1.3.0')

# Throw error if ligand forcefield not available.
    if ligff_query not in ligff_dict:
        raise NameError(
            f"Force field not supported: {ligff_query}. Please use either of [GAFF1, GAFF2, Parsely, Sage]")
        
# define functions for the lig_paramateriserep
def runProcess(system, protocol, engine="AMBER", work_dir=None):
    """
    Given a solvated system (BSS object) and BSS protocol, run a process workflow with either 
    Sander (CPU) or pmemd.cuda (GPU). NPT is typically done with GPU to save computing time.
    Returns the processed system.
    """

    # Create the process passing a working directory.
    if engine == "AMBER":
        process = BSS.Process.Amber(
            system, protocol, work_dir=work_dir, exe=pmemd_path)  # need to define pmemd?

    elif engine == "GROMACS":
        process = BSS.Process.Gromacs(
            system, protocol, work_dir=work_dir)

    # Start the process.
    process.start()

    # Wait for the process to exit.
    process.wait()

    # Check for errors.
    if process.isError():
        print(process.stdout())
        print(process.stderr())
        raise _Exceptions.ThirdPartyError("The process exited with an error!")

    # If it worked, try to get the system. No need to block, since it's already finished.
    system = process.getSystem()

    return system


def lig_paramateriserep(system_solvated, leg):  # times at the top
    """
    Default protocols and running of the lig_paramateriserep.

    system_solvated is a BSS system

    leg is if lig or sys .
    """
    # define all the protocols

    # Minimisations
    # minimisation of solvent
    protocol_min_rest = BSS.Protocol.Minimisation(
        steps=10000,
        restraint="all"
    )
    # minimisation full system
    protocol_min = BSS.Protocol.Minimisation(
        steps=10000
    )

    # NVTs
    # NVT restraining all non solvent atoms
    protocol_nvt_sol = BSS.Protocol.Equilibration(
        runtime=400*BSS.Units.Time.picosecond,
        temperature_start=0*BSS.Units.Temperature.kelvin,
        temperature_end=300*BSS.Units.Temperature.kelvin,
        restraint="all"
    )
    # if leg is lig or sys
    if leg == "sys":
        back_rest = "backbone"
    elif leg == "lig":
        back_rest = "heavy"
    else:
        raise NameError("leg must be either 'sys' or 'lig'.")
    # NVT restraining all backbone/heavy atoms
    protocol_nvt_backbone = BSS.Protocol.Equilibration(
        runtime=400*BSS.Units.Time.picosecond,
        temperature=300*BSS.Units.Temperature.kelvin,
        restraint=back_rest
    )
    # NVT no restraints
    protocol_nvt = BSS.Protocol.Equilibration(
        runtime=400*BSS.Units.Time.picosecond,
        temperature_end=300*BSS.Units.Temperature.kelvin
    )

    # NPTs
    # NPT restraining all non solvent heavy atoms
    protocol_npt_heavy = BSS.Protocol.Equilibration(
        runtime=400*BSS.Units.Time.picosecond,
        pressure=1*BSS.Units.Pressure.atm,
        temperature=300*BSS.Units.Temperature.kelvin,
        restraint="heavy"
    )
    # NPT with gradual release of restraints
    protocol_npt_heavy_lighter = BSS.Protocol.Equilibration(
        runtime=400*BSS.Units.Time.picosecond,
        pressure=1*BSS.Units.Pressure.atm,
        temperature=300*BSS.Units.Temperature.kelvin,
        restraint="heavy",
        force_constant=5
    )
    # NPT no restraints
    protocol_npt = BSS.Protocol.Equilibration(
        runtime=1000*BSS.Units.Time.picosecond,
        pressure=1*BSS.Units.Pressure.atm,
        temperature=300*BSS.Units.Temperature.kelvin,
    )

    # run all the protocols
    minimised1 = runProcess(system_solvated, protocol_min_rest, engine=engine)
    minimised2 = runProcess(minimised1, protocol_min, engine=engine)
    equil1 = runProcess(minimised2, protocol_nvt_sol, engine=engine)
    equil2 = runProcess(equil1, protocol_nvt_backbone, engine=engine)
    equil3 = runProcess(equil2, protocol_nvt, engine=engine)
    equil4 = runProcess(equil3, protocol_npt_heavy, engine=engine)
    equil5 = runProcess(equil4, protocol_npt_heavy_lighter, engine=engine)
    sys_equil_fin = runProcess(equil5, protocol_npt, engine=engine)

    return sys_equil_fin

def min_solv(system, solvent, boxtype, box_axis_length, box_axis_unit_query=None):  # times at the top
    """
    Default solvation for minimum size box.
    """

    if not box_axis_unit_query:
        box_axis_unit = BSS.Units.Length.angstrom
    elif box_axis_unit_query.lower() == "nm" or box_axis_unit_query.lower() == "nanometer":
        box_axis_unit = BSS.Units.Length.nanometer
    elif box_axis_unit_query.lower() == "a" or box_axis_unit_query.lower() == "angstrom":
        box_axis_unit = BSS.Units.Length.angstrom
    else:
        raise NameError("Input unit not recognised. Please use any of ['nanometer', 'angstrom', 'nm', 'a']")
                         
    # type of solvation models available in BSS
    boxtype_dict = {"cubic": BSS.Box.cubic,
                "truncatedoctahedron": BSS.Box.truncatedOctahedron,
                "octahedral": BSS.Box.truncatedOctahedron}
    
    # Throw error if box type not available.
    if boxtype not in boxtype_dict:
        raise NameError("Input box type not recognised. Please use any of ['cubic', 'truncatedoctahedron', 'octahedral']") 
    
    # define the box sizes based on the sizes of what is being solvated
    box_min, box_max = system.getAxisAlignedBoundingBox()
    # calcualte the minimum box size needed
    box_size = [y - x for x, y in zip(box_min, box_max)]
    # add the user defined box size around the min system size
    box_sizes = [x + int(box_axis_length) * box_axis_unit for x in box_size]

    # for amber22 currently, eq fails if the overall box is less than 41 A
    # check the box size and adjust if needed
    min_size = 42*_angstrom   
    if max(box_sizes) < min_size:
        print(f"max box size {max(box_sizes)} is below the min size {min_size}. This will be replaced.")
        new_max_size = max(box_sizes) + (min_size - max(box_sizes))
        # replace the max box size in the list with the new max box size
        # so that in the next part, this max box size is used for the box type 
        for index, size in enumerate(box_sizes):
            if size == max(box_sizes):
                box_sizes[index] = new_max_size

    # Solvate based on the boxtype query
    # this also adds ions to balance the charge
    boxtype_func = boxtype_dict[boxtype]
    box, angles = boxtype_func(max(box_sizes))
    mol_solvated = BSS.Solvent.solvate(solvent, molecule=system,
                                           box=box, angles=angles, ion_conc=0.15)

    nmols = mol_solvated.nMolecules()

    print(f"box dimensions for {box_axis_length} {box_axis_unit_query} {boxtype} are:")
    print(f"box_min : {box_min}")
    print(f"box_max : {box_max}")
    print(f"box_size : {box_size}")
    print(f"box_sizes : {box_sizes}")
    print(f"with the final box : {box} with angles as : {angles}")
    print(f"The total no of molecules is : {nmols}")

    return mol_solvated


#################
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
