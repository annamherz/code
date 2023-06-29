import BioSimSpace as BSS
from BioSimSpace import _Exceptions
import sys
import os
import glob
import csv
from cmath import nan

BSS.setVerbose(True)

# box size
box_axis_length_list = [10, 15, 20, 25, 30, 35, 40]
box_axis_unit_query = "angstrom"
box_axis_unit = BSS.Units.Length.angstrom

boxtype_query_list = ["cubic", "truncatedoctahedron"]
# dictionary of functions
boxtype_dict = {
    "cubic": BSS.Box.cubic,
    "truncatedoctahedron": BSS.Box.truncatedOctahedron,
    "octahedral": BSS.Box.truncatedOctahedron,
}


# ligand ff
ligff_query = "sage"


def lig_paramaterise(molecule, ligff_query):
    ligff_dict = {
        "gaff": BSS.Parameters.gaff,
        "gaff2": BSS.Parameters.gaff2,
        "sage": BSS.Parameters.parameterise,
        "parsely": BSS.Parameters.parameterise,
    }
    func = ligff_dict[ligff_query]
    if ligff_query == "gaff" or ligff_query == "gaff2":
        return func(molecule)
    if ligff_query == "sage":
        return func(molecule, "openff_unconstrained-2.0.0")
    if ligff_query == "parsely":
        return func(molecule, "openff_unconstrained-1.3.0")


# protein ff
protff_query = "ff14SB"

# solvent ff
solvent_query = "TIP3P"


# set files to check
for protein in ["mcl1", "tyk2"]:
    ligands_folder = f"../{protein}/ligands"
    protein_file = f"../{protein}/protein/{protein}_parameterised"

    # make folder for solvated strucrues
    solvated_folder = f"solvated/{protein}"
    bound_folder = f"{solvated_folder}/sys"
    free_folder = f"{solvated_folder}/lig"

    # folders that need to be grown
    folders = [solvated_folder, bound_folder, free_folder]
    for folder in folders:
        if not os.path.isdir(folder):
            os.mkdir(folder)
            print(f"made dir {folder}")

    # folders for each box type and size
    for leg in [bound_folder, free_folder]:
        for boxtype in boxtype_query_list:
            folder = f"{leg}/{boxtype}"
            if not os.path.isdir(folder):
                os.mkdir(folder)
                print(f"made dir {folder}")
            for length in box_axis_length_list:
                folder = f"{leg}/{boxtype}/{str(length)}"
                if not os.path.isdir(folder):
                    os.mkdir(folder)
                    print(f"made dir {folder}")

    ligands = []
    ligand_files = sorted(os.listdir(ligands_folder))
    for lig in ligand_files:
        name = lig.split(".")[0]
        ligands.append(name)
    print(ligands)
    print(len(ligands))

    prot_wat = BSS.IO.readMolecules([f"{protein_file}.rst7", f"{protein_file}.prm7"])

    # ligand_dict make
    ligand_dict = {}
    system_dict = {}

    for lig in ligands:
        try:
            ligand = BSS.IO.readMolecules(f"{ligands_folder}/{lig}.sdf")[0]
            # paramaterise the ligand
            lig_p = lig_paramaterise(ligand, ligff_query).getMolecule()
            # Combine protein, ligand and crystallographic waters.
            system = lig_p + prot_wat
            ligand_dict[lig] = lig_p
            system_dict[lig] = system
            print(f"parameterised {lig}")
        except:
            print(f"could not parameterise {lig}")

    with open(f"{protein}_solvated_molecules.csv", "w") as file:
        writer = csv.writer(file, delimiter=";")
        writer.writerow(
            [
                "ligand",
                "protein",
                "leg",
                "boxtype",
                "boxlength",
                "box",
                "angles",
                "nMolecules",
            ]
        )

        for boxtype_query in boxtype_query_list:
            for box_axis_length in box_axis_length_list:
                for lig in ligands:
                    lig_p = ligand_dict[lig]
                    system = system_dict[lig]
                    legs_mols = [lig_p, system]
                    legs = ["lig", "sys"]

                    # zip together the molecules in that leg with the name for that leg
                    for leg, leg_mol in zip(legs, legs_mols):
                        try:
                            # define the box sizes based on the sizes of what is being solvated
                            box_min, box_max = leg_mol.getAxisAlignedBoundingBox()
                            # calcualte the minimum box size needed
                            box_size = [y - x for x, y in zip(box_min, box_max)]
                            # add the user defined box size around the min system size
                            box_sizes = [
                                x + int(box_axis_length) * box_axis_unit
                                for x in box_size
                            ]

                            # Solvate based on the boxtype query
                            # this also adds ions to balance the charge
                            boxtype_func = boxtype_dict[boxtype_query]
                            print(
                                f"Solvating {leg} for {lig}, {box_axis_length} {box_axis_unit_query} {boxtype_query}..."
                            )
                            box, angles = boxtype_func(max(box_sizes))
                            leg_mol_solvated = BSS.Solvent.solvate(
                                solvent_query,
                                molecule=leg_mol,
                                box=box,
                                angles=angles,
                                ion_conc=0.15,
                            )

                            nmols = leg_mol_solvated.nMolecules()

                            BSS.IO.saveMolecules(
                                f"{solvated_folder}/{leg}/{boxtype_query}/{str(box_axis_length)}/{lig}_solv",
                                leg_mol_solvated,
                                ["PRM7", "RST7", "PDB"],
                            )

                        except:
                            print(f"error when solvating.")
                            box = nan
                            angles = nan
                            nmols = nan

                        writer.writerow(
                            [
                                lig,
                                protein,
                                leg,
                                boxtype_query,
                                box_axis_length,
                                box,
                                angles,
                                nmols,
                            ]
                        )
