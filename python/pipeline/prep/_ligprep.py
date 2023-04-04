import BioSimSpace as BSS
from BioSimSpace.Units.Length import angstrom as _angstrom

from ..utils._validate import *

class ligprep():
    """class to store lig prep functions.
    """

    def __init__(self):
        pass

    @staticmethod
    def lig_paramaterise(molecule, ligff_query):
        # dicitonary of functions available
        validate.lig_ff(ligff_query)

        param_molecule = BSS.Parameters.parameterise(molecule, ligff_query).getMolecule()
        
        return param_molecule
            

    @staticmethod
    def minimum_solvation(system, solvent, box_type, box_edges, box_edges_unit="angstrom", verbose=True):
        """
        Default solvation for minimum size box.
        """

        # validate inputs
        try:
            solvent = validate.solvent_ff(solvent)
            box_edges = validate.integer(box_edges)
            box_edges_unit = validate.box_edges_unit(box_edges_unit)
            box_type = validate.box_type(box_type)
        except Exception as e:
            print(f"The provided arguments could not be validated.\n Exception is:\n {e}")    
                            
        # type of solvation models available in BSS
        boxtype_dict = {"cubic": BSS.Box.cubic,
                    "truncatedOctahedron": BSS.Box.truncatedOctahedron,
                    "octahedral": BSS.Box.truncatedOctahedron}

        # define the box sizes based on the sizes of what is being solvated
        box_min, box_max = system.getAxisAlignedBoundingBox()
        # calcualte the minimum box size needed
        box_size = [y - x for x, y in zip(box_min, box_max)]
        # add the user defined box size around the min system size
        box_sizes = [x + int(box_edges) * box_edges_unit for x in box_size]

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

        # Solvate based on the box_type query
        # this also adds ions to balance the charge
        boxtype_func = boxtype_dict[box_type]
        box, angles = boxtype_func(max(box_sizes))
        mol_solvated = BSS.Solvent.solvate(solvent, molecule=system,
                                            box=box, angles=angles, ion_conc=0.15)

        nmols = mol_solvated.nMolecules()

        if verbose == True:
            print(f"box dimensions for {box_edges} {box_edges_unit} {box_type} are:")
            print(f"box_min : {box_min}")
            print(f"box_max : {box_max}")
            print(f"box_size : {box_size}")
            print(f"box_sizes : {box_sizes}")
            print(f"with the final box : {box} with angles as : {angles}")
            print(f"The total no of molecules is : {nmols}")

        return mol_solvated
