import BioSimSpace as BSS
from BioSimSpace import _Exceptions

from ..utils._validate import *

class merge():
    """class of static methods for merging ligands and systems
    """

    @staticmethod
    def merge_ligands( ligand_0, ligand_1, **kwargs):
        """Merges two ligands in preperation for FEP run.

        Args:
            ligand_0 (BioSimSpace._SireWrappers._molecule.Molecule): BSS molecule
            ligand_1 (BioSimSpace._SireWrappers._molecule.Molecule): BSS molecule
            kwargs (dict): 'allow ring breaking' as True or False

        Returns:
            BioSimSpace._SireWrappers._molecule.Molecule: merged ligands as BSS object
        """

        # default ring breaking is not allowed
        complete_rings = True
        allow_ring_breaking = False
        allow_ring_size_change = False
        prune_perturbed_constraints=None
        prune_crossing_constraints=None
        align_to = "lig0" 
        scoring_function = "rmsd_align"
        prematch = {}

        for key,value in kwargs.items():
            key = key.replace("_","").replace(" ","").upper()
            if key == "COMPLETERINGSONLY":
                complete_rings = validate.boolean(value)
            if key == "ALLOWRINGBREAKING":
                allow_ring_breaking = validate.boolean(value)
            if key == "ALLOWRINGSIZECHANGE":
                allow_ring_size_change = validate.boolean(value)
            if key == "PRUNEPERTURBEDCONSTRAINTS":
                prune_perturbed_constraints = validate.boolean(value)
            if key == "PRUNECROSSINGRESTRAINTS":
                prune_crossing_constraints = validate.boolean(value)           
            if key == "ALIGNTO":
                align_to = validate.string(value)  
            if key == "SCORINGFUNCTION":
                scoring_function = validate.string(value) # rmsd, rmsd_align, rmsd_flex_align
            if key == "PREMATCH":
                prematch = validate.dictionary(value)

        # Align ligand2 on ligand1
        # get the mapping of ligand0 to atoms in ligand1
        mapping = BSS.Align.matchAtoms(
                                    ligand_0, ligand_1,
                                    scoring_function=scoring_function,
                                    complete_rings_only=complete_rings,
                                    prematch=prematch,
                                    prune_perturbed_constraints=prune_perturbed_constraints,
                                    prune_crossing_constraints=prune_crossing_constraints
                                    )           
        inv_mapping = {v: k for k, v in mapping.items()}

        # function for aligning depends on scoring function
        func_dict = {"rmsd_align": BSS.Align.rmsdAlign,
                     "rmsd_flex_align": BSS.Align.flexAlign}
        
        align_func = func_dict[scoring_function]

        if align_to == "lig0":
            # need inverse mapping to align
            # aligns atoms in first argument to atoms in second argument
            ligand_1_a = align_func(ligand_1, ligand_0, inv_mapping)
            # Generate merged molecule.
            merged_ligands = BSS.Align.merge(
                ligand_0, ligand_1_a, mapping,
                allow_ring_breaking=allow_ring_breaking,
                allow_ring_size_change=allow_ring_size_change
            )
            
        elif align_to == "lig1":
            # align ligand 0 to ligand 1
            ligand_0_a = align_func(ligand_0, ligand_1, mapping)
            # Generate merged molecule.
            merged_ligands = BSS.Align.merge(
                ligand_0_a, ligand_1, mapping, 
                allow_ring_breaking=allow_ring_breaking,
                allow_ring_size_change=allow_ring_size_change
            )
            
        else:
            raise ValueError("must align to 'lig0' or 'lig1'")

        return merged_ligands
    

    @staticmethod
    def extract_ligand(system):
        """extracts the ligand from a BSS system

        Args:
            system (BioSimSpace._SireWrappers._system.System): system containing a ligand

        Returns:
            BioSimSpace._SireWrappers._molecule.Molecule: ligand as a molecule object
        """

        system = validate.system(system)

        # Extract ligands. Do this based on nAtoms and nResidues, as sometimes
        # the order of molecules is switched, so we can't use index alone.
        ligand = None
        n_residues = [mol.nResidues() for mol in system]
        n_atoms = [mol.nAtoms() for mol in system]
        for i, (n_resi, n_at) in enumerate(zip(n_residues, n_atoms)):
            if n_resi == 1 and n_at > 5:
                ligand = system.getMolecule(i)
            else:
                pass
            if ligand:
                break
        
        return ligand

    @staticmethod
    def merge_system(system0=None, system1=None, **kwargs):
        """merges to BSS systems for FEP.

        Args:
            system0 (BioSimSpace._SireWrappers._system.System): system for the ligand at lambda 0. Is the base system.
            system1 (BioSimSpace._SireWrappers._system.System): system from which the ligand at lmabda 1 will be obtained
            kwargs (dict): extra options for merging the ligands

        Raises:
            _Exceptions.AlignmentError: If the ligands could not be aligned

        Returns:
            BioSimSpace._SireWrappers._system.System: system0 with the ligand replaced by the merged ligand
        """

        # default arguments
        align_to = "lig0" 

        # check kwargs
        for key,value in kwargs.items():
            if key == "align to":
                align_to = validate.string(value)

        system_0 = validate.system(system0).copy()
        system_1 = validate.system(system1).copy()

        ligand_0 = merge.extract_ligand(system_0)
        ligand_1 = merge.extract_ligand(system_1)

        if ligand_0 and ligand_1:
            pass
        else:
            raise _Exceptions.AlignmentError(
                "Could not extract ligands from input systems. Check that your ligands/proteins are properly prepared!")

        # merge the ligands based on the engine.
        print("mapping, aligning and merging the ligands...")
        merged_trans = merge.merge_ligands(ligand_0, ligand_1, **kwargs)

        if align_to == "lig0":
            # put the merged ligand into the system_0 in place of ligand_0
            system_0.removeMolecules(ligand_0)
            system_final = merged_trans + system_0

        elif align_to == "lig1":
            # put the merged ligand into the system_0 in place of ligand_0
            system_1.removeMolecules(ligand_1)
            system_final = merged_trans + system_1

        else:
            raise ValueError("must align to 'lig0' or 'lig1'")
        
        return system_final


    @staticmethod
    def atom_mappings(system0, system1, **kwargs):
        """get the atoms and mappings for ligands in two systems

        Args:
            system0 (BioSimSpace._SireWrappers.System): unmerged system at lambda 0.0
            system1 (BioSimSpace._SireWrappers.System): unmerged system at lambda 1.0

        Raises:
            _Exceptions.AlignmentError: can't extract ligands from input system

        Returns:
            dict: ligand_0 atoms, ligand_1 atoms, mapping dictionary)
        """

        prune_perturbed_constraints=None
        prune_crossing_constraints=None

        for key,value in kwargs.items():
            if key == "prune perturbed constraints":
                prune_perturbed_constraints = validate.boolean(value)
            if key == "prune crossing constraints":
                prune_crossing_constraints = validate.boolean(value)           

        system_0 = validate.system(system0)
        system_1 = validate.system(system1)

        ligand_0 = merge.extract_ligand(system_0)
        ligand_1 = merge.extract_ligand(system_1)

        if ligand_0 and ligand_1:
            pass
        else:
            raise _Exceptions.AlignmentError(
                "Could not extract ligands from input systems. Check that your ligands/proteins are properly prepared!")


        # Align ligand1 on ligand0
        mapping = BSS.Align.matchAtoms(
                                    ligand_0, ligand_1,
                                    complete_rings_only=True,
                                    prune_perturbed_constraints=prune_perturbed_constraints,
                                    prune_crossing_constraints=prune_crossing_constraints
                                    )      

        return (ligand_0.getAtoms(), ligand_1.getAtoms(), mapping)