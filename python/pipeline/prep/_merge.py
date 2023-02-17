import BioSimSpace as BSS
from BioSimSpace import _Exceptions

from ..utils._validate import *

class merge():
    """class of static methods for merging ligands and systems
    """

    @staticmethod
    def merge_ligands( ligand_0, ligand_1, engine_query):
        """Merges two ligands in preperation for FEP run.

        Args:
            ligand_0 (_type_): BSS molecule
            ligand_1 (_type_): BSS molecule
            engine_query (str): must be either AMBER, SOMD, GROMACS

        Returns:
            BioSimSpace._SireWrappers._molecule.Molecule: merged ligands as BSS object
        """

        engine = validate.engine(engine_query)

        # Align ligand2 on ligand1
        mapping = BSS.Align.matchAtoms(
            ligand_0, ligand_1, complete_rings_only=True)           
        inv_mapping = {v: k for k, v in mapping.items()}
        ligand_2_a = BSS.Align.rmsdAlign(ligand_1, ligand_0, inv_mapping)

        # Generate merged molecule.
        merged_ligands = BSS.Align.merge(
            ligand_0, ligand_2_a, mapping, allow_ring_breaking=True)

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
    def merge_system(system0, system1, engine_query):
        """merges to BSS systems for FEP.

        Args:
            system0 (BioSimSpace._SireWrappers._system.System): system for the ligand at lambda 0. Is the base system.
            system1 (BioSimSpace._SireWrappers._system.System): system from which the ligand at lmabda 1 will be obtained
            engine_query (BioSimSpace.FreeEnergy.engines()): an engine allowed for FEP

        Raises:
            _Exceptions.AlignmentError: If the ligands could not be aligned

        Returns:
            BioSimSpace._SireWrappers._system.System: system0 with the ligand replaced by the merged ligand
        """

        engine = validate.engine(engine_query)
        system_0 = validate.system(system0)
        system_1 = validate.system(system1)

        ligand_0 = merge.extract_ligand(system_0)
        ligand_1 = merge.extract_ligand(system_1)

        if ligand_0 and ligand_1:
            pass
        else:
            raise _Exceptions.AlignmentError(
                "Could not extract ligands from input systems. Check that your ligands/proteins are properly prepared!")

        # merge the ligands based on the engine.
        print("mapping, aligning and merging the ligands...")
        merged_trans = merge.merge_ligands(ligand_0, ligand_1, engine)

        # put the merged ligand into the system_0 in place of ligand_0
        system_0.removeMolecules(ligand_0)
        system_final = merged_trans + system_0

        return system_final


    @staticmethod
    def atom_mappings(system0, system1, engine_query):

        engine = validate.engine(engine_query)
        system_0 = validate.system(system0)
        system_1 = validate.system(system1)

        ligand_0 = merge.extract_ligand(system_0)
        ligand_1 = merge.extract_ligand(system_1)

        if ligand_0 and ligand_1:
            pass
        else:
            raise _Exceptions.AlignmentError(
                "Could not extract ligands from input systems. Check that your ligands/proteins are properly prepared!")

        # Align ligand2 on ligand1
        mapping = BSS.Align.matchAtoms(
            ligand_0, ligand_1, engine=engine, complete_rings_only=True)

        return (ligand_0.getAtoms(), ligand_1.getAtoms(), mapping)