

def mergeLigands(ligand_1, ligand_2, engine_query):
    """Merges two ligands in preperation for FEP run.

    Args:
        ligand_1 (_type_): BSS molecule
        ligand_2 (_type_): BSS molecule
        engine_query (_type_): must be either amber somd or groamcs

    Returns:
        _type_: merged ligands as BSS object
    """
    # Align ligand2 on ligand1
    mapping = BSS.Align.matchAtoms(
        ligand_1, ligand_2, engine=engine_query, complete_rings_only=True)
    inv_mapping = {v: k for k, v in mapping.items()}
    ligand_2_a = BSS.Align.rmsdAlign(ligand_2, ligand_1, inv_mapping)

    # Generate merged molecule.
    merged_ligands = BSS.Align.merge(
        ligand_1, ligand_2_a, mapping, allow_ring_breaking=True)

    return merged_ligands