import BioSimSpace as BSS
from BioSimSpace import _Exceptions
from BioSimSpace.Units.Length import angstrom as _angstrom
from distutils.dir_util import copy_tree

from ..utils._validate import *



def run_process(system, protocol, engine="AMBER", pmemd_path=None, work_dir=None):
    """
    Given a solvated system (BSS object) and BSS protocol, run a process workflow with either 
    AMBER or GROMACS. Returns the processed system.
    """

    # Create the process passing a working directory.
    if engine == "AMBER":
        if pmemd_path:
            process = BSS.Process.Amber(
                system, protocol, work_dir=work_dir, exe=pmemd_path)
        else:
            process = BSS.Process.Amber(
                system, protocol, work_dir=work_dir)

    elif engine == "GROMACS":
        process = BSS.Process.Gromacs(
            system, protocol, work_dir=work_dir)

    # Start the process, wait for it to exit
    process.start()
    process.wait()

    # Check for errors.
    if process.isError():
        print(process.stdout())
        print(process.stderr())
        raise _Exceptions.ThirdPartyError("The process exited with an error!")

    # If it worked, try to get the system. No need to block, since it's already finished.
    system = process.getSystem()

    return system


class ligprep():
    """class to store lig prep functions
    """


    def lig_paramaterise(molecule, ligff_query):
        # dicitonary of functions available
        validate.lig_ff(ligff_query)
        
        return BSS.Parameters.parameterise(molecule, ligff_query)
            

    def minimise_equilibrate_leg(system_solvated, leg, engine="AMBER", pmemd=None):  # times at the top
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
        minimised1 = run_process(system_solvated, protocol_min_rest, engine, pmemd)
        minimised2 = run_process(minimised1, protocol_min, engine, pmemd)
        equil1 = run_process(minimised2, protocol_nvt_sol, engine, pmemd)
        equil2 = run_process(equil1, protocol_nvt_backbone, engine, pmemd)
        equil3 = run_process(equil2, protocol_nvt, engine, pmemd)
        equil4 = run_process(equil3, protocol_npt_heavy, engine, pmemd)
        equil5 = run_process(equil4, protocol_npt_heavy_lighter, engine, pmemd)
        sys_equil_fin = run_process(equil5, protocol_npt, engine, pmemd)

        return sys_equil_fin

    def minimum_solvation(system, solvent, box_type, box_edges, box_edges_unit=None, verbose=True):
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


class fepprep():
    """class for fepprep
    """

    def __init__(self, free_system, bound_system, protocol):
        # instantiate the class with the system and pipeline_protocol
        self._free_system = validate.system(free_system)
        self._bound_system = validate.system(bound_system)
        self._pipeline_protocol = validate.pipeline_protocol(protocol, fepprep=True)
        # generate the BSS protocols from the pipeline protocol
        fepprep._generate_bss_protocols(self)


    def _generate_bss_protocols(self):

        protocol = self._pipeline_protocol

        if protocol.engine == 'AMBER' or protocol.engine == 'GROMACS':
            min_protocol = BSS.Protocol.FreeEnergyMinimisation(num_lam=protocol.num_lambda,
                                                            steps=protocol.min_steps
                                                            )
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

        elif protocol.engine == 'SOMD':
            eq_protocol = BSS.Protocol.FreeEnergy(timestep=protocol.timestep*protocol.timestep_unit,
                                                num_lam=protocol.num_lambda,
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

        # set the new protocols to self as well
        self._min_protocol = min_protocol
        self._heat_protocol = heat_protocol
        self._eq_protocol = eq_protocol
        self._freenrg_protocol = freenrg_protocol


    def generate_folders(self, work_dir):

        work_dir = validate.folder_path(work_dir)

        system_free = self._free_system
        system_bound = self._bound_system
        protocol = self._pipeline_protocol
        min_protocol = self._min_protocol
        heat_protocol = self._heat_protocol
        eq_protocol = self._eq_protocol
        freenrg_protocol = self._freenrg_protocol


        print(f"setting up FEP run in {work_dir}...")

        if protocol.engine == 'AMBER' or protocol.engine == 'GROMACS':

            # set up for each the bound and the free leg
            for leg, system in zip(["bound", "free"], [system_bound, system_free]):

                # repartition the hydrogen masses
                if protocol.hmr == True:
                    print(f"repartitioning hydrogen masses for 4fs timestep for {leg}...")
                    if protocol.engine == "AMBER":
                        system.repartitionHydrogenMass(factor=3)
                    elif protocol.engine == "GROMACS":
                        system.repartitionHydrogenMass(factor=4)
                elif protocol.hmr == False:
                    pass

                BSS.FreeEnergy.Relative(
                    system,
                    min_protocol,
                    engine=f"{protocol.engine}",
                    work_dir=f"{work_dir}/{leg}_0/min"
                )

                BSS.FreeEnergy.Relative(
                    system,
                    heat_protocol,
                    engine=f"{protocol.engine}",
                    work_dir=f"{work_dir}/{leg}_0/heat",
                    extra_options={}
                )

                BSS.FreeEnergy.Relative(
                    system,
                    eq_protocol,
                    engine=f"{protocol.engine}",
                    work_dir=f"{work_dir}/{leg}_0/eq",
                    extra_options={}
                )

                BSS.FreeEnergy.Relative(
                    system,
                    freenrg_protocol,
                    engine=f"{protocol.engine}",
                    work_dir=f"{work_dir}/{leg}_0",
                    extra_options={}
                )

        if protocol.engine == "SOMD":
            for leg, system in zip(["bound", "free"], [system_bound, system_free]):

                BSS.FreeEnergy.Relative(
                    system,
                    eq_protocol,
                    engine=f"{protocol.engine}",
                    work_dir=f"{work_dir}/{leg}_0/eq",
                    extra_options={'minimise': True, 'minimise maximum iterations': protocol.min_steps, 'equilibrate': False}
                )

                BSS.FreeEnergy.Relative(
                    system,
                    freenrg_protocol,
                    engine=f"{protocol.engine}",
                    work_dir=f"{work_dir}/{leg}_0",
                    extra_options={'minimise': False, 'equilibrate': False}
                )

        # default folder is with no integer.
        # for the sake of analysis , doesnt matter as finds folders w names of leg
        more_repeats = list(range(1, protocol.repeats))

        print(f"there are {protocol.repeats} folder(s) being made for each leg...")
        for r in more_repeats:
            for leg in ["bound", "free"]:
                copy_tree(f"{work_dir}/{leg}_0", f"{work_dir}/{leg}_{r}")

        print("done.")


    @staticmethod
    def merge_ligands( ligand_0, ligand_1, engine_query):
        """Merges two ligands in preperation for FEP run.

        Args:
            ligand_0 (_type_): BSS molecule
            ligand_1 (_type_): BSS molecule
            engine_query (str): must be either AMBER, SOMD, GROMACS

        Returns:
            _type_: merged ligands as BSS object
        """

        engine = validate.engine(engine_query)

        # Align ligand2 on ligand1
        mapping = BSS.Align.matchAtoms(
            ligand_0, ligand_1, engine=engine, complete_rings_only=True)
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

        ligand_0 = fepprep.extract_ligand(system_0)
        ligand_1 = fepprep.extract_ligand(system_1)

        if ligand_0 and ligand_1:
            pass
        else:
            raise _Exceptions.AlignmentError(
                "Could not extract ligands from input systems. Check that your ligands/proteins are properly prepared!")

        # merge the ligands based on the engine.
        print("mapping, aligning and merging the ligands...")
        merged_trans = fepprep.merge_ligands(ligand_0, ligand_1, engine)

        # put the merged ligand into the system_0 in place of ligand_0
        system_0.removeMolecules(ligand_0)
        system_final = merged_trans + system_0

        return system_final
