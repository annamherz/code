import BioSimSpace as BSS
from distutils.dir_util import copy_tree, remove_tree

from ..utils._validate import *
from ._merge import *
from ._ligprep import *
from ._equilibrate import *

class fepprep():
    """class for fepprep
        makes all the protocols and writes the folders
    """

    def __init__(self, free_system=None, bound_system=None, protocol=None):

        # instantiate the class with the system and pipeline_protocol
        if free_system:
            self._merge_free_system = validate.system(free_system).copy()
        else:
            self._merge_free_system = None
            self._free_system_0 = None
            self._free_system_1 = None
            print("please add a system for free with lig0 and free with lig1 and merge these")

        if bound_system:
            self._merge_bound_system = validate.system(bound_system).copy()
        else:
            self._merge_bound_system = None
            self._bound_system_0 = None
            self._bound_system_1 = None
            print("please add a system for bound with lig0 and bound with lig1 and merge these")

        self._pipeline_protocol = validate.pipeline_protocol(protocol, fepprep=True)
        # generate the BSS protocols from the pipeline protocol
        fepprep._generate_bss_protocols(self)

    def add_system(self, system, free_bound=None, start_end=None):

        if free_bound not in ["free", "bound"]:
            raise ValueError("free_bound must be free or bound.")
        if start_end not in ["start", "end"]:
            raise ValueError("start_end must be start or end.")
                
        if free_bound == "free" and start_end == "start":
            self._free_system_0 = validate.system(system).copy()
        if free_bound == "free" and start_end == "end":
            self._free_system_1 = validate.system(system).copy()
        if free_bound == "bound" and start_end == "start":
            self._bound_system_0 = validate.system(system).copy()
        if free_bound == "bound" and start_end == "end":
            self._bound_system_1 = validate.system(system).copy()


    def _merge_systems(self, align_to):
            
        free_system = merge.merge_system(self._free_system_0, self._free_system_1, **{"align to": align_to})
        bound_system = merge.merge_system(self._bound_system_0, self._bound_system_1, **{"align to": align_to})
        
        return free_system, bound_system
    
    def merge_systems(self, align_to="lig0"):

        self._merge_free_system, self._merge_bound_system = self._merge_systems(align_to)

    def _generate_bss_protocols(self):

        protocol = self._pipeline_protocol

        if protocol.engine() == 'AMBER' or protocol.engine() == 'GROMACS':
            
            min_protocol = BSS.Protocol.FreeEnergyMinimisation(
                                                            num_lam=protocol.num_lambda(),
                                                            steps=protocol.min_steps(),
                                                            )
            heat_protocol = BSS.Protocol.FreeEnergyEquilibration(timestep=protocol.timestep()*protocol.timestep_unit(),
                                                                num_lam=protocol.num_lambda(),
                                                                runtime=protocol.eq_runtime()*protocol.eq_runtime_unit(),
                                                                pressure=None,
                                                                temperature_start=protocol.start_temperature()*protocol.temperature_unit(),
                                                                temperature_end=protocol.end_temperature()*protocol.temperature_unit(),
                                                                hmr_factor=protocol.hmr_factor()
                                                                )
            eq_protocol = BSS.Protocol.FreeEnergyEquilibration(timestep=protocol.timestep()*protocol.timestep_unit(),
                                                            num_lam=protocol.num_lambda(),
                                                            runtime=protocol.eq_runtime()*protocol.eq_runtime_unit(),
                                                            temperature=protocol.temperature()*protocol.temperature_unit(),
                                                            pressure=protocol.pressure()*protocol.pressure_unit(),
                                                            restart=True,
                                                            hmr_factor=protocol.hmr_factor()
                                                            )
            freenrg_protocol = BSS.Protocol.FreeEnergy(timestep=protocol.timestep()*protocol.timestep_unit(),
                                                        num_lam=protocol.num_lambda(),
                                                        runtime=protocol.sampling()*protocol.sampling_unit(),
                                                        temperature=protocol.temperature()*protocol.temperature_unit(),
                                                        pressure=protocol.pressure()*protocol.pressure_unit(),
                                                        restart=True,
                                                        hmr_factor=protocol.hmr_factor()
                                                    )

        elif protocol.engine() == 'SOMD':

            min_protocol = None
            heat_protocol = None

            eq_protocol = BSS.Protocol.FreeEnergy(timestep=protocol.timestep()*protocol.timestep_unit(),
                                                num_lam=protocol.num_lambda(),
                                                temperature=protocol.temperature()*protocol.temperature_unit(),
                                                runtime=(protocol.eq_runtime()*2)*protocol.eq_runtime_unit(),
                                                pressure=protocol.pressure()*protocol.pressure_unit(),
                                                hmr_factor=protocol.hmr_factor()
                                                )
            freenrg_protocol = BSS.Protocol.FreeEnergy(timestep=protocol.timestep()*protocol.timestep_unit(),
                                                        num_lam=protocol.num_lambda(),
                                                        runtime=protocol.sampling()*protocol.sampling_unit(),
                                                        temperature=protocol.temperature()*protocol.temperature_unit(),
                                                        pressure=protocol.pressure()*protocol.pressure_unit(),
                                                        hmr_factor=protocol.hmr_factor()
                                                       )

        # set the new protocols to self as well
        self._min_protocol = min_protocol
        self._heat_protocol = heat_protocol
        self._eq_protocol = eq_protocol
        self._freenrg_protocol = freenrg_protocol

    def prep_system_middle(self, pmemd_path, work_dir=None):

        if self._pipeline_protocol.fepprep() == "middle":
            # Solvate and run each the bound and the free system.
            legs_mols, legs = [self._merge_free_system, self._merge_bound_system], ["lig", "sys"]

            # zip together the molecules in that leg with the name for that leg
            for leg, leg_mol in zip(legs, legs_mols):
                print(f"carrying out for {leg}")
                leg_equil_final = minimise_equilibrate_leg(leg_mol, "AMBER", pmemd_path, lig_fep="fepprep", work_dir=work_dir)
                if leg == "lig":
                    self._merge_free_system = leg_equil_final
                if leg == "sys":
                    self._merge_bound_system = leg_equil_final
        
        return self._merge_free_system, self._merge_bound_system 

    def _generate_folders(self, system_free, system_bound, work_dir):
        
        protocol = self._pipeline_protocol
        min_protocol = self._min_protocol
        heat_protocol = self._heat_protocol
        eq_protocol = self._eq_protocol
        freenrg_protocol = self._freenrg_protocol

        print(f"setting up FEP run in {work_dir}...")

        if protocol.engine() == 'AMBER' or protocol.engine() == 'GROMACS':

            # set up for each the bound and the free leg
            for leg, system in zip(["bound", "free"], [system_bound, system_free]):

                # repartition the hydrogen masses, so only needs to be done once during the setup
                if protocol.hmr() == True:
                    print(f"repartitioning hydrogen masses for 4fs timestep for {leg}...")
                    if protocol.hmr_factor() == "auto":
                        print("using default factors...")
                        if protocol.engine() == "AMBER":
                            system.repartitionHydrogenMass(factor=3)
                        elif protocol.engine() == "GROMACS":
                            system.repartitionHydrogenMass(factor=4)
                    else:
                        print(f"using {protocol.hmr_factor()} as a factor...")
                        system.repartitionHydrogenMass(factor=protocol.hmr_factor())
                elif protocol.hmr() == False:
                    pass

                BSS.FreeEnergy.Relative(
                    system,
                    min_protocol,
                    engine=f"{protocol.engine()}",
                    work_dir=f"{work_dir}/{leg}_0/min"
                )

                BSS.FreeEnergy.Relative(
                    system,
                    heat_protocol,
                    engine=f"{protocol.engine()}",
                    work_dir=f"{work_dir}/{leg}_0/heat",
                    extra_options={}
                )

                BSS.FreeEnergy.Relative(
                    system,
                    eq_protocol,
                    engine=f"{protocol.engine()}",
                    work_dir=f"{work_dir}/{leg}_0/eq",
                    extra_options={}
                )

                BSS.FreeEnergy.Relative(
                    system,
                    freenrg_protocol,
                    engine=f"{protocol.engine()}",
                    work_dir=f"{work_dir}/{leg}_0",
                    extra_options={}
                )

        if protocol.engine() == "SOMD":
            for leg, system in zip(["bound", "free"], [system_bound, system_free]):

                BSS.FreeEnergy.Relative(
                    system,
                    eq_protocol,
                    engine=f"{protocol.engine()}",
                    work_dir=f"{work_dir}/{leg}_0/eq",
                    extra_options={'minimise': True, 'minimise maximum iterations': protocol.min_steps, 'equilibrate': False}
                )

                BSS.FreeEnergy.Relative(
                    system,
                    freenrg_protocol,
                    engine=f"{protocol.engine()}",
                    work_dir=f"{work_dir}/{leg}_0",
                    extra_options={'minimise': False, 'equilibrate': False}
                )


    def generate_folders(self, work_dir):

        work_dir = validate.folder_path(work_dir, create=True)

        if self._pipeline_protocol.fepprep() == "both":
            ligs = ["lig0", "lig1"]
            for lig in ligs:
                free_system, bound_system = self._merge_systems(align_to=lig)
                self._generate_folders(free_system, bound_system, f"{work_dir}/{lig}")

            lambdas_list = self._min_protocol.getLambdaValues()
            middle_index=len(lambdas_list)//2        
            first_half=lambdas_list[:middle_index]
            sec_half=lambdas_list[middle_index:]
            
            # copy files to main folder
            print("copying generated folders for the endstates into a combined folder, so first half is lig0 and second half is lig1")
            for lig, lam_list in zip(ligs, [first_half, sec_half]):
                for leg in ["bound", "free"]:
                    for part in ["min/","heat/", "eq/", ""]:
                        for lam in lam_list:
                            copy_tree(f"{work_dir}/{lig}/{leg}_0/{part}lambda_{lam:.4f}", f"{work_dir}/{leg}_0/{part}lambda_{lam:.4f}")
                
                # remove the dir
                print(f"removing directory for {lig} as copied...")
                remove_tree(f"{work_dir}/{lig}")

        else:
            if not self._merge_free_system or not self._merge_bound_system:
                print("no merged systems, merging....")
                self.merge_systems()
            self._generate_folders(self._merge_free_system, self._merge_bound_system, work_dir)

        # default folder is with no integer.
        # for the sake of analysis , doesnt matter as finds folders w names of leg
        more_repeats = list(range(1, self._pipeline_protocol.repeats()))

        print(f"there are {self._pipeline_protocol.repeats()} folder(s) being made for each leg...")
        for r in more_repeats:
            for leg in ["bound", "free"]:
                copy_tree(f"{work_dir}/{leg}_0", f"{work_dir}/{leg}_{r}")

        print("done.")
