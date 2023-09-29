import BioSimSpace as BSS
from BioSimSpace.MD._md import _find_md_engines
from BioSimSpace import _Exceptions
import logging

from typing import Union, Optional

import pipeline
from ..utils._validate import *


def run_process(
    system: BSS._SireWrappers.System,
    protocol: str,
    engine: str = "AMBER",
    exe: Optional[str] = None,
    work_dir: Optional[str] = None,
) -> BSS._SireWrappers.System:
    """Given a solvated system (BSS object) and BSS protocol, run a process workflow with either
    AMBER or GROMACS. Returns the processed system.

    Args:
        system (BioSimSpace._SireWrappers.System): a BSS system
        protocol (pipeline.prep.pipeline_protocol): a pipeline protocol
        engine (str, optional): engine used for the process. Defaults to "AMBER".
        exe (str, optional): path to executable to run. Defaults to None.
        work_dir (str, optional): work dir for the BSS process. Defaults to None.

    Raises:
        _Exceptions.ThirdPartyError: If the process is not able to run.

    Returns:
        BioSimSpace._SireWrappers.System: the system after the process is completed.
    """

    eng_proc_dict = {
        "AMBER": BSS.Process.Amber,
        "GROMACS": BSS.Process.Gromacs,
        "SOMD": BSS.Process.Somd,
    }

    func = eng_proc_dict[engine]  # validate.engine(engine)

    if exe:
        pass
    else:
        found_engines, found_exes = _find_md_engines(
            system, protocol, engine, gpu_support=True
        )
        exe = found_exes[0]

    logging.info(f"using {exe} as an engine for the process...")

    # Create the process passing a working directory.
    process = func(system, protocol, work_dir=work_dir, exe=exe)

    # Start the process, wait for it to exit
    process.start()
    process.wait()

    # Check for errors.
    if process.isError():
        logging.error(process.stdout())
        logging.error(process.stderr())
        raise _Exceptions.ThirdPartyError("The process exited with an error!")

    # If it worked, try to get the system. No need to block, since it's already finished.
    system = process.getSystem()

    return system


def min_prots(lig_fep: str) -> (BSS.Protocol, BSS.Protocol):
    """define the minimisation protocols for the equilibration. For fepprep is at lambda 0.5

    Args:
        lig_fep (str): 'ligprep' or 'fepprep' ; which the protocols are for

    Returns:
        BSS.Protocol: the minimisation restrained and minimisation no restraints protocols.
    """

    prot_dict = {
        "ligprep": (BSS.Protocol.Minimisation, {}),
        "fepprep": (
            BSS.Protocol.FreeEnergyMinimisation,
            {"lam_vals": [0.5], "lam": 0.5},
        ),
    }

    func = prot_dict[lig_fep][0]
    args = prot_dict[lig_fep][1]

    # minimisation of solvent
    protocol_min_rest = func(steps=10000, restraint="heavy", **args)
    # minimisation full system
    protocol_min = func(steps=10000, **args)

    return protocol_min_rest, protocol_min


def nvt_prots(lig_fep: str) -> (BSS.Protocol, BSS.Protocol, BSS.Protocol):
    """define the nvt protocols for the equilibration. For fepprep is at lambda 0.5

    Args:
        lig_fep (str): 'ligprep' or 'fepprep' ; which the protocols are for

    Returns:
        BSS.Protocol: the nvt restrained all, restrained heavy, no restraints protocols.
    """

    prot_dict = {
        "ligprep": (BSS.Protocol.Equilibration, {}),
        "fepprep": (
            BSS.Protocol.FreeEnergyEquilibration,
            {"lam_vals": [0.5], "lam": 0.5},
        ),
    }

    temperature = 300

    func = prot_dict[lig_fep][0]
    args = prot_dict[lig_fep][1]

    # NVTs
    # NVT restraining all non solvent atoms
    protocol_nvt_sol = func(
        runtime=400 * BSS.Units.Time.picosecond,
        temperature_start=0 * BSS.Units.Temperature.kelvin,
        temperature_end=temperature * BSS.Units.Temperature.kelvin,
        restraint="all",
        timestep=1 * BSS.Units.Time.femtosecond,
        force_constant=100,
        restart=False,
        **args,
    )

    # NVT restraining all backbone/heavy atoms
    protocol_nvt_heavy = func(
        runtime=200 * BSS.Units.Time.picosecond,
        temperature=temperature * BSS.Units.Temperature.kelvin,
        restraint="heavy",
        force_constant=25,
        #    restart=True,
        **args,
    )
    # NVT no restraints
    protocol_nvt = func(
        runtime=200 * BSS.Units.Time.picosecond,
        # temperature_start=0 * BSS.Units.Temperature.kelvin,
        temperature=temperature * BSS.Units.Temperature.kelvin,  # temperature
        #  restart=True,
        **args,
    )

    return protocol_nvt_sol, protocol_nvt_heavy, protocol_nvt


def npt_prots(lig_fep: str) -> (BSS.Protocol, BSS.Protocol, BSS.Protocol):
    """define the npt protocols for the equilibration. For fepprep is at lambda 0.5

    Args:
        lig_fep (str): 'ligprep' or 'fepprep' ; which the protocols are for

    Returns:
        BSS.Protocol: the npt restrained heavy 10, restrained heavy 5, no restraints protocols.
    """

    prot_dict = {
        "ligprep": (BSS.Protocol.Equilibration, {}),
        "fepprep": (
            BSS.Protocol.FreeEnergyEquilibration,
            {"lam_vals": [0.5], "lam": 0.5},
        ),
    }

    temperature = 300

    func = prot_dict[lig_fep][0]
    args = prot_dict[lig_fep][1]

    # NPTs
    # NPT restraining all non solvent atoms
    protocol_npt_heavy = func(
        runtime=400 * BSS.Units.Time.picosecond,
        pressure=1 * BSS.Units.Pressure.atm,
        temperature=temperature * BSS.Units.Temperature.kelvin,
        restraint="all",
        force_constant=50,
        restart=False,
        **args,
    )
    # NPT with gradual release of restraints
    protocol_npt_heavy_lighter = func(
        runtime=400 * BSS.Units.Time.picosecond,
        pressure=1 * BSS.Units.Pressure.atm,
        temperature=temperature * BSS.Units.Temperature.kelvin,
        restraint="heavy",
        force_constant=10,
        restart=False,
        **args,
    )
    # NPT no restraints
    protocol_npt = func(
        runtime=400 * BSS.Units.Time.picosecond,
        pressure=1 * BSS.Units.Pressure.atm,
        temperature=temperature * BSS.Units.Temperature.kelvin,
        restart=False,
        **args,
    )

    return protocol_npt_heavy, protocol_npt_heavy_lighter, protocol_npt


def md_prots(runtime: int, **kwargs) -> BSS.Protocol:
    """define the npt protocols for the MD. Equilibrate before.

    Args:
        runtime (int): runtime of the MD simulation in nanoseconds.

    Returns:
        BSS.Protocol: MD protocol.
    """

    runtime = validate.integer(runtime)

    prot_dict = {
        "ligprep": (BSS.Protocol.Equilibration, {}),
        "fepprep": (
            BSS.Protocol.FreeEnergyEquilibration,
            {"lam_vals": [0.5], "lam": 0.5},
        ),
    }

    func = BSS.Protocol.Production
    args = kwargs

    temperature = 300

    # NPT no restraints
    protocol_md = func(
        runtime=runtime * BSS.Units.Time.nanosecond,
        pressure=1 * BSS.Units.Pressure.atm,
        temperature=temperature * BSS.Units.Temperature.kelvin,
        restart=False,
        **args,
    )

    return protocol_md


def minimise_equilibrate_leg(
    system_solvated: BSS._SireWrappers.System,
    engine: str = "AMBER",
    pmemd: Optional[str] = None,
    lig_fep: str = "ligprep",
    work_dir: Optional[str] = None,
    timestep: int = 2,
) -> BSS._SireWrappers.System:  # times at the top
    """minimse and equilibrate for the given leg of the pipeline

    Args:
        system_solvated (BioSimSpace._SireWrappers.System): the solvated BSS system
        engine (str, optional): engine used for the process. Defaults to "AMBER".
        pmemd (str, optional): pmemd path for use with AMBER. Defaults to None.
        lig_fep (str, optional):'ligprep' or 'fepprep'. Defaults to "ligprep".
        work_dir (str, optional): location for the runs. Defaults to None.

    Returns:
        BioSimSpace._SireWrappers.System: the final minimised and equilibrated system
    """

    if work_dir:
        work_dir = validate.folder_path(work_dir, create=True)

    # define all the protocols
    logging.info("defining min and eq protocols...")
    # Minimisations
    protocol_min_rest, protocol_min = min_prots(lig_fep)

    # NVTs
    protocol_nvt_sol, protocol_nvt_heavy, protocol_nvt = nvt_prots(lig_fep)

    # NPTs
    protocol_npt_heavy, protocol_npt_heavy_lighter, protocol_npt = npt_prots(lig_fep)

    protocol_nvt_sol.setTimeStep(timestep * BSS.Units.Time.femtosecond)
    protocol_nvt_heavy.setTimeStep(timestep * BSS.Units.Time.femtosecond)
    protocol_nvt.setTimeStep(timestep * BSS.Units.Time.femtosecond)
    protocol_npt_heavy.setTimeStep(timestep * BSS.Units.Time.femtosecond)
    protocol_npt_heavy_lighter.setTimeStep(timestep * BSS.Units.Time.femtosecond)
    protocol_npt.setTimeStep(timestep * BSS.Units.Time.femtosecond)

    # run all the protocols
    if lig_fep == "ligprep":
        logging.info("minimising...")
        minimised1 = run_process(
            system_solvated,
            protocol_min_rest,
            engine,
            pmemd,
            work_dir=f"{work_dir}/min1",
        )
        minimised2 = run_process(
            minimised1, protocol_min, engine, pmemd, work_dir=f"{work_dir}/min2"
        )
        logging.info("equilibrating NVT...")
        equil1 = run_process(
            minimised2, protocol_nvt_sol, engine, pmemd, work_dir=f"{work_dir}/nvt1"
        )

        # logging.info("nvt2")
        # equil2 = run_process(
        #     equil1, protocol_nvt_heavy, engine, pmemd, work_dir=f"{work_dir}/nvt2"
        # )
        # logging.info("nvt3")
        # equil3 = run_process(equil2, protocol_nvt, engine, pmemd, work_dir=f"{work_dir}/nvt3")

        logging.info("equilibrating NPT...")
        equil4 = run_process(
            equil1, protocol_npt_heavy, engine, pmemd, work_dir=f"{work_dir}/npt1"
        )
        equil5 = run_process(
            equil4,
            protocol_npt_heavy_lighter,
            engine,
            pmemd,
            work_dir=f"{work_dir}/npt2",
        )
        sys_equil_fin = run_process(
            equil5, protocol_npt, engine, pmemd, work_dir=f"{work_dir}/npt3"
        )

    if lig_fep == "fepprep":
        logging.info("minimising...")
        minimised1 = run_process(system_solvated, protocol_min_rest, engine, pmemd)
        minimised2 = run_process(minimised1, protocol_min, engine, pmemd)
        logging.info("equilibrating NVT...")
        equil_nvt = run_process(minimised2, protocol_nvt, engine, pmemd)
        logging.info("equilibrating NPT...")
        sys_equil_fin = run_process(equil_nvt, protocol_npt, engine, pmemd)

    return sys_equil_fin


def md_leg(
    system_equil: BSS._SireWrappers.System,
    engine: str = "AMBER",
    pmemd: Optional[str] = None,
    runtime: int = 10,
    work_dir: Optional[str] = None,
    timestep: int = 2,
) -> BSS._SireWrappers.System:  # times at the top
    """minimse and equilibrate for the given leg of the pipeline

    Args:
        system_equil (BioSimSpace._SireWrappers.System): the equilibrated BSS system .
        engine (str, optional): engine used for the process. Defaults to "AMBER".
        pmemd (str, optional): pmemd path for use with AMBER. Defaults to None.
        runtime (int, optional): Runtime of the MD simulation. Defaults to 10.
        work_dir (str, optional): location for the runs. Defaults to None.

    Returns:
        BioSimSpace._SireWrappers.System: the final system
    """

    if work_dir:
        work_dir = validate.folder_path(work_dir, create=True)

    # define all the protocols

    # MD
    protocol_md = md_prots(runtime)
    protocol_md.setTimeStep(timestep * BSS.Units.Time.femtosecond)

    sys_md_fin = run_process(
        system_equil, protocol_md, engine, pmemd, work_dir=work_dir
    )

    return sys_md_fin
