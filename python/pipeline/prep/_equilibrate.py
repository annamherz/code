import BioSimSpace as BSS
from BioSimSpace import _Exceptions
from BioSimSpace.Units.Length import angstrom as _angstrom

from ..utils._validate import *


def run_process(system, protocol, engine="AMBER", pmemd_path=None, work_dir=None):
    """    Given a solvated system (BSS object) and BSS protocol, run a process workflow with either 
    AMBER or GROMACS. Returns the processed system.

    Args:
        system (BioSimSpace._SireWrappers.System): a BSS system
        protocol (pipeline.prep.pipeline_protocol): a pipeline protocol
        engine (str, optional): engine used for the process. Defaults to "AMBER".
        pmemd_path (str, optional): path to pmemd for amber. Defaults to None.
        work_dir (str, optional): work dir for the BSS process. Defaults to None.

    Raises:
        _Exceptions.ThirdPartyError: If the process is not able to run.

    Returns:
        BioSimSpace._SireWrappers.System: the system after the process is completed.
    """

    # TODO fix so always runs pmemd

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

def min_prots(lig_fep):
    """define the minimisation protocols for the equilibration. For fepprep is at lambda 0.5

    Args:
        lig_fep (str): 'ligprep' or 'fepprep' ; which the protocols are for

    Returns:
        BSS.Protocol: the minimisation restrained and minimisation no restraints protocols.
    """

    prot_dict = {"ligprep": (BSS.Protocol.Minimisation, {}),
                 "fepprep": (BSS.Protocol.FreeEnergyMinimisation, {"lam_vals":[0.5], "lam":0.5})}
    
    func = prot_dict[lig_fep][0]
    args = prot_dict[lig_fep][1]

    # minimisation of solvent
    protocol_min_rest = func(
        steps=10000,
        restraint="all",
        **args
    )
    # minimisation full system
    protocol_min = func(
        steps=10000,
        **args
    )

    return protocol_min_rest, protocol_min

def nvt_prots(lig_fep):
    """define the nvt protocols for the equilibration. For fepprep is at lambda 0.5 

    Args:
        lig_fep (str): 'ligprep' or 'fepprep' ; which the protocols are for

    Returns:
        BSS.Protocol: the nvt restrained all, restrained heavy, no restraints protocols.
    """

    prot_dict = {"ligprep": (BSS.Protocol.Equilibration, {}),
                 "fepprep": (BSS.Protocol.FreeEnergyEquilibration, {"lam_vals":[0.5], "lam":0.5})}
    
    func = prot_dict[lig_fep][0]
    args = prot_dict[lig_fep][1]

    # NVTs
    # NVT restraining all non solvent atoms
    protocol_nvt_sol = func(
        runtime=400*BSS.Units.Time.picosecond,
        temperature=300*BSS.Units.Temperature.kelvin,
        restraint="all",
        restart=True,
        **args
    )

    # NVT restraining all backbone/heavy atoms
    protocol_nvt_heavy = func(
        runtime=400*BSS.Units.Time.picosecond,
        temperature=300*BSS.Units.Temperature.kelvin,
        restraint="heavy",
        restart=True,
        **args
    )
    # NVT no restraints
    protocol_nvt = func(
        runtime=400*BSS.Units.Time.picosecond,
        temperature_start=0*BSS.Units.Temperature.kelvin,
        temperature_end=300*BSS.Units.Temperature.kelvin,
        restart=True,
        **args
    )

    return protocol_nvt_sol, protocol_nvt_heavy, protocol_nvt

def npt_prots(lig_fep):
    """define the npt protocols for the equilibration. For fepprep is at lambda 0.5 

    Args:
        lig_fep (str): 'ligprep' or 'fepprep' ; which the protocols are for

    Returns:
        BSS.Protocol: the npt restrained heavy 10, restrained heavy 5, no restraints protocols.
    """

    prot_dict = {"ligprep": (BSS.Protocol.Equilibration, {}),
                 "fepprep": (BSS.Protocol.FreeEnergyEquilibration, {"lam_vals":[0.5], "lam":0.5})}
    
    func = prot_dict[lig_fep][0]
    args = prot_dict[lig_fep][1]

    # NPTs
    # NPT restraining all non solvent heavy atoms
    protocol_npt_heavy = func(
        runtime=400*BSS.Units.Time.picosecond,
        pressure=1*BSS.Units.Pressure.atm,
        temperature=300*BSS.Units.Temperature.kelvin,
        restraint="heavy",
        restart=True,
        **args
    )
    # NPT with gradual release of restraints
    protocol_npt_heavy_lighter = func(
        runtime=400*BSS.Units.Time.picosecond,
        pressure=1*BSS.Units.Pressure.atm,
        temperature=300*BSS.Units.Temperature.kelvin,
        restraint="heavy",
        force_constant=5,
        restart=True,
        **args
    )
    # NPT no restraints
    protocol_npt = func(
        runtime=1000*BSS.Units.Time.picosecond,
        pressure=1*BSS.Units.Pressure.atm,
        temperature=300*BSS.Units.Temperature.kelvin,
        restart=True,
        **args
    )

    return protocol_npt_heavy, protocol_npt_heavy_lighter, protocol_npt

def minimise_equilibrate_leg(system_solvated, engine="AMBER", pmemd=None, lig_fep="ligprep", work_dir=None):  # times at the top
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
    print("defining min and eq protocols...")
    # Minimisations
    protocol_min_rest, protocol_min = min_prots(lig_fep)

    # NVTs
    protocol_nvt_sol, protocol_nvt_heavy, protocol_nvt = nvt_prots(lig_fep)

    # NPTs
    protocol_npt_heavy, protocol_npt_heavy_lighter, protocol_npt = npt_prots(lig_fep)

    # run all the protocols
    if lig_fep == "ligprep":
        print("minimising...")
        minimised1 = run_process(system_solvated, protocol_min_rest, engine, pmemd, work_dir=work_dir)
        minimised2 = run_process(minimised1, protocol_min, engine, pmemd, work_dir=work_dir)
        print("equilibrating NVT...")
        equil1 = run_process(minimised2, protocol_nvt_sol, engine, pmemd, work_dir=work_dir)
        equil2 = run_process(equil1, protocol_nvt_heavy, engine, pmemd, work_dir=work_dir)
        equil3 = run_process(equil2, protocol_nvt, engine, pmemd, work_dir=work_dir)
        print("equilibrating NPT...") 
        equil4 = run_process(equil3, protocol_npt_heavy, engine, pmemd, work_dir=work_dir)
        equil5 = run_process(equil4, protocol_npt_heavy_lighter, engine, pmemd, work_dir=work_dir)
        sys_equil_fin = run_process(equil5, protocol_npt, engine, pmemd, work_dir=work_dir)
    
    if lig_fep == "fepprep":
        print("minimising...")
        minimised1 = run_process(system_solvated, protocol_min_rest, engine, pmemd)
        minimised2 = run_process(minimised1, protocol_min, engine, pmemd)
        print("equilibrating NVT...")
        equil_nvt = run_process(minimised2, protocol_nvt, engine, pmemd)
        print("equilibrating NPT...") 
        sys_equil_fin = run_process(equil_nvt, protocol_npt, engine, pmemd)

    return sys_equil_fin
