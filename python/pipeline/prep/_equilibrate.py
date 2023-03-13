import BioSimSpace as BSS
from BioSimSpace import _Exceptions
from BioSimSpace.Units.Length import angstrom as _angstrom

from ..utils._validate import *


def run_process(system, protocol, engine="AMBER", pmemd_path=None, work_dir=None):
    """
    Given a solvated system (BSS object) and BSS protocol, run a process workflow with either 
    AMBER or GROMACS. Returns the processed system.
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

    prot_dict = {"ligprep": (BSS.Protocol.Equilibration, {}),
                 "fepprep": (BSS.Protocol.FreeEnergyEquilibration, {"lam_vals":[0.5], "lam":0.5})}
    
    func = prot_dict[lig_fep][0]
    args = prot_dict[lig_fep][1]

    # NVTs
    # NVT restraining all non solvent atoms
    protocol_nvt_sol = func(
        runtime=400*BSS.Units.Time.picosecond,
        temperature_start=0*BSS.Units.Temperature.kelvin,
        temperature_end=300*BSS.Units.Temperature.kelvin,
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
        temperature_end=300*BSS.Units.Temperature.kelvin,
        restart=True,
        **args
    )

    return protocol_nvt_sol, protocol_nvt_heavy, protocol_nvt

def npt_prots(lig_fep):

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

def minimise_equilibrate_leg(system_solvated, engine="AMBER", pmemd=None, lig_fep="ligprep"):  # times at the top
    """
    Default protocols and running of the lig_paramateriserep.
    system_solvated is a BSS system
    lig_sys is if lig or sys .
    """

    # define all the protocols
    print("defining min and eq protocols...")
    # Minimisations
    protocol_min_rest, protocol_min = min_prots(lig_fep)

    # NVTs
    protocol_nvt_sol, protocol_nvt_heavy, protocol_nvt = nvt_prots(lig_fep)

    # NPTs
    protocol_npt_heavy, protocol_npt_heavy_lighter, protocol_npt = npt_prots(lig_fep)

    # run all the protocols
    print("minimising...")
    minimised1 = run_process(system_solvated, protocol_min_rest, engine, pmemd)
    minimised2 = run_process(minimised1, protocol_min, engine, pmemd)
    print("equilibrating NVT...")
    equil1 = run_process(minimised2, protocol_nvt_sol, engine, pmemd)
    equil2 = run_process(equil1, protocol_nvt_heavy, engine, pmemd)
    equil3 = run_process(equil2, protocol_nvt, engine, pmemd)
    print("equilibrating NPT...") 
    equil4 = run_process(equil3, protocol_npt_heavy, engine, pmemd)
    equil5 = run_process(equil4, protocol_npt_heavy_lighter, engine, pmemd)
    sys_equil_fin = run_process(equil5, protocol_npt, engine, pmemd)

    return sys_equil_fin
