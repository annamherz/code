import sys
import os
import shutil

import MDAnalysis as mda
import MDAnalysis.transformations as trans

from ._validate import *

def add_header_simfile(trans_dir):
    """Adds header to simfiles if failed to generate.

    Args:
        trans_dir (str): path to the trans_dir for which headers will be added.
    """
    trans_dir = validate.folder_path(trans_dir)

    items_in_folder = [leg for leg in sorted(os.listdir(trans_dir))]
    legs = []

    for item in items_in_folder:
        if not "pickle" in item:
            if "bound" in item:
                legs.append(item)
            elif "free" in item:
                legs.append(item)
        else:
            pass

    for leg in legs:

        lambdas = []
        for dir in sorted(os.listdir(f"{trans_dir}/{leg}")):
            if "lambda" in dir:
                lambdas.append(dir.split("_")[1])

        for lam in lambdas:
            direc = f"{trans_dir}/{leg}/lambda_{lam}"

            sim_okay = False
            hash_counter = 0
            line_counter = 0
            # first check if the simfile.dat in the folder needs a new header
            # only check the first few lines
            with open(f"{direc}/simfile.dat") as sfile:
                for line in sfile.readlines():
                    line_counter += 1
                    if (line.startswith('#')):
                        hash_counter += 1
                        if hash_counter == 8:
                            sim_okay = True
                            break
                    if line_counter == 15:
                        break

            if not sim_okay:
                print(f"will write header for simfiles in {direc}")
                try:
                    if not os.path.exists(f"{direc}/old_simfile.dat"):
                        os.rename(f"{direc}/simfile.dat", f"{direc}/old_simfile.dat")

                    if not os.path.exists(f"{direc}/old_old_simfile.dat"):
                        os.rename(f"{direc}/old_simfile.dat", f"{direc}/old_old_simfile.dat")

                    file1 = open(f"{direc}/old_old_simfile.dat",'r')
                    file2 = open(f"{direc}/old_simfile.dat",'w')
                    
                    # reading each line from original
                    for line in file1.readlines():
                        if not (line.startswith('#')):
                            if not (line.startswith('Simulation')):
                                if not (line.startswith('Generating')):
                                    file2.write(line)
                    
                    # close and save the files
                    file2.close()
                    file1.close()

                except:
                    continue
                if os.path.exists(f"{direc}/simfile_header.dat"):
                    file = open(f"{direc}/simfile_header.dat", "w")
                else:
                    file = open(f"{direc}/simfile_header.dat", "a")
                # TODO fix so get info from simfile if possible eg re length
                file.write("#This file was generated to fix the header \n")
                file.write("#Using the somd command, of the molecular library Sire version <2022.3.0> \n")
                file.write("#For more information visit: https://github.com/michellab/Sire \n")
                file.write("# \n")
                file.write("#General information on simulation parameters: \n")
                file.write("#Simulation used 250000 moves, 4 cycles and 4000 ps of simulation time \n")
                file.write(f"#Generating lambda is		 {lam}\n")
                file.write("#Alchemical array is		 (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)\n")
                file.write("#Generating temperature is 	300 t\n")
                file.write("#Energy was saved every 200 steps \n")
                file.write("#\n")
                file.write("#\n")
                file.write("#   [step]      [potential kcal/mol]       [gradient kcal/mol]      [forward Metropolis]     [backward Metropolis]                   [u_kl] \n")
                file.close()

                # merge both files
                simfiles_tomerge = [f"{direc}/simfile_header.dat", f"{direc}/old_simfile.dat"]

                with open(f"{direc}/simfile.dat", "w") as outfile:

                    for names in simfiles_tomerge:
                        with open(names) as infile:
                            outfile.write(infile.read())
            
            else:
                pass


class extract():
    """ class for extracting output data
    """

    def __init__(self, trans_dir, extract_dir=None):

        self.trans_dir = validate.folder_path(trans_dir)
        # need to check that the trans directory contains outputs
        if not "outputs" in self.trans_dir:
            raise ValueError(f"{trans_dir} must be in an 'outputs' folder.")
        # set the extract directory
        if extract_dir:
            self.extract_dir = validate.folder_path(extract_dir, create=True)
        else:
            extract._set_extract_dir(self)
    

    def _set_extract_dir(self):

        extract_dir = f"{self.trans_dir.split('outputs')[0]}outputs_extracted{self.trans_dir.split('outputs')[1]}"

        # make directory
        extract_dir = validate.folder_path(extract_dir, create=True)
        print(f"using {extract_dir} for target file location as none was set.")

        self.extract_dir = extract_dir


    def extract_output(self):
        """Extracts only the output files from the 

        Args:
            main_dir (str): Main directory, outputs directory.
        """
        extract._extract_output(self.trans_dir, self.extract_dir)


    @staticmethod
    def _extract_output(folder, extract_dir):

        dir_list = [dirs[0] for dirs in os.walk(folder)]

        for dirs in dir_list:
            if not "min" in dirs:
                if not "heat" in dirs:
                    if not "eq" in dirs:
                        if "lambda" in dirs:      
                            new_dir = validate.folder_path(f"{extract_dir}{dirs.split(f'{folder}')[1]}", create=True)                  
                            if "AMBER" in dirs:
                                try:
                                    shutil.copyfile(f"{dirs}/amber.out", f"{new_dir}/amber.out")
                                except:
                                    print(f"{dirs} does not have a recognised AMBER input.")
                            elif "GROMACS" in dirs:
                                try:
                                    shutil.copyfile(f"{dirs}/gromacs.xvg", f"{new_dir}/gromacs.xvg")                        
                                except:
                                    print(f"{dirs} does not have a recognised GROMACS input.")
                            elif "SOMD" in dirs:
                                try:
                                    shutil.copyfile(f"{dirs}/simfile.dat", f"{new_dir}/simfile.dat")
                                except:
                                    print(f"{dirs} does not have a recognised SOMD input.")
                            
                            # if cant find the engine in the file path try all
                            else:
                                try:
                                    shutil.copyfile(f"{dirs}/amber.out", f"{new_dir}/amber.out")
                                except:
                                    print(f"{dirs} does not have a recognised AMBER input.")
                                try:
                                    shutil.copyfile(f"{dirs}/gromacs.xvg", f"{new_dir}/gromacs.xvg")                        
                                except:
                                    print(f"{dirs} does not have a recognised GROMACS input.")
                                try:
                                    shutil.copyfile(f"{dirs}/simfile.dat", f"{new_dir}/simfile.dat")
                                except:
                                    print(f"{dirs} does not have a recognised SOMD input.")
        
        
    def extract_frames(self, traj_lambdas=None, rmsd=True, overwrite=False):

        trans_dir = self.trans_dir
        rmsd = validate.boolean(rmsd)
        overwrite = validate.boolean(overwrite)

        # get all lambda windows
        items_in_folder = [leg for leg in sorted(os.listdir(trans_dir))]
        legs = []

        if traj_lambdas:
            traj_lambdas = validate.is_list(traj_lambdas)

        for item in items_in_folder:
            if not "pickle" in item:
                if "bound" in item:
                    legs.append(item)
                elif "free" in item:
                    legs.append(item)
            else:
                pass

        for leg in legs:

            lambdas = []
            for dir in sorted(os.listdir(f"{trans_dir}/{leg}")):
                if "lambda" in dir:
                    lambdas.append(dir.split("_")[1])

            for lam in lambdas:
                # only do it for considered lambdas
                if lam in traj_lambdas:

                    # get target directory for extraction and create
                    direc = f"{trans_dir}/{leg}/lambda_{lam}"
                    traj_extract_dir = f"{direc.replace(self.trans_dir, self.extract_dir)}"
                    traj_extract_dir = validate.folder_path(traj_extract_dir, create=True)

                    # check if the output files already exist in the target directory
                    if not overwrite:
                        file_exists = False
                        onlyfiles = [f for f in os.listdir(traj_extract_dir) if os.path.isfile(f"{traj_extract_dir}/{f}")]
                        if "system_0.pdb" in onlyfiles:
                            file_exists = True
                        if rmsd:
                            if "rmsd_ligand.png" in onlyfiles:
                                file_exists = True
                            else:
                                file_exists = False

                        if file_exists:
                            print("traj extracted already exists, will NOT extract again.")
                            return
                    
                    # create mda universe based on file type
                    if "SOMD" in direc:
                        # must be parm7 for mda
                        shutil.copy(f"{direc}/somd.prm7",f"{direc}/somd.parm7")
                        coord_file = validate.file_path(f"{direc}/somd.parm7")

                        # combine traj files
                        traj_files = []
                        for file in os.listdir(direc):
                            if "dcd" in file:
                                traj_files.append(validate.file_path(f"{direc}/{file}"))
                        if not traj_files:
                            raise ValueError("there are no dcd trajectory files for somd in this folder.")

                    elif "AMBER" in direc:
                        # must be parm7 for mda
                        shutil.copy(f"{direc}/amber.prm7",f"{direc}/amber.parm7")
                        coord_file = validate.file_path(f"{direc}/amber.parm7")

                        traj_files = validate.file_path(f"{direc}/amber.nc")

                    elif "GROMACS" in direc:
                        coord_file = validate.file_path(f"{direc}/gromacs.tpr")

                        traj_files = validate.file_path(f"{direc}/gromacs.trr")

                    u = mda.Universe(coord_file, traj_files)
                    if "free" in leg:
                        u = extract.centre_molecule(u, "resname LIG")
                    elif "bound" in leg:
                        u = extract.centre_molecule(u, "protein")

                    # want to write the frames for the 
                    extract._write_traj_frames(u, traj_extract_dir)
                    if rmsd:
                        extract._rmsd_trajectory(u, traj_extract_dir)
            
                else:
                    pass

    
    @staticmethod
    def centre_molecule(universe, selection):
        
        u = universe

        # want to centre the protein/ligand
        molecules = u.select_atoms(f'{selection}')
        not_molecules = u.select_atoms(f'not {selection}')
        transforms = [trans.unwrap(molecules),
                    trans.center_in_box(molecules, wrap=True),
                    trans.wrap(not_molecules)]
        u.trajectory.add_transformations(*transforms)

        return u


    @staticmethod
    def _write_traj_frames(universe, traj_extract_dir):

        write_dir = traj_extract_dir
        u = universe

        print("starting to write trajectory frames...")
        timesteps = [u.trajectory.ts for ts in u.trajectory]
        timestep_freq = int(len(timesteps)/4)
        timesteps_used = [0, timestep_freq-1, (timestep_freq*2)-1, (timestep_freq*3)-1, int(len(timesteps)-1)]

        for time in timesteps_used:
            u.trajectory[time]
            with mda.Writer(f"{write_dir}/system_{str(time)}.pdb") as W:
                W.write(u)
        
        print("finished writing 5 frames to pdb.")


    @staticmethod
    def _rmsd_trajectory(universe, traj_extract_dir):

        u = universe
        ref = universe
        write_dir = traj_extract_dir

        import MDAnalysis.analysis.rms
        R = MDAnalysis.analysis.rms.RMSD(u, ref,
                select="resname LIG",
                groupselections=["resname LIG and resname LIG",
                                    "backbone and resname LIG"])
        R.run()

        import matplotlib.pyplot as plt
        rmsd = R.rmsd.T   # transpose makes it easier for plotting
        time = rmsd[1]
        fig = plt.figure(figsize=(6,6))
        ax = fig.add_subplot(111)
        ax.plot(time, rmsd[2], 'k-', color="teal")
        ax.legend(loc="best")
        ax.set_xlabel("time (ps)")
        ax.set_ylabel(r"RMSD of ligand ($\AA$)")
        fig.savefig(f"{write_dir}/rmsd_ligand.png")  

        print(f"the maximum rmsd of the ligand is {max(rmsd[2])} and the minimum is {min(rmsd[2])}")


    @staticmethod
    def extract_output_all(main_dir, extract_dir=None):
        """Extracts only the output files if an output directory is given.

        Args:
            main_dir (str): must be the outputs directory.
        """
        # this should be the outputs directory
        main_dir = validate.folder_path(main_dir)
        if main_dir.split("/")[-1] != "outputs":
            raise ValueError(f"{main_dir} must be the outputs directory for extract all.")

        if not extract_dir:
            extract_dir = f"{main_dir}_extracted"

        # make directory
        extract_dir = validate.folder_path(extract_dir, create=True)

        extract._extract_output(main_dir, extract_dir)

        # TODO so does entire extraction for all things