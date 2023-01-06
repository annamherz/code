import sys
import os
import shutil

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


def extract_output_single(trans_dir, extract_dir=None):
    """Extracts only the output files from the 

    Args:
        main_dir (str): Main directory, outputs directory.
    """
    # this should be the outputs directory
    trans_dir = validate.folder_path(trans_dir)

    if not extract_dir:
        extract_dir = f"{trans_dir.split('outputs')[0]}outputs_extracted{trans_dir.split('outputs')[1]}"

    # make directory
    extract_dir = validate.folder_path(extract_dir, create=True)

    _extract_output(trans_dir, extract_dir)


def extract_output_all(main_dir, extract_dir=None):
    """Extracts only the output files from the 

    Args:
        main_dir (str): Main directory, outputs directory.
    """
    # this should be the outputs directory
    main_dir = validate.folder_path(main_dir)
    if main_dir.split("/")[-1] != "outputs":
        raise ValueError(f"{main_dir} must be the outputs directory for extract all.")

    if not extract_dir:
        extract_dir = f"{main_dir}_extracted"

    # make directory
    extract_dir = validate.folder_path(extract_dir, create=True)

    _extract_output(main_dir, extract_dir)


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



# def extract_trajectory_frames(main_dir, extract_dir=None):
