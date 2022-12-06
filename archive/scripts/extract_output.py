# import libraries

import sys
import os
import shutil

# this should be the outputs directory
main_dir = f"{sys.argv[1]}"
extract_dir = f"{main_dir}_extracted"

# make directory
if not os.path.isdir(extract_dir):
    os.mkdir(extract_dir)

dir_list = [dirs[0] for dirs in os.walk(main_dir)]

for dirs in dir_list:
    if not "min" in dirs:
        if not "heat" in dirs:
            if not "eq" in dirs:
                if not os.path.isdir(f"{extract_dir}{dirs.split(f'{main_dir}')[1]}"):
                    os.mkdir(f"{extract_dir}{dirs.split(f'{main_dir}')[1]}")

                if "lambda" in dirs:                        
                    if "AMBER" in dirs:
                        try:
                            shutil.copyfile(f"{dirs}/amber.out", f"{extract_dir}/{dirs.split(f'{main_dir}')[1]}/amber.out")
                        except:
                            print(f"{dirs} does not have a recognised AMBER input.")
                    elif "GROMACS" in dirs:
                        try:
                            shutil.copyfile(f"{dirs}/gromacs.xvg", f"{extract_dir}/{dirs.split(f'{main_dir}')[1]}/gromacs.xvg")                        
                        except:
                            print(f"{dirs} does not have a recognised GROMACS input.")
                    elif "SOMD" in dirs:
                        try:
                            shutil.copyfile(f"{dirs}/simfile.dat", f"{extract_dir}/{dirs.split(f'{main_dir}')[1]}/simfile.dat")
                        except:
                            print(f"{dirs} does not have a recognised SOMD input.")
                        