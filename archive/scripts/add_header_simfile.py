
import sys
import os
import shutil

# folder = f"{sys.argv[1]}"
trans_dir = f"{sys.argv[1]}"
# trans_dir = f"{os.getcwd()}/{folder}"

legs = [leg for leg in sorted(os.listdir(trans_dir))]

for leg in legs:
    if leg == "graphs":
        legs.remove(leg)

lambdas = ["0.0000","0.1000","0.2000","0.3000","0.4000","0.5000","0.6000","0.7000","0.8000","0.9000","1.0000"]

for leg in legs:
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
