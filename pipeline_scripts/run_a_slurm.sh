#!/bin/bash
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --job-name=ligprep
#SBATCH -o ../slurm_logs/slurm_run.out
#SBATCH -e ../slurm_logs/slurm_run.err

export MAINDIRECTORY="/home/anna/Documents/benchmark/tyk2_sage" # Set file path for protein
export PYTHONPATH="/home/anna/BioSimSpace/python" # if using a cloned git branch of BSS - otherwise comment out
export BSS="/home/anna/anaconda3/bin/activate biosimspace-dev" # to use the conda env to make sure sire works correctly - sourced in each sh script
export amber="/home/anna/amber22/amber.sh" # sourced in each script
export gromacs="/usr/local/gromacs/bin/GMXRC" # sourced in each script

# export all execution model files for later scripts
export lig_file="$MAINDIRECTORY/execution_model/ligands.dat"
export net_file="$MAINDIRECTORY/execution_model/network.dat"
export prot_file="$MAINDIRECTORY/execution_model/protocol.dat"

# sourcing - as needed in the othe sh scripts
source $BSS
source $amber
source $gromacs
source $scripts_dir/extract_execution_model_bash.sh

date
echo "Folder for these runs is : $MAINDIRECTORY"

echo "prep for $lig..."
python ligprep.py ejm31