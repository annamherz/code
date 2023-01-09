#!/bin/bash
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=10
#SBATCH --job-name=solv
#SBATCH -o slurm_log.out
#SBATCH -e slurm_log.err

# sourcing
export PYTHONPATH="/home/anna/BioSimSpace/python"
export BSS="/home/anna/anaconda3/bin/activate biosimspace-dev"
export amber="/home/anna/amber22"
export gromacs="/usr/local/gromacs/bin/GMXRC"

source $BSS
source $amber
source $gromacs

date

python solvation_script.py
