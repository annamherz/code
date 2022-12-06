#!/bin/bash
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=5
#SBATCH --job-name=fepprep
#SBATCH -o ../slurm_logs/fepprep_%A_%a.out
#SBATCH -e ../slurm_logs/fepprep_%A_%a.err

# sourcing
source $BSS
source $amber
source $gromacs
source $scripts_dir/extract_execution_model_bash.sh

# cd /home/anna/BioSimSpace
# git checkout feature-amber-pre-2023
# export BSS="/home/anna/anaconda3/bin/activate biosimspace-dev"
# source $BSS
# cd $MAINDIRECTORY

date
start=`date +%s`

echo "Folder for these runs is : $MAINDIRECTORY"
echo "Network file is : $net_file"

trans=${trans_array[$SLURM_ARRAY_TASK_ID]}
eng=${eng_array[$SLURM_ARRAY_TASK_ID]}
win=${win_array[$SLURM_ARRAY_TASK_ID]}

echo "fepprep for $trans using $eng"
python $scripts_dir/fepprep.py $trans $eng $win # no eq during setup

end=`date +%s`
runtime=$((end-start))
echo "runtime was $runtime seconds, which is $((runtime/60)) minutes"