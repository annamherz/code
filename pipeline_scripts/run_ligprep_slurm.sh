#!/bin/bash
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=5
#SBATCH --job-name=ligprep
#SBATCH -o ../slurm_logs/ligprep_%A_%a.out
#SBATCH -e ../slurm_logs/ligprep_%A_%a.err

# sourcing
source $BSS
source $amber
source $gromacs
source $scripts_dir/extract_execution_model_bash.sh

# cd /home/anna/BioSimSpace
# git checkout feature-amber-fep
# export BSS="/home/anna/anaconda3/bin/activate working"
# source $BSS
# cd $MAINDIRECTORY

date
start=`date +%s`

echo "Folder for these runs is : $MAINDIRECTORY"
echo "Ligands file is : $lig_file"

lig=${lig_array[$SLURM_ARRAY_TASK_ID]}

echo "prep for $lig..."
python $scripts_dir/ligprep.py $lig

end=`date +%s`
runtime=$((end-start))
echo "runtime was $runtime seconds, which is $((runtime/60)) minutes"