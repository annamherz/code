#!/bin/bash
#SBATCH -n 1
#SBATCH -p GTX980,GTX1080,RTX3080
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=ana
#SBATCH --time=01:00:00
#SBATCH -o ../slurm_logs/ana_%A.out
#SBATCH -e ../slurm_logs/ana_%A.err

# sourcing
module load cuda/11.6
module load amber/22
module load gromacs/22.2
source $scripts_dir/extract_execution_model_bash.sh
export BSS="/home/anna/anaconda3/bin/activate pipeline"
source $BSS

date
start=`date +%s`

echo "Folder for these runs is : $MAINDIRECTORY"
echo "Analysis for the transformation $1, $2."

python $scripts_dir/analysis.py -pert $1 -eng $2 -mf $MAINDIRECTORY -a $ana_file

end=`date +%s`
runtime=$((end-start))
echo "runtime was $runtime seconds, which is $((runtime/60)) minutes"
