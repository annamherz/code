#!/bin/bash
#SBATCH -n 1
#SBATCH -p GTX980,GTX1080,RTX3080
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=extr
#SBATCH --time=01:00:00
#SBATCH -o ../slurm_logs/extr_%A.out
#SBATCH -e ../slurm_logs/extr_%A.err

# sourcing
module load cuda/11.6
module load amber/22
module load gromacs/22.2
source $scripts_dir/extract_execution_model_bash.sh
export BSS="/home/anna/anaconda3/bin/activate pipeline"
source $BSS

date
start=`date +%s`

folder=$MAINDIRECTORY/outputs/$2/$1

echo "extracting output for $folder"
python $scripts_dir/extract_output.py -f $folder -p $prot_file

echo "done."

end=`date +%s`
runtime=$((end-start))
echo "runtime was $runtime seconds, which is $((runtime/60)) minutes"