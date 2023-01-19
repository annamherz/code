#!/bin/bash
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=5
#SBATCH --job-name=extr
#SBATCH --time=01:00:00
#SBATCH -o ../slurm_logs/extr_%A.out
#SBATCH -e ../slurm_logs/extr_%A.err

# sourcing
source $scripts_dir/extract_execution_model_bash.sh

date
start=`date +%s`

folder=$MAINDIRECTORY/outputs/$2/$1

echo "extracting output for $folder"
python $scripts_dir/extract_output.py $folder

echo "done."

end=`date +%s`
runtime=$((end-start))
echo "runtime was $runtime seconds, which is $((runtime/60)) minutes"