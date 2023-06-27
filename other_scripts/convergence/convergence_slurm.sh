#!/bin/bash
#SBATCH -n 1
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-task=30
#SBATCH --job-name=convergence
#SBATCH --time=48:00:00
#SBATCH -o convergence_%A.out
#SBATCH -e convergence_%A.err

export PYTHONPATH=export PYTHONPATH="/home/anna/BioSimSpace/python:$PYTHONPATH" 
export BSS="/home/anna/anaconda3/bin/activate biosimspace-dev"
source $BSS

date
start=`date +%s`

python truncation_check.py -mf /home/anna/Documents/benchmark/extracted/tyk2 -n /home/anna/Documents/benchmark/extracted/tyk2/execution_model/network_lomap.dat -e MBAR
python truncation_check.py -mf /home/anna/Documents/benchmark/extracted/mcl1 -n /home/anna/Documents/benchmark/extracted/mcl1/execution_model/network_lomap.dat -e MBAR
python truncation_check.py -mf /home/anna/Documents/benchmark/extracted/p38 -n /home/anna/Documents/benchmark/extracted/p38/execution_model/network_lomap.dat -e MBAR

python truncation_check.py -mf /home/anna/Documents/benchmark/extracted/tyk2 -n /home/anna/Documents/benchmark/extracted/tyk2/execution_model/network_lomap.dat -e MBAR -eq False -s True
python truncation_check.py -mf /home/anna/Documents/benchmark/extracted/mcl1 -n /home/anna/Documents/benchmark/extracted/mcl1/execution_model/network_lomap.dat -e MBAR -eq False -s True
python truncation_check.py -mf /home/anna/Documents/benchmark/extracted/p38 -n /home/anna/Documents/benchmark/extracted/p38/execution_model/network_lomap.dat -e MBAR -eq False -s True

python truncation_check.py -mf /home/anna/Documents/benchmark/extracted/tyk2 -n /home/anna/Documents/benchmark/extracted/tyk2/execution_model/network_lomap.dat -e MBAR -eq True -s True
python truncation_check.py -mf /home/anna/Documents/benchmark/extracted/mcl1 -n /home/anna/Documents/benchmark/extracted/mcl1/execution_model/network_lomap.dat -e MBAR -eq True -s True
python truncation_check.py -mf /home/anna/Documents/benchmark/extracted/p38 -n /home/anna/Documents/benchmark/extracted/p38/execution_model/network_lomap.dat -e MBAR -eq True -s True

end=`date +%s`
runtime=$((end-start))
echo "runtime was $runtime seconds, which is $((runtime/60)) minutes"
