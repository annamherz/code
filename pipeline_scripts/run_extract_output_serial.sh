#! /bin/bash 
#

# export important file locations
export CURRENTDIR="$(pwd)"
export MAINDIRECTORY="/home/anna/Documents/benchmark/tyk2_benchmark" # Set file path for protein
export scripts_dir="/home/anna/Documents/code/pipeline_scripts" # choose location of scripts

# export all execution model files for later scripts
export lig_file="/home/anna/Documents/benchmark/tyk2_benchmark/execution_model/ligands.dat"
export net_file="/home/anna/Documents/benchmark/tyk2_benchmark/execution_model/network_combined_reruns.dat"
export prot_file="/home/anna/Documents/benchmark/tyk2_benchmark/execution_model/protocol.dat"

# make sure engines etc are sourced correctly
export PYTHONPATH=export PYTHONPATH="/home/anna/BioSimSpace/python:$PYTHONPATH" # if using a cloned git branch of BSS - otherwise comment out
export BSS="/home/anna/anaconda3/bin/activate biosimspace-dev" # to use the conda env to make sure sire works correctly - sourced in each sh script
export amber="/home/anna/amber22/amber.sh" # sourced in each script
export gromacs="/usr/local/gromacs/bin/GMXRC" # sourced in each script

# sourcing - as needed in the othe sh scripts
source $BSS
source $amber
source $gromacs
source $scripts_dir/extract_execution_model_bash.sh

###########
echo "The folder for all these runs is $MAINDIRECTORY"
echo ${lig_array[@]}
echo ${trans_array[@]}
echo ${eng_array[@]}
echo ${win_array[@]}


# extraction and analysis
for i in "${!trans_array[@]}"; do
echo "running extraction for $i"
/bin/bash $scripts_dir/run_extract_output_slurm.sh ${trans_array[i]} ${eng_array[i]}
echo "analysing $i"
/bin/bash $scripts_dir/run_analysis_slurm.sh ${trans_array[i]} ${eng_array[i]}
done
