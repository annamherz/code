#! /bin/bash 

# export important file locations
export CURRENTDIR="$(pwd)"
export MAINDIRECTORY="/home/anna/Documents/benchmark/XXX" # Set file path for protein
export scripts_dir="/home/anna/Documents/code/other_scripts" # choose location of scripts

# export all execution model files for later scripts
export lig_file="$MAINDIRECTORY/execution_model/ligands.dat"
export net_file="$MAINDIRECTORY/execution_model/network.dat"
export prot_file="$MAINDIRECTORY/execution_model/protocol.dat"

# make sure engines etc are sourced correctly
export PYTHONPATH=export PYTHONPATH="/home/anna/BioSimSpace/python:$PYTHONPATH" # if using a cloned git branch of BSS - otherwise comment out
export BSS="/home/anna/anaconda3/bin/activate biosimspace-dev" # to use the conda env to make sure sire works correctly - sourced in each sh script

# sourcing - as needed in the othe sh scripts
source $BSS
source $scripts_dir/extract_execution_model_bash.sh

for i in "${!trans_array[@]}"; do
python $scripts_dir/get_atom_mappings.py ${trans_array[i]} ${eng_array[i]}
done