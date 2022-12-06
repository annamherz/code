#! /bin/bash 
#
#run all the lig prep, FEP prep, and the production runs
# Prior to this, the execution model should have been set up using the overall_network.ipynb
# This should have also tested if all the required parts are installed.

# export important file locations
export CURRENTDIR="$(pwd)"
export MAINDIRECTORY="/export/users/anna/tyk2_benchmark" # Set file path for protein
export scripts_dir="/export/users/anna/scripts" # choose location of scripts
export protein_file="/export/users/anna/inputs/tyk2/tyk2_parameterised" # this should be the prm7 and rst7 file name. best as input folder 
export ligands_folder="/export/users/anna/inputs/tyk2/ligands"

# export all execution model files for later scripts
export lig_file="$MAINDIRECTORY/execution_model/ligands.dat"
export net_file="$MAINDIRECTORY/execution_model/network_combined.dat"
export prot_file="$MAINDIRECTORY/execution_model/protocol.dat"

###########
echo "The folder for all these runs is $MAINDIRECTORY"
echo ${lig_array[@]}
echo ${trans_array[@]}
echo ${eng_array[@]}
echo ${win_array[@]}

# make output dir for slurm out and err files
if [[ ! -d ../slurm_logs ]]; then
    mkdir ../slurm_logs
fi

# Production runs and analysis for the transformation
for i in "${!trans_array[@]}"; do
jidprod=$(sbatch --parsable --array=0-$((${win_array[i]}-1)) $scripts_dir/run_production_cluster.sh ${trans_array[i]} ${eng_array[i]} ${win_array[i]} $repeats)
echo "Production jobid for ${trans_array[i]}, ${eng_array[i]} is $jidprod"

