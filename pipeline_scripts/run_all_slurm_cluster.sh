#!/bin/bash 
#
#run all the lig prep, FEP prep, and the production runs
# Prior to this, the execution model should have been set up using the overall_network.ipynb
# This should have also tested if all the required parts are installed.

# export important file locations
export CURRENTDIR="$(pwd)"
export MAINDIRECTORY="/export/users/anna/reruns/XXX" # Set file path for protein
export scripts_dir="/export/users/anna/code/pipeline_scripts" # choose location of scripts

# if have a different, already prepped ligands folder, can change this here
export prep_folder="$MAINDIRECTORY/prep"

# export all execution model files for later scripts
export lig_file="$MAINDIRECTORY/execution_model/ligands.dat"
export net_file="$MAINDIRECTORY/execution_model/network.dat"
export prot_file="$MAINDIRECTORY/execution_model/protocol.dat"
export ana_file="$MAINDIRECTORY/execution_model/analysis_protocol.dat"

# remove any ^M from end of file lines
dos2unix "$lig_file"
dos2unix "$net_file"
dos2unix "$prot_file"
dos2unix "$ana_file"

# make sure using the correct BSS branch
# cd /home/anna/BioSimSpace
# git checkout benchmark-2022
# cd $MAINDIRECTORY

# make sure engines etc are sourced correctly
export BSS="/home/anna/anaconda3/bin/activate pipeline"
source $BSS
module load cuda/11.6
module load amber/22
module load gromacs/22.2
source $scripts_dir/extract_execution_model_bash.sh

###########
echo "The folder for all these runs is $MAINDIRECTORY"
echo ${trans_array[@]}
echo ${eng_array[@]}
echo ${win_array[@]}

# make output dir for slurm out and err files
if [[ ! -d ../slurm_logs ]]; then
    mkdir ../slurm_logs
fi

# chmod all files so can be executed by sbatch.
# chmod u+x run_ligprep_slurm.sh
# chmod u+x run_fepprep_slurm.sh
# chmod u+x run_production_slurm.sh
# chmod u+x run_analysis_slurm.sh

# Production runs and analysis for the transformation
for i in "${!trans_array[@]}"; do
jidprod=$(sbatch --parsable --array=0-$((${win_array[i]}-1)) $scripts_dir/run_production_slurm_cluster.sh ${trans_array[i]} ${eng_array[i]} ${win_array[i]} $repeats)
echo "Production jobid for ${trans_array[i]}, ${eng_array[i]} is $jidprod"
jidextract=$(sbatch --dependency=afterany:${jidprod} --parsable $scripts_dir/run_extract_output_slurm_cluster.sh ${trans_array[i]} ${eng_array[i]})
echo "Extraction jobid for ${trans_array[i]}, ${eng_array[i]} is $jidextract"
jidana=$(sbatch --dependency=afterany:${jidextract} --parsable $scripts_dir/run_analysis_slurm_cluster.sh ${trans_array[i]} ${eng_array[i]})
echo "Analysis jobid for ${trans_array[i]}, ${eng_array[i]} is $jidana"
done
