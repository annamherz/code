#!/bin/bash 
#
#run all the lig prep, FEP prep, and the production runs
# Prior to this, the execution model should have been set up using the overall_network.ipynb
# This should have also tested if all the required parts are installed.

# export important file locations
export CURRENTDIR="$(pwd)"
export MAINDIRECTORY="/home/anna/Documents/benchmark/XXX" # Set file path for protein
export scripts_dir="/home/anna/Documents/code/pipeline_scripts" # choose location of scripts
export protein_file="/home/anna/Documents/benchmark/inputs/XXX_parameterised" # this should be the prm7 and rst7 file name. best as input folder 
export ligands_folder="/home/anna/Documents/benchmark/inputs/XXX"

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
export PYTHONPATH=export PYTHONPATH="/home/anna/BioSimSpace/python:$PYTHONPATH" # if using a cloned git branch of BSS - otherwise comment out
export BSS="/home/anna/anaconda3/bin/activate biosimspace-dev" # to use the conda env to make sure sire works correctly - sourced in each sh script
export amber="/home/anna/amber22" # sourced in each script
export gromacs="/usr/local/gromacs/bin/GMXRC" # sourced in each script

# sourcing - as needed in the othe sh scripts
source $BSS
source $amber/amber.sh
source $gromacs
source $scripts_dir/extract_execution_model_bash.sh

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

# chmod all files so can be executed by sbatch.
# chmod u+x run_ligprep_slurm.sh
# chmod u+x run_fepprep_slurm.sh
# chmod u+x run_production_slurm.sh
# chmod u+x run_analysis_slurm.sh

# Run the runs
# ligand prep
jidlig=$(sbatch --parsable --array=0-$((${#lig_array[@]}-1)) $scripts_dir/run_ligprep_slurm.sh)
echo "ligand prep jobid is $jidlig"

# FEP prep
jidfep=$(sbatch --dependency=afterany:${jidlig} --parsable --array=0-$((${#trans_array[@]}-1)) $scripts_dir/run_fepprep_slurm.sh)
echo "FEP prep jobid is $jidfep"

# Production runs and analysis for the transformation
for i in "${!trans_array[@]}"; do
jidprod=$(sbatch --dependency=afterany:${jidfep} --parsable --array=0-$((${win_array[i]}-1)) $scripts_dir/run_production_slurm.sh ${trans_array[i]} ${eng_array[i]} ${win_array[i]} $repeats)
echo "Production jobid for ${trans_array[i]}, ${eng_array[i]} is $jidprod"
jidextract=$(sbatch --dependency=afterany:${jidprod} --parsable $scripts_dir/run_extract_output_slurm.sh ${trans_array[i]} ${eng_array[i]})
echo "Extraction jobid for ${trans_array[i]}, ${eng_array[i]} is $jidextract"
jidana=$(sbatch --dependency=afterany:${jidextract} --parsable $scripts_dir/run_analysis_slurm.sh ${trans_array[i]} ${eng_array[i]})
echo "Analysis jobid for ${trans_array[i]}, ${eng_array[i]} is $jidana"
done

# python $scripts_dir/analysis_network.py
