# running the pipeline

`python initialise_pipeline.py -lf path_to_ligands_folder -pf path_to_protein_file -mf path_to_main_folder_to_setup`

alternatively:

`python initialise_pipeline.py`
and it will also ask for the above folder/file paths.

Additional input that will be asked for during setup:
Sampling time (per window) - Q can be excluded?
MD engines to use
Whether HMR should be applied - Q can be excluded?
How many repeats to run - Q can be excluded?
Whether to save the trajectories

Then, the network can be edited as follows:
Remove any ligands to generate the network for
Add perturbations
Remove perturbations
Whether to also run in the reverse direction.

Finally, it will ask for a path to a sources file that includes all the module exports/pythonpaths.
eg:
```
export BSS="/home/anna/anaconda3/bin/activate /home/anna/mambaforge/envs/openbiosim-dev"
export amber="/home/anna/amber20"
export gromacs="/usr/local/gromacs-23.1/bin/GMXRC"

# sourcing
source $BSS
source $amber/amber.sh
source $gromacs

export PYTHONPATH="$HOME/BioSimSpace/python:$PYTHONPATH"
export PYTHONPATH="$HOME/Documents/code/python:$PYTHONPATH"

```
can also for example include lines such as:
```
module load amber/22
```
for the slurm scheduler.

The setup can also be carried out using the setup_notebook.ipynb in this folder. This provides more extensive protocol editing and visualisation options.

In the main folder, the following folders will be generated:
execution_model (contains: network.dat, ligands.dat, protocol.dat, analysis_protocol.dat)
scripts (contains: all scripts for runnning pipeline)
run_all_slurm.sh

The run_* files in the scripts folder need to be adjusted for the slurm scheduler. - Q incl this in the setup ? can change for branch in 

running the run_all_slurm.sh should start all the runs. - Q keep analysis and extract seperately? as not using GPUs

The results can then be computed using the analyse_pipeline.py, but it is better to use the notebook (analysis_notebook.ipynb). - Q I have used branch of cinnabar, but can also find the actual version?