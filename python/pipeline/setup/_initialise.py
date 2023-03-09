import glob
import BioSimSpace as BSS

from ..utils._validate import *
from ..utils._files import *
from ..prep._protocol import *
from ..analysis._network import *

class initialise_pipeline():

    def __init__(self):

        # need a main folder
        # execuction folder automatically made in main folder
        # protein file path

        self.ligands = None
        self.ligand_names = None
        self.ligands_dict = None
        self.perturbations = None

        self._ligands_folder = None
        self._protein_path = None

        self._main_folder = None
        self._exec_folder = None
        # TODO set ligands file and network file etc so can change name of this also change then when writing bash script

        self.protocol = None
        self.analysis_protocol = None
        
        self._is_ligands_setup = False
        self._is_network_setup = False


    def ligands_folder(self, folder_path=None):

        if folder_path:
            folder_path = validate.folder_path(folder_path)

            ligand_files = glob.glob(f"{folder_path}/*.sdf")
            if len(ligand_files) < 3:
                raise ValueError("need atleast two *.sdf ligands in folder for setting up the pipeline.")

            self._ligands_folder = folder_path
        
        else:
            pass

        return self._ligands_folder

    def main_folder(self, folder_path=None, create_exec_folder=True):

        if folder_path:
            # make main folder
            self._main_folder = validate.folder_path(folder_path, create=True)

            if create_exec_folder:
                self.exec_folder(f"{self._main_folder}/execution_model")
        
        else:
            pass

        return self._main_folder


    def exec_folder(self, folder_path=None):

        if folder_path:
            self._exec_folder = validate.folder_path(folder_path, create=True)
        else:
            pass

        return self._exec_folder

    def protein_path(self, protein_path=None):

        if protein_path:
            try:
                validate.file_path(f"{protein_path}.rst7")
                validate.file_path(f"{protein_path}.prm7")
                self._protein_path = protein_path
            except:
                raise ValueError("protein path must be that to the rst7 and prm7 parameterised protein")
        else:
            pass

        return self._protein_path

    # TODO copy protein and ligands to main folder so have inputs all together

    @staticmethod
    def _setup_ligands(path_to_ligands, exec_folder):

        #generate transformation network based on ligands
        ligand_files = sorted(glob.glob(f"{path_to_ligands}/*.sdf"))

        ligands = []
        ligand_names = []
        ligands_dict = {}

        for filepath in ligand_files:
            # append the molecule object to a list.
            lig = BSS.IO.readMolecules(filepath)[0]
            ligands.append(lig)
            
            # append the molecule name to another list so that we can use the name of each molecule in our workflow.
            lig_name = filepath.split("/")[-1].replace(".sdf","")
            ligand_names.append(lig_name)

            ligands_dict[lig_name] = lig

        print(f"there were {len(ligand_names)} ligands found. These are:\n")
        for lig in ligand_names:
            print(lig)

        # write ligands file.
        write_ligands(ligand_names, f"{exec_folder}/ligands.dat")         

        return ligands, ligand_names, ligands_dict           

    def setup_ligands(self):

        if not self._ligands_folder:
            raise ValueError("please provide a ligands folder first using .ligands_folder(path_to_ligands)")
        if not self._exec_folder:
            raise ValueError("please provide an execution model folder first using .exec_folder(path_to_exec) or create it using the .main_folder(path, create_exec_folder=True)")

        ligands, ligand_names, ligands_dict = initialise_pipeline._setup_ligands(self._ligands_folder, self._exec_folder)

        self.ligands = ligands
        self.ligand_names = ligand_names
        self.ligands_dict = ligands_dict

        self._is_ligands_setup = True

    def remove_ligand(self, lig):

        if self._is_ligands_setup:
            if lig in self.ligand_names:
                self.ligand_names.remove(lig)
                self.ligands.remove(self.ligands_dict[lig])
                del self.ligands_dict[lig]

            # write ligands file again as updted
            write_ligands(self.ligand_names, f"{self._exec_folder}/ligands.dat")  

        else:
            print("please setup ligands first before removing any")

    @staticmethod
    def _setup_network(ligands, ligand_names, folder, links_file=None):

        transformations, lomap_scores = BSS.Align.generateNetwork(ligands, plot_network=True,
                                                                  names=ligand_names,
                                                                  work_dir=f"{folder}/visualise_network",
                                                                  links_file=links_file
                                                                  )

        # dict of perts
        pert_network_dict = {}
        transformations_named = [(ligand_names[transf[0]], ligand_names[transf[1]]) for transf in transformations]
        for score, transf in sorted(zip(lomap_scores, transformations_named)):
            pert_network_dict[transf] = score
            print(transf, score)

        return pert_network_dict
            

    def setup_network(self, folder="LOMAP", links_file=None):

        if not self.ligand_names or not self.ligands:
            print("please run setup_ligands before setting up the network")
            return

        pert_network_dict = initialise_pipeline._setup_network(self.ligands,
                                           self.ligand_names,
                                           f"{self.exec_folder()}/{folder}",
                                           links_file=links_file)
        
        self.pert_network_dict = pert_network_dict
        self.perturbations = [f"{key[0]}~{key[1]}" for key in pert_network_dict.keys()]

        write_lomap_scores(pert_network_dict, f"{self.exec_folder()}/lomap_scores.dat")

        self._is_network_setup = True

    def remove_perturbation(self, pert):

        if self._is_network_setup:
            if pert in self.perturbations:
                self.perturbations.remove(pert)
                del self.pert_network_dict[(f"{pert.split('~')[0]}", f"{pert.split('~')[1]}")]

            write_lomap_scores(self.pert_network_dict, f"{self.exec_folder()}/lomap_scores.dat")

        else:
            print("please setup network first before removing any perturbations")

    def add_perturbation(self, pert):

        if self._is_network_setup:

            lig0 = f"{pert.split('~')[0]}"
            lig1 = f"{pert.split('~')[1]}"
            
            # regenerate the network just for that perturbation
            single_transformation, single_lomap_score = BSS.Align.generateNetwork([self.ligands_dict[lig0], self.ligands_dict[lig1]],
                                                                                  names=[lig0, lig1],
                                                                                  plot_network=False)

            self.pert_network_dict[(lig0, lig1)] = single_lomap_score[0]
            self.perturbations.append(pert)

            write_lomap_scores(self.pert_network_dict, f"{self.exec_folder()}/lomap_scores.dat")

        else:
            print("please setup network first before adding any perturbations")

    def draw_network(self, folder=None):

        if not folder:
            folder = self.exec_folder()
        
        if not self._is_network_setup:
            print("please setup network first")
            return
        
        graph = net_graph(self.ligand_names, self.perturbations)
        graph.draw_graph(file_dir=folder)


    def setup_protocols(self, protocol_dictionary=None):

        protocol = pipeline_protocol(protocol_dictionary, auto_validate=True)
        protocol.rewrite_protocol(file_path=f"{self.exec_folder()}/protocol.dat")
        self.protocol = protocol

        # create an analysis protocol
        ana_protocol = analysis_protocol(auto_validate=True)
        ana_protocol.rewrite_protocol(file_path=f"{self.exec_folder()}/analysis_protocol.dat")
        self.analysis_protocol = ana_protocol

    def add_pipeline_protocol(self, protocol):

        self.protocol = validate.pipeline_protocol(protocol)
        self.protocol.rewrite_protocol(file_path=f"{self.exec_folder()}/protocol.dat")

    def add_analysis_protocol(self, protocol):

        self.analysis_protocol = validate.analysis_protocol(protocol)
        self.analysis_protocol.rewrite_protocol(file_path=f"{self.exec_folder()}/analysis_protocol.dat")


    def write_network(self, file_path=None):

        if self._is_network_setup:
            if not file_path:
                file_path = f"{self.exec_folder()}/network.dat"
            write_network(self.pert_network_dict, self.protocol, file_path)

        else:
            print("please setup network first before writing the network file.")


    def write_run_all(self):

        # rewrite all the ligands, network files

        with open(f"{self._main_folder}/run_all_slurm.sh", "w") as rsh:
            rsh.write(f'''\
#!/bin/bash 
#
# generated to run all the lig prep, FEP prep, and the production runs
# please check all file paths and MD engines are sourced to correct path between the hashtags

###########

# make sure engines etc are sourced correctly
export PYTHONPATH=export PYTHONPATH="/home/anna/BioSimSpace/python:$PYTHONPATH" # if using a cloned git branch of BSS - otherwise comment out
export BSS="/home/anna/anaconda3/bin/activate biosimspace-dev" # to use the conda env to make sure sire works correctly - sourced in each sh script
export amber="/home/anna/amber22/amber.sh" # sourced in each script
export gromacs="/usr/local/gromacs/bin/GMXRC" # sourced in each script

# export important file locations
export MAINDIRECTORY="{self._main_folder}"
export protein_file="{self._protein_path}"
export ligands_folder="{self._ligands_folder}"

# export all execution model files for later scripts
export lig_file="$MAINDIRECTORY/execution_model/ligands.dat"
export net_file="$MAINDIRECTORY/execution_model/network.dat"
export prot_file="$MAINDIRECTORY/execution_model/protocol.dat"
export ana_file="$MAINDIRECTORY/execution_model/analysis_protocol.dat"

# this should be the location of in the pipeline module
export scripts_dir="/home/anna/Documents/code/pipeline_scripts" # choose location of scripts

# remove any ^M from end of file lines
dos2unix "$lig_file"
dos2unix "$net_file"
dos2unix "$prot_file"
dos2unix "$ana_file"

# chmod all files so can be executed by sbatch.
# chmod u+x $scripts_dir/run_ligprep_slurm.sh
# chmod u+x $scripts_dir/run_fepprep_slurm.sh
# chmod u+x $scripts_dir/run_production_slurm.sh
# chmod u+x $scripts_dir/run_extract_output_slurm.sh
# chmod u+x $scripts_dir/run_analysis_slurm.sh

# sourcing - as needed in the othe sh scripts
source $BSS
source $amber
source $gromacs
source $scripts_dir/extract_execution_model_bash.sh

###########

''')

            with open(f"{self._main_folder}/run_all_slurm.sh", "w") as rsh:
                rsh.write('''\

echo "The folder for all these runs is $MAINDIRECTORY"
echo ${lig_array[@]}
echo ${trans_array[@]}
echo ${eng_array[@]}
echo ${win_array[@]}

# make output dir for slurm out and err files
if [[ ! -d ../slurm_logs ]]; then
    mkdir ../slurm_logs
fi

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

''')