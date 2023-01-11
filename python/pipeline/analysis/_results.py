import BioSimSpace as BSS
import os
import itertools as it

from ..utils import *
from ._network import *
from ._plotting import *
from ._dictionaries import *

class analysis_engine():
    """class to analyse results files and plot
    """

    def __init__(self, results_directory, exp_file=None, engines=None, net_file=None, res_folder=None, extra_options=None):
        
        # TODO some way to consider different extensions for analysis methods
        # get engines for analysis
        if not engines:
            self.engines = BSS.FreeEnergy.engines()
        else:
            try:
                try:
                    engines = validate.is_list(engines)
                    val_engines = []
                    for engine in engines:
                        engine_val = validate.engine(engine)
                        val_engines.append(engine_val)
                    self.engines = val_engines
                # if single engine string, put into list
                except:
                    engines = validate(engines)
                    self.engines = [engines]
            except:
                print("engine input not recognised. Will use all engines.")
                self.engines = BSS.FreeEnergy.engines()
        
        # get files from results directory
        self._results_directory = validate.folder_path(results_directory)
        self._results_repeat_files = analysis_engine._get_results_repeat_files(self)  
        self._results_files = analysis_engine._get_results_files(self)

        if not exp_file:
            print("please set an experimental yml file so this can be used, eg using .get_experimental(exp_file). ")
            self.exp_file = None
        else:
            self.exp_file = validate.file_path(exp_file)

        if not net_file:
            print("no network file, will use all perturbations found in results files from the results dir.")
            self._net_file = None
        else:
            self._net_file = validate.file_path(net_file)
        
        if not res_folder:
            print("no results output folder provided, writing all output to the 'results_directory'.")
            self.results_folder = validate.folder_path(f"{results_directory}")
            self.graph_dir = validate.folder_path(f"{results_directory}/graphs", create=True)
        else:
            self.results_folder = validate.folder_path(res_folder, create=True)
            self.graph_dir = validate.folder_path(f"{res_folder}/graphs", create=True)



        # get info from the network
        self.perturbations = None
        self.ligands = None  
        analysis_engine._set_network(self) # get network info

        # read the extra options
        analysis_engine._set_default_options(self)
        if extra_options:
            analysis_engine.set_options(extra_options)
        
        # as not yet computed, set this to false
        self._is_computed = False


    def _get_results_repeat_files(self):
        res_dir = self._results_directory
        all_files = os.listdir(res_dir)
        rep_files = []
        for file in all_files:
            if "repeat" in file:
                rep_files.append(f"{res_dir}/{file}")
        
        files_dict = {}
        
        for eng in self.engines:
            eng_files = []
            for file in rep_files:
                if eng in file:
                    eng_files.append(file)
            files_dict[eng] = eng_files
        
        return files_dict

    def _get_results_files(self):
        res_dir = self._results_directory
        all_files = os.listdir(res_dir)
        sum_files = []
        for file in all_files:
            if "summary" in file:
                sum_files.append(f"{res_dir}/{file}")
        
        files_dict = {}
        
        for eng in self.engines:
            eng_files = []
            for file in sum_files:
                if eng in file:
                    eng_files.append(file)
            files_dict[eng] = eng_files
        
        return files_dict

    
    def _set_network(self):

        # get perturbations and ligands for the network
        # if network file, get from that
        # if not, use all results files for this

        # get results files from dict into a list, flatten the list
        results_lists = list(self._results_files.values()) + list(self._results_repeat_files.values())
        results_files = [res_file for sublist in results_lists for res_file in sublist]

        values = get_info_network(results_files=results_files,
                                net_file=self._net_file,
                                extra_options = {"engines":self.engines}
                                )

        self.perturbations = values[0]
        self.ligands = values[1]   

    def _set_default_options(self):

        self._compute_experimental = True
        self._compute_per_ligand = True
        self._comupute_cycle_closure = True

        self._save_pickle = True
        self._try_pickle = True


    def set_options(self, options_dict):

        options_dict = validate.dictionary(options_dict)

        if "compute_experimental" in options_dict:
            self._compute_experimental = validate.boolean(options_dict["compute_experimental"])
        if "compute_per_ligand" in options_dict:
            self._compute_per_ligand = validate.boolean(options_dict["compute_per_ligand"])
        if "compute_cycle_closure" in options_dict:
            self._compute_cycle_closure = validate.boolean(options_dict["compute_cycle_closure"])
        if "XXXX" in options_dict:
            self._XXXX = validate.boolean(options_dict["XXXX"])

        if "save_pickle" in options_dict:
            self._save_pickle = validate.boolean(options_dict["save_pickle"])
        if "try_pickle" in options_dict:
            self._try_pickle = validate.boolean(options_dict["try_pickle"])


    def get_experimental(self, exp_file=None):

        if not exp_file:
            exp_file = self.exp_file
        else:
            self.exp_file = validate.file_path(exp_file)

        # check if yml, as this is needed for 
        if exp_file.split(".")[-1] != "yml":
            raise ValueError("the provided experimental file should be in yml format")
        
        # TODO more checks for type of data
        exper_val_dict = convert.yml_into_exper_dict(exp_file) # this output is in kcal/mol

        new_exper_val_dict = make_dict._from_cinnabar_experimental_val(exper_val_dict, self.ligands)

        self.exper_val_dict = new_exper_val_dict

        return new_exper_val_dict

    def get_experimental_pert(self, exper_val_dict=None):

        if exper_val_dict:
            exper_val_dict = validate.dictionary(exper_val_dict)
        else:
            exper_val_dict = self.exper_val_dict
        
        pert_dict = make_dict._from_cinnabar_experimental_diff(exper_val_dict, self.perturbations)

        self.exper_pert_dict = pert_dict

        return pert_dict

    def compute(self):

        # compute the experimental for perturbations
        analysis_engine.get_experimental(self) # get experimental val dict
        analysis_engine.get_experimental_pert(self) # from cinnabar expeirmental diff ? make_dict class

        # compute the per ligand for the network

        # compute the cycle closure

        # make a dict of the computed results
        comp_dict = make_dict.comp_results(results_files=self._values_dict["results_files"], 
                                            perturbations=self._values_dict["perts"],
                                            output_file=f"{res_folder}/{comp_pert_file_name}_{eng}",
                                            engine=eng)
        self._values_dict["pert_results"] = comp_dict

        self._is_computed = True

        # to use cinnabar, first get the results files, exper dict, and convert this using cinnabar_file method written in _convert
        # make so have cinnabar file ready to go for plotting and analysis w it


    # for all have engine options if only want to plot a single engine
    def plot_bar_pert(self, engine=None):

        if not engine:
            engines = self.engines
        else:
            engines = [validate.engine(engine)]

        # plot one can also just use the cinnabar plotting if just one

    def plot_bar_lig(self, engines=None):

        if not engine:
            engines = self.engines
        else:
            engines = [validate.engine(engine)]

    def plot_scatter_pert(self, engines=None):

        if not engine:
            engines = self.engines
        else:
            engines = [validate.engine(engine)]

    def plot_scatter_lig(self, engines=None):

        if not engine:
            engines = self.engines
        else:
            engines = [validate.engine(engine)]


    # def plot_convergence()
    # def statistical_analysis()
    # # check if can extract stats from the cinnabar stuff
    # def cycle_closure()

        
        # dict to df for plotting, check if df in variables if not remake



# for analysis between diff engines
#         # calc mae between
#         # stat compare convergence
#         # compare how fast reach concurrent results for length of runs?

        # self._values_dict = {}
        # self._values_dict["pert_results"] = None
        # self._values_dict["val_results"] = None
        # self._values_dict["freenrg_df_pert"] = None
        # self._values_dict["freenrg_df_val"] = None
        # self._values_dict["freenrgworkflows_ouput"] = None
        # self._values_dict["cycle_closures"] = None