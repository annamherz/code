import BioSimSpace as BSS
import os
import itertools as it
import sys

from ..utils import *
from ._network import *
from ._plotting import *
from ._dictionaries import *

import cinnabar
from cinnabar import wrangle,plotting
class analysis_engines():
    """class to analyse results files and plot
    """

    def __init__(self, results_directory, exp_file=None, engines=None, net_file=None, res_folder=None, file_ext=None, extra_options=None):
        
        # TODO some way to consider different extensions for analysis methods
        # use same style of file ext as written by the analysis, maybe using same extra options dict format?

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
        self._results_repeat_files = analysis_engines._get_results_repeat_files(self)  
        self._results_files = analysis_engines._get_results_files(self)

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
        analysis_engines._set_network(self) # get network info

        # read the extra options
        analysis_engines._set_default_options(self)
        if extra_options:
            analysis_engines.set_options(extra_options)
        
        # as not yet computed, set this to false
        self._is_computed = False

        # set all the dicts for analysis
        # per engine dicts
        self.calc_pert_dict = {} # diff from the results repeat files, average
        self.cinnabar_calc_val_dict = {}  # from the cinnabar network analysis
        self.normalised_exper_val_dict = {} # normalised from the cinnabar network analysis
        self.cinnabar_calc_pert_dict = {} # from cinnabar network edges
        self.cinnabar_exper_pert_dict = {} # from cinnabar network edges

        # solo dicts for exper
        self.exper_val_dict = None # yml converted into experimental values, actual, for ligands in object
        self.exper_pert_dict = None # yml converted into experimental values, actual, for perturbations in object

        # storing the nx digraphs, per engine
        self._cinnabar_networks = {}
        # overall graph
        self.network_graph = None
        # cycles
        self.cycle_dict = {}

        # for checking against free energy workflows
        self._fwf_experimental_DDGs = None
        self._fwf_computed_relative_DDGs = {}

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

        # TODO set other default self outputs to None


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

        new_exper_val_dict = make_dict._exper_from_ligands(exper_val_dict, self.ligands)

        self.exper_val_dict = new_exper_val_dict

        return new_exper_val_dict

    def get_experimental_pert(self, exper_val_dict=None):

        if exper_val_dict:
            exper_val_dict = validate.dictionary(exper_val_dict)
        else:
            exper_val_dict = self.exper_val_dict
        
        pert_dict = make_dict._exper_from_perturbations(exper_val_dict, self.perturbations)

        self.exper_pert_dict = pert_dict

        return pert_dict


    def compute(self):

        # get all the dictionaries needed for plotting
        analysis_engines._compute_dicts(self)
        
        # compute the cycle closure
        self._make_graph()
        analysis_engines._compute_cycle_closure(self)

        # get statistics

        self._is_computed = True


    def _compute_dicts(self):

        # compute the experimental for perturbations
        analysis_engines.get_experimental(self) # get experimental val dict
        analysis_engines.get_experimental_pert(self) # from cinnabar expeirmental diff ? make_dict class

        # get the files into cinnabar format for analysis
        for eng in self.engines:
            results_files = self._results_repeat_files[eng]
            convert.cinnabar_file(results_files, self.exper_val_dict, f"{self.results_folder}/cinnabar_{eng}", perturbations=self.perturbations)
            # TODO some way to incl extension for the files in the naming here, or alternatively own folder is good
        
            # compute the per ligand for the network
            network = wrangle.FEMap(f"{self.results_folder}/cinnabar_{eng}.csv")
            self._cinnabar_networks.update({eng:network})

            # for self plotting of per ligand
            self.cinnabar_calc_val_dict.update({eng: make_dict.from_cinnabar_network_node(network, "calc")})
            self.normalised_exper_val_dict.update({eng: make_dict.from_cinnabar_network_node(network, "exp", normalise=True)})

            # for self plotting of per pert
            calc_diff_dict = make_dict.comp_results(self._results_repeat_files[eng], self.perturbations, eng) # older method
            self.calc_pert_dict.update({eng:calc_diff_dict})
            # from cinnabar graph
            self.cinnabar_calc_pert_dict.update({eng: make_dict.from_cinnabar_network_edges(network, "calc", self.perturbations)})
            self.cinnabar_exper_pert_dict.update({eng: make_dict.from_cinnabar_network_edges(network, "exp", self.perturbations)})


    def _make_graph(self):
        
        graph = net_graph(self.ligands, self.perturbations)
        self.network_graph = graph
   

    def draw_graph(self, output_dir=None, use_cinnabar=False, engine=None):
        
        if use_cinnabar:
            if engine:
                engines = [engine]
            else:
                engines = self.engines

            for eng in self.engines:
                if output_dir:
                    file_name = f"{output_dir}/cinnabar_network_{eng}.png"
                else:
                    file_name = None
                self._cinnabar_networks[eng].draw_graph(file_name=file_name)

        else:
            if not self.network_graph:
                self._make_graph()

            self.network_graph.draw_graph(file_dir=output_dir)


        # TODO also incl cinnabar graph drawing functionality?
        # TODO eng specific graph drawing or this doesnt matter


    def _compute_cycle_closure(self):

        # cycle closures

        if not self.network_graph:
            self._make_graph()
        
        network_graph = self.network_graph

        for eng in self.engines:

            pert_dict = self.calc_pert_dict[eng] # TODO can also use pert dict from cinnabar?

            cycle_closures = network_graph.cycle_closures()

            cycles = make_dict.cycle_closures(pert_dict, cycle_closures)    

            # print(f"{eng} cycle vals is {cycles[1]}")
            # print(f"{eng} cycle mean is {cycles[2]}")
            # print(f"{eng} cycle deviation is {cycles[3]}")

            self.cycle_dict.update({eng:(cycles[0], cycles[1], cycles[2], cycles[3])}) # the cycles dict

    def _compute_statistics(self):
        pass
    # def statistical_analysis()
    # # check if can extract stats from the cinnabar stuff

    def get_stats_cinnabar(self):
        pass

        # do for all networks in dict ie engines, also save



    # for all have engine options if only want to plot a single engine
    def plot_bar_pert(self, engine=None):

        if not engine:
            engines = self.engines
        else:
            engines = [validate.engine(engine)]

        # plot one can also just use the cinnabar plotting if just one

    def plot_bar_lig(self, engines=None):
        
        # TODO add file path, if cinnabar or other (if just one engine)
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




    # freenergworkflows stuff for comparison

    def _get_exp_fwf(self, fwf_path=None):
        # using freenergworkflows
        if not fwf_path:
            raise ValueError("pls incl the path to freenergworkflows")
        sys.path.insert(1, fwf_path)
        import experiments

        # first need to convert the yml file into one useable by freenergworkflows
        exp_file_dat = f"{self.exp_file.split('.')[0]}_exp_dat.dat"
        convert.yml_into_freenrgworkflows(self.exp_file, exp_file_dat)
        experiments = experiments.ExperimentalData()
        experiments.compute_affinities(exp_file_dat, data_type="IC50", comments="#", delimiter=",")
        experimental_DDGs = experiments.freeEnergiesInKcal

        exp_pert_dict,exp_lig_dict = make_dict.experimental_from_freenrgworkflows(experimental_DDGs, self.ligands, self.perturbations)
        self._fwf_experimental_DDGs = experimental_DDGs

        return exp_lig_dict, exp_pert_dict


    def _get_ana_fwf(self, fwf_path=None, engine=None):
        # using freenergworkflows 
        if not fwf_path:
            raise ValueError("pls incl the path to freenergworkflows")
        sys.path.insert(1, fwf_path)
        import networkanalysis

        if not engine:
            raise ValueError("please incl an engine")

        # using the network analyser
        nA = networkanalysis.NetworkAnalyser()

        first_file = False
        for file_name in self._results_repeat_files[engine]:
            if first_file is False:
                nA.read_perturbations_pandas(file_name, comments='#')
                first_file = True
            else:
                # add more replicates to the graph. FreeNrgWorkflows will take care of averaging 
                # the free energies as well as propagating the error.
                nA.add_data_to_graph_pandas(file_name)

        computed_relative_DDGs = nA.freeEnergyInKcal

        freenrg_dict = make_dict.from_freenrgworkflows_network_analyser(computed_relative_DDGs)
        self._fwf_computed_relative_DDGs.update({engine : computed_relative_DDGs})

        return freenrg_dict


    def _get_stats_fwf(self, fwf_path=None, engine=None):
        # using freenergworkflows 
        if not fwf_path:
            raise ValueError("pls incl the path to freenergworkflows")
        sys.path.insert(1, fwf_path)
        import stats

        if not engine:
            raise ValueError("please incl an engine")

        computed_relative_DDGs = self._fwf_computed_relative_DDGs[engine]
        experimental_DDGs = self._fwf_experimental_DDGs

        _stats = stats.freeEnergyStats()
        _stats.generate_statistics(computed_relative_DDGs,experimental_DDGs,repeats=10000)
        r_confidence = _stats.R_confidence
        tau_confidence = _stats.tau_confidence
        mue_confidence = _stats.mue_confidence
        print ("R confidence is:   %.2f < %.2f < %.2f" %(r_confidence[1], r_confidence[0], r_confidence[2]))
        print ("MUE confidence is: %.2f < %.2f < %.2f" %(mue_confidence[1], mue_confidence[0], mue_confidence[2]))
        print ("Tau confidence is: %.2f < %.2f < %.2f" %(tau_confidence[1], tau_confidence[0], tau_confidence[2]))

        return r_confidence, tau_confidence, mue_confidence 



# TODO new class that inherits from above, so can compare different methods



# for analysis between diff engines
#         # calc mae between
#         # stat compare convergence
#         # compare how fast reach concurrent results for length of runs?