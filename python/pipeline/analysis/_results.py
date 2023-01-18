import BioSimSpace as BSS
import os
import itertools as it
import sys
import re

from ..utils import *
from ._network import *
from ._plotting import *
from ._dictionaries import *

import cinnabar
from cinnabar import wrangle,plotting,stats

class analysis_engines():
    """class to analyse results files and plot
    """

    def __init__(self, results_directory, exp_file=None, engines=None, net_file=None, res_folder=None, file_ext=None, extra_options=None):
        
        self._results_directory = validate.folder_path(results_directory)

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
                    engines = validate.string(engines)
                    engines = validate.engine(engines)
                    self.engines = [engines]
            except:
                print("engine input not recognised. Will use all engines.")
                self.engines = BSS.FreeEnergy.engines()
        
        if not exp_file:
            print("please set an experimental yml file so this can be used, eg using .get_experimental(exp_file). ")
            self.exp_file = None
        else:
            self.exp_file = validate.file_path(exp_file)

        if not net_file:
            print("no network file, will use all perturbations found in results files from the results dir.")
            self._net_file = None
            self.net_ext = "network"
        else:
            self._net_file = validate.file_path(net_file)
            self.net_ext = validate.string(f"{net_file.split('/')[-1].split('.')[0]}")
        
        if not res_folder:
            print("no results output folder provided, writing all output to the 'results_directory'.")
            self.results_folder = validate.folder_path(f"{results_directory}")
            self.graph_dir = validate.folder_path(f"{results_directory}/graphs", create=True)
        else:
            self.results_folder = validate.folder_path(res_folder, create=True)
            self.graph_dir = validate.folder_path(f"{res_folder}/graphs", create=True)

        if not file_ext:
            self.file_ext = ".+" # wildcard, all files in the folder included
        else:
            self.file_ext = validate.string(file_ext)
            # TODO add so can also read in dictionary like analysis??
        
        # TODO add a file ext extra option for when writing out the files re naming eg so can incl diff networks
        # Actually add network info

        # get files from results directory
        self._results_repeat_files = self._get_results_repeat_files()  
        self._results_files = self._get_results_files()

        # get info from the network
        self.perturbations = None
        self.ligands = None  
        self._set_network() # get network info

        # read the extra options
        self._set_default_options()
        if extra_options:
            self.set_options(extra_options)
        
        # as not yet computed, set this to false
        self._is_computed = False

        # set all the dicts for analysis
        # per engine dicts
        self.calc_pert_dict = {} # diff from the results repeat files, average
        self.cinnabar_calc_val_dict = {}  # from the cinnabar network analysis
        self.cinnabar_exper_val_dict = {} # normalised from the cinnabar network analysis
        self.cinnabar_calc_pert_dict = {} # from cinnabar network edges
        self.cinnabar_exper_pert_dict = {} # from cinnabar network edges

        # solo dicts for exper
        self.exper_val_dict = None # yml converted into experimental values, actual, for ligands in object
        self.normalised_exper_val_dict = None # yml converted into experimental values, then normalised
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

        # for plotting
        self._plotting_object = None

    def _get_results_repeat_files(self):
        res_dir = self._results_directory
        all_files = os.listdir(res_dir)
        rep_files = []
        for file in all_files:
            if "repeat" in file:
                if re.search(self.file_ext, file):
                    rep_files.append(f"{res_dir}/{file}")
        
        files_dict = {}
        
        for eng in self.engines:
            eng_files = []
            for file in rep_files:
                if eng in file:
                    if re.search(self.file_ext, file):
                        eng_files.append(file)
            files_dict[eng] = eng_files
        
        return files_dict

    def _get_results_files(self):
        res_dir = self._results_directory
        all_files = os.listdir(res_dir)
        sum_files = []
        for file in all_files:
            if "summary" in file:
                if re.search(self.file_ext, file):
                    sum_files.append(f"{res_dir}/{file}")
        
        files_dict = {}
        
        for eng in self.engines:
            eng_files = []
            for file in sum_files:
                if eng in file:
                    if re.search(self.file_ext, file):
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
        self._comupute_cycles = True

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
            self._compute_cycles = validate.boolean(options_dict["compute_cycles"])
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

        # experimental value dict
        new_exper_val_dict = make_dict._exper_from_ligands(exper_val_dict, self.ligands)
        self.exper_val_dict = new_exper_val_dict

        # normalise the experimental values
        normalised_exper_val_dict = make_dict._exper_from_ligands(exper_val_dict, self.ligands, normalise=True)
        self.normalised_exper_val_dict = normalised_exper_val_dict

        return new_exper_val_dict, normalised_exper_val_dict

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
        self._compute_dicts()
        
        # compute the cycle closure
        self._make_graph()
        self._compute_cycle_closures()

        # get statistics
        self.pert_statistics = {}
        self.val_statistics = {}
        # for eng in self.engines:
            # self.pert_statistics.update({eng: self.compute_statistics(pert_val="pert", engine=eng)})
            # self.val_statistics.update({eng: self.compute_statistics(pert_val="val", engine=eng)})

        self._is_computed = True


    def _compute_dicts(self):

        # compute the experimental for perturbations
        self.get_experimental() # get experimental val dict and normalised dict
        self.get_experimental_pert() # from cinnabar expeirmental diff ? make_dict class

        # get the files into cinnabar format for analysis
        for eng in self.engines:
            results_files = self._results_repeat_files[eng]
            convert.cinnabar_file(results_files, self.exper_val_dict, f"{self.results_folder}/cinnabar_{eng}_{self.file_ext}_{self.net_ext}", perturbations=self.perturbations)
        
            # compute the per ligand for the network
            network = wrangle.FEMap(f"{self.results_folder}/cinnabar_{eng}_{self.file_ext}_{self.net_ext}.csv")
            self._cinnabar_networks.update({eng:network})

            # for self plotting of per ligand
            self.cinnabar_calc_val_dict.update({eng: make_dict.from_cinnabar_network_node(network, "calc")})
            self.cinnabar_exper_val_dict.update({eng: make_dict.from_cinnabar_network_node(network, "exp", normalise=True)})

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

            for eng in engines:
                if output_dir:
                    file_name = f"{output_dir}/cinnabar_network_{eng}_{self.file_ext}_{self.net_ext}.png"
                else:
                    file_name = None
                self._cinnabar_networks[eng].draw_graph(file_name=file_name)

        else:
            if not self.network_graph:
                self._make_graph()

            self.network_graph.draw_graph(file_dir=output_dir)


        # TODO also incl cinnabar graph drawing functionality?
        # TODO eng specific graph drawing or this doesnt matter


    def _compute_cycle_closures(self):

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


    def _initalise_plotting_object(self, check=False):

        # if not checking, always make
        if not check:
            self._plotting_object = plotting_engines(analysis_object=self)

        # if checking, first see if it exists and if not make
        elif check:
            if not self._plotting_object:
                self._plotting_object = plotting_engines(analysis_object=self)


    def plot_all(self):

        self._initalise_plotting_object(check=True)
        plot_obj = self._plotting_object


    def plot_bar_pert(self, engine=None):

        self._initalise_plotting_object(check=True)
        plot_obj = self._plotting_object
        plot_obj.bar(pert_val="pert", engines=engine)

    def plot_bar_lig(self, engine=None):

        self._initalise_plotting_object(check=True)
        plot_obj = self._plotting_object
        plot_obj.bar(pert_val="val", engines=engine)


    def plot_scatter_pert(self, engine=None, use_cinnabar=False):
        
        if use_cinnabar:
            try:
                engine = validate.engine(engine)
            except:
                print("for cinnabar plotting, can only have one engine. Please use the engine keyword to define.")
                return
            
            plotting.plot_DDGs(self._cinnabar_networks[engine].graph,
                              filename=f"{self.graph_dir}/DDGs_{engine}_{self.file_ext}_{self.net_ext}.png",
                              title=f"DDGs for {engine} with {self.file_ext}, {self.net_ext}")

        else:
            self._initalise_plotting_object(check=True)
            plot_obj = self._plotting_object
            plot_obj.scatter(pert_val="pert", engines=engine)

    def plot_scatter_lig(self, engine=None, use_cinnabar=False):
        
        if use_cinnabar:
            try:
                engine = validate.engine(engine)
            except:
                print("for cinnabar plotting, can only have one engine. Please use the engine keyword to define.")
                return
            
            plotting.plot_DGs(self._cinnabar_networks[engine].graph,
                              filename=f"{self.graph_dir}/DGs_{engine}_{self.file_ext}_{self.net_ext}.png",
                              title=f"DGs for {engine} with {self.file_ext}, {self.net_ext}")

        else:
            self._initalise_plotting_object(check=True)
            plot_obj = self._plotting_object
            plot_obj.scatter(pert_val="val", engines=engine)


    def plot_outliers(self, engine=None, outliers=5, pert_val="pert"):

        self._initalise_plotting_object(check=True)
        plot_obj = self._plotting_object
        plot_obj.outlier(pert_val=pert_val, engines=engine, outliers=outliers)  

    
    def plot_histogram_pert(self, engine=None):

        self._initalise_plotting_object(check=True)
        plot_obj = self._plotting_object
        plot_obj.histogram(engines=engine, pert_val="pert")

    def plot_histogram_lig(self, engine=None):

        self._initalise_plotting_object(check=True)
        plot_obj = self._plotting_object
        plot_obj.histogram(engines=engine, pert_val="val")
    
    # TODO histogram plot for each repeat


    # def plot_convergence(self):

    #     self._initalise_plotting_object(check=True)
    #     plot_obj = self._plotting_object
        
    #     output_dir = f"{self.results_folder}/convergence_plots/{self.file_ext}"
        
    #     for pert in self.perturbations:
    #         plot_obj.outlier(perturbation=pert, engines=self.engines, output_dir=output_dir)  





# TODO make a static method?
    def compute_statistics(self, pert_val=None, engine=None, x=None, y=None, xerr=None, yerr=None):

        eng = validate.engine(engine)

        if eng:
            if pert_val == "pert":
                x = [val[0] for val in self.cinnabar_exper_pert_dict[eng]]
                y = [val[0] for val in self.cinnabar_calc_pert_dict[eng]]
                xerr = np.asarray([val[1] for val in self.cinnabar_exper_pert_dict[eng]])
                yerr = np.asarray([val[1] for val in self.cinnabar_calc_pert_dict[eng]])
            elif pert_val == "val":
                x = [val[0] for val in self.cinnabar_exper_val_dict[eng]]
                y = [val[0] for val in self.cinnabar_calc_val_dict[eng]]
                xerr = np.asarray([val[1] for val in self.cinnabar_exper_val_dict[eng]])
                yerr = np.asarray([val[1] for val in self.cinnabar_calc_val_dict[eng]])          
            else:
                raise ValueError("pert_val must be 'pert' for perturbations or 'val' for values")      
        else:
            x = x
            y = y
            xerr = xerr
            yerr = yerr

        statistics  = ["RMSE", "MUE"]
        statistics_string = ""

        for statistic in statistics:
            s = stats.bootstrap_statistic(x, y, xerr, yerr, statistic=statistic)
            string = f"{statistic}:   {s['mle']:.2f} [95%: {s['low']:.2f}, {s['high']:.2f}] " + "\n"
            statistics_string += string
            
        return statistics_string

    # # check if can extract stats from the cinnabar stuff

    def get_stats_cinnabar(self):
        pass

        # do for all networks in dict ie engines, also save


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