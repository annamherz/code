import BioSimSpace as BSS
import os
import itertools as it
import sys
import re

from ..utils import *
from ._network import *
from ._analysis import *
from ._plotting import *
from ._statistics import *
from ._convergence import *
from ._dictionaries import *

import cinnabar
from cinnabar import wrangle,plotting,stats

class analysis_network():
    """class to analyse results files and plot
    """

    def __init__(self, results_directory, exp_file=None, engines=None, net_file=None, output_folder=None, analysis_ext=None, extra_options=None):
        
        self._results_directory = validate.folder_path(results_directory)

        # get engines for analysis
        self.engines = validate.engines(engines)
        
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
        
        if not output_folder:
            print("no output folder provided, writing all output to the 'results_directory'.")
            self.output_folder = validate.folder_path(f"{results_directory}")
            self.graph_dir = validate.folder_path(f"{results_directory}/graphs", create=True)
        else:
            self.output_folder = validate.folder_path(output_folder, create=True)
            self.graph_dir = validate.folder_path(f"{output_folder}/graphs", create=True)

        if not analysis_ext:
            self.file_ext = ".+" # wildcard, all files in the folder included
        else:
            try:
                # see if it is just a file extension name
                self.file_ext = validate.string(file_ext)
            except:
                # read in a dictionary or analysis protocol file
                try:
                    analysis_options = analysis_protocol(analysis_ext, auto_validate=True)
                    self.file_ext = analyse.file_ext(analysis_options)
                except:
                    raise TypeError(f"{analysis_ext} analysis ext must be either a string or an analysis protocol file/dictionary")

        
        # TODO add a file ext extra option for when writing out the files re naming eg so can incl diff networks
        # Actually add network info

        # set defaults
        self.temperature = 300

        # overwrite if in extra options
        if extra_options:
            extra_options = validate.dictionary(extra_options)

            # temperature for converting experimental values
            if "temperature" in extra_options.keys():
                self.temperature = [validate.is_float(extra_options["temperature"])]
                

        # get files from results directory
        self._results_repeat_files = self._get_results_repeat_files()  
        self._results_free_repeat_files = self._get_results_repeat_files(leg="free")  
        self._results_bound_repeat_files = self._get_results_repeat_files(leg="bound")  
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
        # per engine dicts (also used for other results)
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

        # for other results
        self.other_results_names = []

        # for checking against free energy workflows
        self._fwf_experimental_DDGs = None
        self._fwf_computed_relative_DDGs = {}

        # for plotting
        self._plotting_object = None
        # for stats
        self._stats_object = None

    def _get_results_repeat_files(self, leg=None):
        res_dir = self._results_directory
        all_files = os.listdir(res_dir)

        files_for_dict = []

        if leg: # leg should be free or bound
            for file in all_files:
                if f"{leg}_" in file:
                    if re.search(self.file_ext, file):
                        files_for_dict.append(f"{res_dir}/{file}")

        else: # search for the freenrg
            for file in all_files:
                if "freenrg" in file:
                    files_for_dict.append(f"{res_dir}/{file}")

        files_dict = {}
        
        for eng in self.engines:
            eng_files = []
            for file in files_for_dict:
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
        file_name = [res_file for sublist in results_lists for res_file in sublist]

        values = get_info_network(results_files=file_name,
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
            if not self.exp_file:
                warnings.warn("need an experimental file to proceed with most of the calculations. please set using self.get_experimental(file)")
            else:
                exp_file = self.exp_file
        else:
            self.exp_file = validate.file_path(exp_file)

        # check if yml, as this is needed for 
        if exp_file.split(".")[-1] != "yml":
            raise ValueError("the provided experimental file should be in yml format")
        
        # TODO more checks for type of data
        exper_val_dict = convert.yml_into_exper_dict(exp_file, temp=self.temperature) # this output is in kcal/mol

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

    def remove_perturbations(self, pert):
        
        perts = validate.is_list(pert, make_list=True)
        for pert in perts:
            self.perturbations.remove(pert)
        # remove plotting object as needs to be reintialised with new perturbations
        self._plotting_object = None
        self._stats_object = None

    def compute(self, cycle_closure=True, statistics=True):

        # get all the dictionaries needed for plotting
        self._compute_dicts()
        
        # TODO have these entirely seperate?
        # compute the cycle closure
        if cycle_closure:
            self._make_graph()
            self._compute_cycle_closures()

        if statistics:
            # TODO get statistics
            self.pert_statistics = {}
            self.val_statistics = {}
            # for eng in self.engines:
                # self.pert_statistics.update({eng: self.compute_statistics(pert_val="pert", engine=eng)})
                # self.val_statistics.update({eng: self.compute_statistics(pert_val="val", engine=eng)})
        
        # initialise plotting and stats objects
        self._initalise_plotting_object()
        self._initalise_stats_object()

        self._is_computed = True


    def _compute_dicts(self):

        # compute the experimental for perturbations
        self.get_experimental() # get experimental val dict and normalised dict
        self.get_experimental_pert() # from cinnabar expeirmental diff ? make_dict class

        # get the files into cinnabar format for analysis
        for eng in self.engines:
            file_name = self._results_repeat_files[eng]
            cinnabar_file_name = f"{self.output_folder}/cinnabar_{eng}_{self.file_ext}_{self.net_ext}"
            convert.cinnabar_file(file_name, self.exper_val_dict, cinnabar_file_name, perturbations=self.perturbations)

            try:
                # compute the per ligand for the network
                network = wrangle.FEMap(f"{cinnabar_file_name}.csv")
                self._cinnabar_networks.update({eng:network})
                
                # for self plotting of per ligand
                self.cinnabar_calc_val_dict.update({eng: make_dict.from_cinnabar_network_node(network, "calc")})
                self.cinnabar_exper_val_dict.update({eng: make_dict.from_cinnabar_network_node(network, "exp", normalise=True)})

                # from cinnabar graph
                self.cinnabar_calc_pert_dict.update({eng: make_dict.from_cinnabar_network_edges(network, "calc", self.perturbations)})
                self.cinnabar_exper_pert_dict.update({eng: make_dict.from_cinnabar_network_edges(network, "exp", self.perturbations)})

            except Exception as e:
                print(e)
                print(f"could not create cinnabar network for {eng}")

            # for self plotting of per pert
            calc_diff_dict = make_dict.comp_results(self._results_repeat_files[eng], self.perturbations, eng) # older method
            self.calc_pert_dict.update({eng:calc_diff_dict})


    def compute_other_results(self, file_name=None, name=None):
        
        if not isinstance(file_name, list):
            file_name = [file_name]
        for file in file_name:
            validate.file_path(file)
        name = validate.string(name)
        self.other_results_names.append(name)
        new_file_path = f"{file_name[0].replace(file_name[0].split('/')[-1], '')[:-1]}/{name}_results_file"
        # TODO write list of results files if multiple are passed in

        results_dict = make_dict.comp_results(file_name, perturbations=None, engine=None, output_file=new_file_path)
        perturbations,ligands = get_info_network_from_dict(results_dict)

        # first convert the other results into a value dict, do this using cinnabar but no experimental
        cinnabar_file_name = f"{self.output_folder}/cinnabar_{name}_{self.file_ext}_{self.net_ext}"
        convert.cinnabar_file([f"{new_file_path}.csv"], self.exper_val_dict, cinnabar_file_name, perturbations=self.perturbations)
    
        # compute the per ligand for the network
        network = wrangle.FEMap(f"{cinnabar_file_name}.csv")
        self._cinnabar_networks.update({name:network})

        self.cinnabar_calc_val_dict.update({name: make_dict.from_cinnabar_network_node(network, "calc")})
        self.cinnabar_calc_pert_dict.update({name: make_dict.from_cinnabar_network_edges(network, "calc", perturbations)})

        # initialise plotting and stats objects again so theyre added
        self._initalise_plotting_object()
        self._initalise_stats_object()

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

            pert_dict = self.cinnabar_calc_pert_dict[eng]

            cycle_closures = network_graph.cycle_closures()

            cycles = make_dict.cycle_closures(pert_dict, cycle_closures)    

            # print(f"{eng} cycle vals is {cycles[1]}")
            # print(f"{eng} cycle mean is {cycles[2]}")
            # print(f"{eng} cycle deviation is {cycles[3]}")

            self.cycle_dict.update({eng:(cycles[0], cycles[1], cycles[2], cycles[3])}) # the cycles dict

    def cycle_closures(self):
        # TODO 
        # if is computed
        # take data from cycle closure dict and plot
        pass

    def _initalise_plotting_object(self, check=False):

        # if not checking, always make
        if not check:
            self._plotting_object = plotting_engines(analysis_object=self)

        # if checking, first see if it exists and if not make
        elif check:
            if not self._plotting_object:
                self._plotting_object = plotting_engines(analysis_object=self)
        
        return self._plotting_object

    def plot_bar_pert(self, engine=None):

        plot_obj = self._initalise_plotting_object(check=True)
        plot_obj.bar(pert_val="pert", engines=engine)

    def plot_bar_lig(self, engine=None):

        plot_obj = self._initalise_plotting_object(check=True)
        plot_obj.bar(pert_val="val", engines=engine)


    def plot_scatter_pert(self, engine=None, use_cinnabar=False):
        
        if use_cinnabar:
            if engine:
                try:
                    engine = validate.engine(engine)
                    engines = [engine]
                except:
                    print("for cinnabar plotting, can only have one engine. Please use the engine keyword to define, or leave blank to autoidentify available engines.")
                    return
            
            else:
                print("no engine specified for plotting, will plot seperate graphs for each self.engines ")
                engines = self.engines
            
            for eng in engines:
                plotting.plot_DDGs(self._cinnabar_networks[eng].graph,
                                filename=f"{self.graph_dir}/DDGs_{eng}_{self.file_ext}_{self.net_ext}.png",
                                title=f"DDGs for {eng}, {self.net_ext}") #with {self.file_ext}

        else:
            plot_obj = self._initalise_plotting_object(check=True)
            plot_obj.scatter(pert_val="pert", engines=engine)

    def plot_scatter_lig(self, engine=None, use_cinnabar=False):
        
        if use_cinnabar:
            if engine:
                try:
                    engine = validate.engine(engine)
                    engines = [engine]
                except:
                    print("for cinnabar plotting, can only have one engine. Please use the engine keyword to define, or leave blank to autoidentify available engines.")
                    return
            
            else:
                print("no engine specified for plotting, will plot seperate graphs for each self.engines ")
                engines = self.engines
            
            for eng in engines:
                plotting.plot_DGs(self._cinnabar_networks[eng].graph,
                                filename=f"{self.graph_dir}/DGs_{eng}_{self.file_ext}_{self.net_ext}.png",
                                title=f"DGs for {eng}, {self.net_ext}") #with {self.file_ext}

        else:
            plot_obj = self._initalise_plotting_object(check=True)
            plot_obj.scatter(pert_val="val", engines=engine)

    def plot_eng_vs_eng(self, engine_a=None, engine_b=None, pert_val="pert"):
        
        plot_obj = self._initalise_plotting_object(check=True)

        if pert_val == "pert":
            binding = "$\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$"
        elif pert_val == "val":
            binding = "$\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$"
        plotting_dict = {"title":f"{engine_a} vs {engine_b}\n for {self.file_ext}, {self.net_ext}",
                            "y label":f"{engine_a} "+ binding,
                            "x label":f"{engine_b} " + binding,
                            "key": False
                            }
        
        plot_obj.scatter(pert_val=pert_val, engines=engine_a, name=engine_b, **plotting_dict)


    def plot_other_results(self, name=None, engine=None, pert_val=None, outliers=None):

        name = validate.string(name)
        
        # first, check if name exists in other dict
        if name in self.other_results_names:
            pass
        else:
            raise NameError(f"{name} does not exist as an added other results.")

        plot_obj = self._initalise_plotting_object(check=True)

        if outliers:
            plot_obj.outlier(pert_val=pert_val, engines=engine, outliers=outliers, name=name)  
        else:
            plot_obj.scatter(pert_val=pert_val, engines=engine, name=name)


    def plot_outliers(self, engine=None, outliers=5, pert_val="pert"):

        plot_obj = self._initalise_plotting_object(check=True)
        plot_obj.outlier(pert_val=pert_val, engines=engine, outliers=outliers)  

    
    def plot_histogram_runs(self, engine=None):
        
        plot_obj = self._initalise_plotting_object(check=True)

        # get list of the errors from the repeat files
        error_dict = {}

        for eng in self.engines:
            error_list = make_dict.error_list_from_files(self._results_repeat_files[eng])
            error_dict.update({eng:error_list})

        plot_obj.histogram(engines=engine, error_dict=error_dict)

    def plot_histogram_sem(self, engine=None, pert_val="val"):

        plot_obj = self._initalise_plotting_object(check=True)
        plot_obj.histogram(engines=engine, pert_val=pert_val)
        
    # TODO histogram plot for the error in the free and the error in the bound legs

    def plot_histogram_legs(self, engine=None, free_bound=None):

        free_bound = validate.string(free_bound)
        if free_bound.lower() not in ["free","bound"]:
            raise ValueError("free_bound must be either 'free' or 'bound' to denote one of the legs of the run.")
        
        plot_obj = self._initalise_plotting_object(check=True)

        error_dict = {}
        
        for eng in self.engines:
            
            if free_bound == "free":
                error_list = make_dict.error_list_from_files(self._results_free_repeat_files[eng])
            if free_bound == "bound":
                error_list = make_dict.error_list_from_files(self._results_bound_repeat_files[eng])

            error_dict.update({eng:error_list})

        plot_obj.histogram(engines=engine, error_dict=error_dict, file_ext=free_bound)

    def _initalise_stats_object(self, check=False):

        # if not checking, always make
        if not check:
            self._stats_object = stats_engines(analysis_object=self)

        # if checking, first see if it exists and if not make
        elif check:
            if not self._stats_object:
                self._stats_object = stats_engines(analysis_object=self)
        
        return self._stats_object

    def calc_mae(self, pert_val=None):

        stats_obj = self._initalise_stats_object(check=True)

        pv = validate.pert_val(pert_val)

        values_dict = stats_obj.values_dict
        engines = self.engines

        mae_pert_df = pd.DataFrame(columns=engines,index=engines)
        mae_pert_df_err = pd.DataFrame(columns=engines,index=engines)

        # iterate over all possible combinations
        for combo in it.product(engines, engines):
            eng1 = combo[0]
            eng2 = combo[1]

            values = stats_obj.compute_mue(pv, x=eng1, y=eng2)
            mean_absolute_error = values[0] # the computed statistic    
            mae_err = values[1] # the stderr from bootstrapping

            mae_pert_df.loc[eng1,eng2]=mean_absolute_error
            mae_pert_df_err.loc[eng1,eng2]=mae_err

        mae_pert_df.to_csv(f"{self.output_folder}/mae_pert_{self.file_ext}.csv", sep=" ")
        mae_pert_df_err.to_csv(f"{self.output_folder}/mae_pert_err_{self.file_ext}.csv", sep=" ")

        return mae_pert_df, mae_pert_df_err

    def compute_mue(self, pert_val=None, engines=None):
        
        stats_obj = self._initalise_stats_object(check=True)
        pert_val = validate.pert_val(pert_val)

        if not engines:
            engines = self.engines
        else:
            engines = validate.engines(engines)

        for eng in engines:
            values = stats_obj.compute_mue(pert_val=pert_val, y=eng) # TODO fix

        return values

    # TODO compute all stats function
    # TODO write an output file for the dictionary               

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