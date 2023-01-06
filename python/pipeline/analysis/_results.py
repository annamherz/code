from ..utils import *
from ._network import *

class analysis():
    """class to analyse results files and plot
    """

    def __init__(self, results_files, exp_file, engine=None, net_file=None, res_folder=None, extra_options=None):

        # validate the input
        if not engine:
            raise ValueError("Must supply an engine")
        else:
            self.engine = validate.engine(engine)
        
        self._results_files = []
        validate.is_list(results_files)
        for file in results_files:
            validate.file_path(file)
        self._exp_file = validate.file_path(exp_file)
        # TODO this should be a yml file

        if not net_file:
            raise ValueError("need a network file currently")
        else:
            self._net_file = validate.file_path(net_file)
        
        if not res_folder:
            print("no results folder provided, writing all output to 'results' in the current directory.")
            self.results_folder = validate.file_path("results")
        else:
            self.results_folder = validate.file_path(res_folder)

        # create all the dictionaries needed for analysis
        self._results_files = None
        self.perturbations = None
        self.ligands = None  
        self._values_dict = {}
        self._values_dict["pert_results"] = None
        self._values_dict["val_results"] = None
        self._values_dict["freenrg_df_pert"] = None
        self._values_dict["freenrg_df_val"] = None
        self._values_dict["freenrgworkflows_ouput"] = None
        self._values_dict["cycle_closures"] = None

        # get info from the network
        analysis._set_network(self)

        # read the extra options
        # create a series of default options for the analysis
        # eg if want to use freenergworkflows
        # setting for each if do freenergworkflows and if do cinnabar     

        # convert into the necessary format


        # make a dict of the computed results
        comp_dict = make_dict.comp_results(results_files=self._values_dict["results_files"], 
                                            perturbations=self._values_dict["perts"],
                                            output_file=f"{res_folder}/{comp_pert_file_name}_{eng}",
                                            engine=eng)
        self._values_dict["pert_results"] = comp_dict

    
    def _set_network(self):

        # first get info from the network
        values = get_info_network(results_files=self._results_files,
                                net_file=self._net_file,
                                output_folder=self.results_folder,
                                extra_options = {"engine":self.engine}
                                )

        self._results_files = values[2]
        self.perturbations = values[0]
        self.ligands = values[1]   


        # incooperate cinnabar. get stats from here?
        # diff parts here for if plotting pert or per ligand
        # dict to df for plotting, check if df in variables if not remake
        # plot eg from one bar
        # one scatter
        # plot one can also just use the cinnabar plotting if just one
        # possible to change the colour of the cinnabar plotting for me here?
        # check if can extract stats from the cinnabar stuff
        # change these plots so also take engines as list, but if length of list is one just plot one eg ehere
        # do convergence 

class multiple():
    """ class for comparing multiple engines
    """

    def __init__(self, engine_dict):
        # dict of engines with analysis object for each engine
    

        # collate ligands and perturbations list for all
        all_perturbations = []
        all_ligands = []
        for pert in values[0]:
            if pert not in all_perturbations:
                all_perturbations.append(pert)
        for lig in values[1]:
            if lig not in all_ligands:
                all_ligands.append(lig)
        

        # calc mae between
        # stat compare convergence
        # compare how fast reach concurrent results for length of runs?