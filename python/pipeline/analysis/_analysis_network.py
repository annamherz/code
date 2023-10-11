import os
import itertools as it
import sys
import re

from ..utils import *
from ._network import *
from ._analysis import *
from ._plotting import *
from ._statistics import *
from ._dictionaries import *
from ._convert import *
from ..prep import *

from typing import Union, Optional

from cinnabar import wrangle, plotting, stats

from math import isnan


class analysis_network:
    """class to analyse results files and plot"""

    def __init__(
        self,
        output_folder: Optional[str] = None,
        exp_file: Optional[str] = None,
        engines: Optional[str] = None,
        net_file: Optional[str] = None,
        results_folder: Optional[str] = None,
        analysis_prot: Optional[analysis_protocol] = None,
        method: Optional[str] = None,
        extra_options: Optional[dict] = None,
    ):
        """analyses the network for a certain system

        Args:
            output_folder (str, optional): directory where all the outputs are located. Should be outputs_extracted in the main folder. This should contain a results folder where the outputs from the analysis was saved. Defaults to None.
            exp_file (str, optional): file path to the experimental results file. Defaults to None.
            engines (list, optional): engines to consider (to comapre multiple). Defaults to None.
            net_file (str, optional): file path to the network of perturbations to analyse. Defaults to None.
            results_folder (str, optional): folder path to where to save all the outputs of the analysis. Defaults to None.
            analysis_prot (pipeline.protocol.analysis_protocol, optional): analysis protocol to make file extension to look for for the results. Defaults to None.
            method (str, optional): Method to consider in the method column of the results files. Defaults to None (will consider all results in the file).
            extra_options (dict, optional): extra options (eg temperature). Defaults to None.

        Raises:
            TypeError: analysis ext must be the correct type.
        """

        if method:
            self.method = validate.string(method)
        else:
            self.method = None

        # get engines for analysis
        if engines:
            self.engines = validate.engines(engines)
        else:
            self.engines = validate.engines("ALL")

        if not exp_file:
            logging.error(
                "please set an experimental yml/csv file so this can be used, eg using .get_experimental(exp_file). "
            )
            self.exp_file = None
        else:
            self.exp_file = validate.file_path(exp_file)

        if not net_file:
            logging.info(
                "no network file, will use all perturbations found in results files from the results dir."
            )
            self._net_file = None
            self.net_ext = "network"
        else:
            self._net_file = validate.file_path(net_file)
            self.net_ext = validate.string(f"{net_file.split('/')[-1].split('.')[0]}")

        if not analysis_prot:
            self.file_ext = ".+"  # wildcard, all files in the folder included
            self.analysis_options = analysis_protocol(file=None, auto_validate=True)
        else:
            # read in a dictionary or analysis protocol file
            try:
                self.analysis_options = analysis_protocol(
                    analysis_prot, auto_validate=True
                )
                self.file_ext = analyse.file_ext(self.analysis_options)
            except:
                raise TypeError(
                    f"{analysis_prot} analysis protocol must be an analysis protocol file/dictionary"
                )

        if output_folder:
            self.output_folder = validate.folder_path(output_folder)
            # get files from results directory
            self._results_repeat_files = self._get_results_repeat_files()
            self._results_free_repeat_files = self._get_results_repeat_files(leg="free")
            self._results_bound_repeat_files = self._get_results_repeat_files(
                leg="bound"
            )
            self._results_files = self._get_results_files()
            self._results_value_files = {}
        else:
            logging.critical(
                "There is no provided results directory. There are no results to analyse. This will probably create an issue for many other functions. please reinstantiate the object with a results directory."
            )
            self.output_folder = None
            self._results_repeat_files = None
            self._results_free_repeat_files = None
            self._results_bound_repeat_files = None
            self._results_files = None
            self._results_value_files = None

        if not results_folder:
            if self.output_folder:
                logging.info(
                    "no output folder provided, writing all output to the 'output_folder/analysis'."
                )
                self.results_folder = f"{self.output_folder}/analysis"
                self.graph_dir = validate.folder_path(
                    f"{self.results_folder}/graphs", create=True
                )
            else:
                logging.info(
                    "no output or results directory, so writing files to current folder..."
                )
                self.results_folder = os.getcwd()
        else:
            self.results_folder = validate.folder_path(results_folder, create=True)
            self.graph_dir = validate.folder_path(
                f"{results_folder}/graphs", create=True
            )

        # set defaults
        self.temperature = 300

        # overwrite if in extra options
        if extra_options:
            extra_options = validate.dictionary(extra_options)

            # temperature for converting experimental values
            if "temperature" in extra_options.keys():
                self.temperature = [validate.is_float(extra_options["temperature"])]

        # get info from the network
        self.perturbations = None
        self._perturbations_dict = {}
        self.ligands = None
        self._ligands_dict = {}
        self._set_network()  # get network info

        # as not yet computed, set this to false
        self._is_computed_dicts = False

        self._set_dictionary_outputs()

    def _get_results_repeat_files(self, leg: Optional[str] = None) -> dict:
        """get the files of all the repeats for a specific leg. Used during init to set free and bound repeat files.

        Args:
            leg (str, optional): Which leg to get the repeats for, ie 'free' or 'bound'. Defaults to None.

        Returns:
            dict: dict of engines as keys and their repeat files for the defined leg.
        """
        res_dir = f"{self.output_folder}/results"
        all_files = os.listdir(res_dir)

        files_for_dict = []

        if leg:  # leg should be free or bound
            for file in all_files:
                if f"{leg}_" in file:
                    if re.search(self.file_ext, file):
                        files_for_dict.append(f"{res_dir}/{file}")

        else:  # search for the freenrg
            for file in all_files:
                if "freenrg" in file:
                    files_for_dict.append(f"{res_dir}/{file}")

        files_dict = {}

        for engine in self.engines:
            eng_files = []
            for file in files_for_dict:
                if engine in file:
                    if re.search(self.file_ext, file):
                        eng_files.append(file)
            files_dict[engine] = eng_files

        return files_dict

    def _get_results_files(self) -> dict:
        """get the summary results files

        Returns:
            dict: dict of engines as keys and their summary files as values.
        """
        res_dir = f"{self.output_folder}/results"
        all_files = os.listdir(res_dir)
        sum_files = []
        for file in all_files:
            if "summary" in file:
                if re.search(self.file_ext, file):
                    sum_files.append(f"{res_dir}/{file}")

        files_dict = {}

        for engine in self.engines:
            eng_files = []
            for file in sum_files:
                if engine in file:
                    if re.search(self.file_ext, file):
                        eng_files.append(file)
            files_dict[engine] = eng_files

        return files_dict

    def _set_network(self):
        """set the network based on the network file or based on all the found files."""

        # get perturbations and ligands for the network
        # if network file, get from that
        # if not, use all results files for this

        if not self.output_folder:
            if not self._net_file:
                logging.error(
                    "As there is no provided results directory or network file, please set perturbations and ligands manually."
                )
                return
            else:
                file_names = None

        else:
            # get results files from dict into a list, flatten the list
            results_lists = list(self._results_files.values()) + list(
                self._results_repeat_files.values()
            )
            file_names = [res_file for sublist in results_lists for res_file in sublist]

        # for all
        values = get_info_network(
            results_files=file_names,
            net_file=self._net_file,
            extra_options={"engines": self.engines},
        )

        self.perturbations = values[0]
        self.ligands = values[1]

        # for individual engines
        for engine in self.engines:
            values = get_info_network(
                results_files=file_names,
                net_file=self._net_file,
                extra_options={"engines": engine},
            )

            # check if there are any values for this engine, if not assume it is the entire network from before
            # as dont want this to be none in later functions
            if not values[0]:
                self._perturbations_dict[engine] = self.perturbations
                self._ligands_dict[engine] = self.ligands
            else:
                self._perturbations_dict[engine] = values[0]
                self._ligands_dict[engine] = values[1]

    def _set_dictionary_outputs(self):
        # set all the dicts for analysis
        # per engine dicts (also used for other results)
        self.calc_pert_dict = {}  # diff from the results repeat files, average
        self.calc_bound_dict = {}
        self.calc_free_dict = {}
        self.cinnabar_calc_val_dict = {}  # from the cinnabar network analysis
        self.cinnabar_exper_val_dict = (
            {}
        )  # normalised from the cinnabar network analysis
        self.cinnabar_calc_pert_dict = {}  # from cinnabar network edges
        self.cinnabar_exper_pert_dict = {}  # from cinnabar network edges

        # solo dicts for exper
        self.exper_val_dict = None  # yml converted into experimental values, actual, for ligands in object
        self.normalised_exper_val_dict = (
            None  # yml converted into experimental values, then normalised
        )
        self.exper_pert_dict = None  # yml converted into experimental values, actual, for perturbations in object

        # for convergence
        self.spert_results_dict = {}
        self.spert_bound_dict = {}
        self.spert_free_dict = {}
        self.epert_results_dict = {}
        self.epert_bound_dict = {}
        self.epert_free_dict = {}

        # storing the nx digraphs, per engine
        self._cinnabar_networks = {}
        # overall graph
        self.ligands_folder = None
        # cycles
        self.cycle_dict = {}

        # for other results
        self.other_results_names = []

        # for checking against free energy workflows
        self._fwf_experimental_DDGs = None
        self._fwf_computed_relative_DDGs = {}
        self._fwf_path = None

        # for plotting
        self._plotting_object = None
        self._histogram_object = None
        # for stats
        self._stats_object = None

    def get_experimental(
        exp_file: str = None, pert_val: str = None, normalised: bool = False
    ) -> dict:
        """get experimental values converted. Can be used to add the experimental values to the object.

        Args:
            exp_file (str, optional): path to experimental yml or csv file. Defaults to None.
            pert_val (str, optional): Whether to return dict for the pert or the vals. Defaults to None.
            normalised (bool, optional): Whether the vals should be normalised. Defaults to False.

        Returns:
            dict: experimental dictionary. Only if pert_val is provided, else returns nothing.
        """

        if pert_val:
            pert_val = validate.pert_val(pert_val)
        normalised = validate.boolean(normalised)

        if not exp_file:
            if not self.exp_file:
                logging.error(
                    "need an experimental file to proceed with most of the calculations. please set using self.get_experimental(file)"
                )
            else:
                exp_file = self.exp_file
        else:
            self.exp_file = validate.file_path(exp_file)

        self._get_experimental()

        if not pert_val:
            return

        elif pert_val == "pert":
            return self.exper_pert_dict

        elif pert_val == "val":
            if normalised:
                return self.normalised_exper_val_dict

            else:
                return self.exper_val_dict

    def _get_experimental(self):
        """get the experimental value dictionaries from a given yml file."""

        exp_file = self.exp_file

        if exp_file.split(".")[-1] == "yml":
            exper_val_dict = convert.yml_into_exper_dict(
                exp_file, temperature=self.temperature
            )  # this output is in kcal/mol
        elif exp_file.split(".")[-1] == "csv":
            exper_val_dict = convert.csv_into_exper_dict(
                exp_file, temperature=self.temperature
            )  # this output is in kcal/mol
        else:
            logging.error(
                "file type for experimental must be yml or csv. No experimental results will be analysed"
            )
            return

        # experimental value dict
        self.exper_val_dict = make_dict.exper_from_ligands(exper_val_dict, self.ligands)

        # normalise the experimental values
        self.normalised_exper_val_dict = make_dict.exper_from_ligands(
            exper_val_dict, self.ligands, normalise=True
        )

        self.exper_pert_dict = make_dict.exper_from_perturbations(
            self.exper_val_dict, self.perturbations
        )

    def _validate_in_names_list(self, name: str, make_list: bool = False) -> list:
        """validate if the name is in the names list

        Args:
            name (str): the name to validate
            make_list (bool, optional): whether to make into a list, default is False.

        Raises:
            ValueError: if not in names list

        Returns:
            str: the validated name or names list.
        """

        names = validate.is_list(name, make_list=True)

        for name in names:
            name = validate.string(name)
            if name not in (self.engines + self.other_results_names):
                raise ValueError(
                    f"{name} must be in {self.engines + self.other_results_names}"
                )

        if not make_list:
            names = names[0]

        return names

    def remove_perturbations(self, perts: list, name: Optional[str] = None):
        """remove perturbations from the network used.

        Args:
            perts (list): list of perturbations to remove.
        """

        perts = validate.is_list(perts, make_list=True)

        if not name:
            for pert in perts:
                if pert in self.perturbations:
                    self.perturbations.remove(pert)
        else:
            name = self._validate_in_names_list(name)
            for pert in perts:
                if pert in self._perturbations_dict[name]:
                    self._perturbations_dict[name].remove(pert)

        if self._is_computed_dicts:
            # recompute dicts
            self._compute_dicts()
        # remove plotting object as needs to be reintialised with new perturbations
        self._plotting_object = None
        self._histogram_object = None
        self._stats_object = None

    def remove_ligands(self, ligs: list, name: Optional[str] = None):
        """remove ligand and assosciated perturbations from the network used.

        Args:
            ligs (list): list of ligands to remove.
        """

        ligs = validate.is_list(ligs, make_list=True)

        if not name:
            for lig in ligs:
                if lig in self.ligands:
                    self.ligands.remove(lig)

                for pert in self.perturbations:
                    if lig in pert:
                        self.perturbations.remove(pert)
        else:
            name = self._validate_in_names_list(name)
            for lig in ligs:
                if lig in self._ligands_dict[name]:
                    self._ligands_dict[name].remove(lig)

                for pert in self._perturbations_dict[name]:
                    if lig in pert:
                        self._perturbations_dict[name].remove(pert)

        if self._is_computed_dicts:
            self._compute_dicts()
        # remove plotting object as needs to be reintialised with new perturbations
        self._plotting_object = None
        self._histogram_object = None
        self._stats_object = None

    def change_name(self, old_name: str, new_name: str):
        """change the name of the data. Can be used for engine or other result names. Will update the self.dicts with this new name.

        Args:
            old_name (str): old name to replace
            new_name (str): new name to replace it with
        """

        dict_list = [
            self._perturbations_dict,
            self._ligands_dict,
            self.calc_pert_dict,
            self.calc_bound_dict,
            self.calc_free_dict,
            self.cinnabar_calc_val_dict,
            self.cinnabar_exper_val_dict,
            self.cinnabar_calc_pert_dict,
            self.cinnabar_exper_pert_dict,
            self._cinnabar_networks,
            self.spert_results_dict,
            self.spert_bound_dict,
            self.spert_free_dict,
            self.epert_results_dict,
            self.epert_bound_dict,
            self.epert_free_dict,
            self._fwf_computed_relative_DDGs,
            self._results_value_files,
        ]

        for adict in dict_list:
            try:
                adict[new_name] = adict.pop(old_name)
            except:
                logging.error(
                    f"could not rename one of the dicts, as it does not have this key as one of its keys."
                )

        if old_name in self.other_results_names:
            self.other_results_names.remove(old_name)
            self.other_results_names.append(new_name)
        elif old_name in self.engines:
            self.engines.remove(old_name)
            self.other_results_names.append(new_name)

        # remove plotting object as needs to be reintialised with new name
        self._plotting_object = None
        self._histogram_object = None
        self._stats_object = None

    def compute_results(self):
        """compute the dictionaries for analysis and those passed to the plotting object.

        Args:
            cycle_closure (bool, optional): whether to compute the cycle closures. Defaults to True.
            statistics (bool, optional): whether to calculate the statistics. Defaults to True.
        """

        # get all the dictionaries needed for plotting
        self._compute_dicts()

        # initialise plotting and stats objects
        self._initialise_plotting_object()

        self._is_computed_dicts = True

    def compute_cycle_closures(self, engines: list = None):
        """Compute cycle closures average and error.

        Args:
            engines (list): MD engines to compute cycle closures for.
        """
        if engines:
            engines = self._validate_in_names_list(engines)
        else:
            engines = self.engines

        for engine in engines:
            self._compute_cycle_closures(engine)
            print(
                f"cycle closure average is {self.cycle_dict[engine][1]} +/- {self.cycle_dict[engine][2]} kcal/mol"
            )

    def _compute_cycle_closures(self, engine: str):
        """compute the cycle closures and their stats for each engine for the network.

        Returns:
            dict: self.cycle_dict (engine: cycles_dict, cycle_vals, np.mean(cycle_vals), np.std(cycle_vals) )
        """

        engine = self._validate_in_names_list(engine)

        graph = network_graph(
            self._ligands_dict[engine],
            self._perturbations_dict[engine],
            self.calc_pert_dict[engine],
        )

        cycles = graph.cycle_closure_dict()
        avg_cycle_closures = graph.average_cycle_closures()

        self.cycle_dict.update(
            {engine: (cycles, avg_cycle_closures[0], avg_cycle_closures[1])}
        )  # the cycles dict : cycles_dict(cycle:(val,err)), mean, deviation

    def _compute_dicts(self):
        """calculate the perturbation dicts from the previously passed repeat files."""

        # compute the experimental for perturbations
        self._get_experimental()  # get experimental val dict and normalised dict

        # for self plotting of per pert
        for engine in (
            self.engines + self.other_results_names
        ):  # other results will only be added after already computed once in function below
            if not self._results_files[engine]:
                files = self._results_repeat_files[engine]
            else:
                files = self._results_files[engine]

            calc_diff_dict = make_dict.comp_results(
                files, self._perturbations_dict[engine], engine, method=self.method
            )  # older method
            self.calc_pert_dict.update({engine: calc_diff_dict})

            # calc the free and bound leg dicts for the engine
            try:
                calc_bound_dict = make_dict.comp_results(
                    self._results_bound_repeat_files[engine],
                    self._perturbations_dict[engine],
                    engine,
                    method=self.method,
                )  # older method
                self.calc_bound_dict.update({engine: calc_bound_dict})
                calc_free_dict = make_dict.comp_results(
                    self._results_free_repeat_files[engine],
                    self._perturbations_dict[engine],
                    engine,
                    method=self.method,
                )  # older method
                self.calc_free_dict.update({engine: calc_free_dict})
            except Exception as e:
                logging.error(e)
                logging.error("Could not calculate dicts for bound/free legs.")

            self._compute_cinnabar_dict(files, engine, method=self.method)

    def _compute_cinnabar_dict(
        self, files: list, engine: str, method: Optional[str] = None
    ):
        """compute cinnabar and get the dictionaries from it."""

        perts, ligs = get_info_network_from_dict(self.calc_pert_dict[engine])

        # get the files into cinnabar format for analysis
        cinnabar_file_name = (
            f"{self.results_folder}/cinnabar_{engine}_{self.file_ext}_{self.net_ext}"
        )

        convert.cinnabar_file(
            files,
            self.exper_val_dict,
            cinnabar_file_name,
            perturbations=perts,
            method=method,
        )

        try:
            # compute the per ligand for the network
            network = wrangle.FEMap(f"{cinnabar_file_name}.csv")
            self._cinnabar_networks.update({engine: network})

            # from cinnabar graph
            self.cinnabar_calc_pert_dict.update(
                {engine: make_dict.from_cinnabar_network_edges(network, "calc", perts)}
            )
            self.cinnabar_exper_pert_dict.update(
                {engine: make_dict.from_cinnabar_network_edges(network, "exp", perts)}
            )

            # for self plotting of per ligand
            self.cinnabar_calc_val_dict.update(
                {engine: make_dict.from_cinnabar_network_node(network, "calc")}
            )
            self.cinnabar_exper_val_dict.update(
                {
                    engine: make_dict.from_cinnabar_network_node(
                        network, "exp", normalise=True
                    )
                }
            )

            write_vals_file(
                self.cinnabar_exper_val_dict[engine],
                f"{self.results_folder}/lig_values_{engine}_{self.file_ext}_{self.net_ext}",
                engine,
                self.file_ext,
                self.method,
            )
            self._results_value_files[engine] = [
                f"{self.results_folder}/lig_values_{engine}_{self.file_ext}_{self.net_ext}.csv"
            ]

        except Exception as e:
            logging.error(e)
            logging.error(f"could not create cinnabar network for {engine}")
            self._cinnabar_networks.update({engine: None})
            self.cinnabar_calc_pert_dict.update({engine: None})
            self.cinnabar_exper_pert_dict.update({engine: None})
            self.cinnabar_calc_val_dict.update({engine: None})
            self.cinnabar_exper_val_dict.update({engine: None})

    def compute_other_results(
        self,
        file_names: Optional[list] = None,
        name: Optional[str] = None,
        method: Optional[str] = None,
        bound_files: Optional[list] = None,
        free_files: Optional[list] = None,
    ):
        """compute other results in a similar manner to the engine results.

        Args:
            file_names (list, optional): list of other results. Defaults to None.
            name (str, optional): name of these other results (for files and graphs and identification). Defaults to None.
            method (str, optional): method in the input files to include only. Defaults to None.
        """

        file_names = validate.is_list(file_names, make_list=True)
        for file in file_names:
            validate.file_path(file)

        # add identifier for the other results
        name = validate.string(name)
        if name in self.other_results_names:
            logging.error(
                f"{name} is already in the other results. please use a different name!"
            )
            return
        else:
            self.other_results_names.append(name)

        # add files to file list
        self._results_repeat_files[name] = file_names

        new_file_path = f"{file_names[0].replace(file_names[0].split('/')[-1], '')[:-1]}/{name}_results_file"

        # for self plotting of per pert
        calc_diff_dict = make_dict.comp_results(
            file_names,
            perturbations=None,
            engine=None,
            name=name,
            method=method,
            output_file=new_file_path,
        )
        # set info to dicts etc
        self.calc_pert_dict.update({name: calc_diff_dict})
        self._results_files[name] = [f"{new_file_path}.csv"]
        perts, ligs = get_info_network_from_dict(calc_diff_dict)
        self._perturbations_dict[name] = perts
        self._ligands_dict[name] = ligs

        if bound_files and free_files:
            bound_files = validate.is_list(bound_files, make_list=True)
            free_files = validate.is_list(free_files, make_list=True)
            for file in bound_files + free_files:
                validate.file_path(file)
            calc_bound_dict = make_dict.comp_results(
                bound_files, perts, engine=None, name=name, method=method
            )
            self.calc_bound_dict.update({name: calc_bound_dict})
            self._results_bound_repeat_files.update({name: bound_files})
            calc_free_dict = make_dict.comp_results(
                free_files, perts, engine=None, name=name, method=method
            )
            self.calc_free_dict.update({name: calc_free_dict})
            self._results_free_repeat_files.update({name: free_files})

        else:
            self.calc_free_dict.update({name: None})
            self.calc_bound_dict.update({name: None})

        self._compute_cinnabar_dict(files=f"{new_file_path}.csv", engine=name)

        # initialise plotting and stats objects again so theyre added
        self._initialise_plotting_object(check=False)
        self._initialise_stats_object(check=False)
        
    def compute_convergence(self, compute_missing: bool = False):
        """compute the convergence of the results.

        Args:
            compute_missing (bool, optional): Compute the missing convergence results. This can take awhile! Defaults to False.
        """

        compute_missing = validate.boolean(compute_missing)

        for engine in self.engines:
            self.spert_results_dict[engine] = {}
            self.spert_bound_dict[engine] = {}
            self.spert_free_dict[engine] = {}
            self.epert_results_dict[engine] = {}
            self.epert_bound_dict[engine] = {}
            self.epert_free_dict[engine] = {}

            for pert in self.perturbations:
                # find correct path, use extracted if it exists
                if self.method:
                    name = f"_{self.method}"
                else:
                    name = ""
                path_to_dir = f"{self.output_folder}/{engine}/{pert}{name}/pickle"
                try:
                    validate.folder_path(path_to_dir)
                except:
                    logging.error(
                        f"cannot find pickle directory for {pert} in {engine}, does '{path_to_dir}' exist?"
                    )
                    path_to_dir = None

                try:
                    pickle_ext = analyse.pickle_ext(
                        self.analysis_options.dictionary(), pert, engine
                    ).split("truncate")[0]

                    with open(
                        f"{path_to_dir}/spert_results_dict_{pickle_ext}.pickle", "rb"
                    ) as file:
                        sresults_dict = pickle.load(file)
                    with open(
                        f"{path_to_dir}/epert_results_dict_{pickle_ext}.pickle", "rb"
                    ) as file:
                        eresults_dict = pickle.load(file)
                    with open(
                        f"{path_to_dir}/spert_bound_dict_{pickle_ext}.pickle", "rb"
                    ) as file:
                        sbound_dict = pickle.load(file)
                    with open(
                        f"{path_to_dir}/epert_bound_dict_{pickle_ext}.pickle", "rb"
                    ) as file:
                        ebound_dict = pickle.load(file)
                    with open(
                        f"{path_to_dir}/spert_free_dict_{pickle_ext}.pickle", "rb"
                    ) as file:
                        sfree_dict = pickle.load(file)
                    with open(
                        f"{path_to_dir}/epert_free_dict_{pickle_ext}.pickle", "rb"
                    ) as file:
                        efree_dict = pickle.load(file)

                    self.spert_results_dict[engine][pert] = sresults_dict
                    self.spert_bound_dict[engine][pert] = sbound_dict
                    self.spert_free_dict[engine][pert] = sfree_dict
                    self.epert_results_dict[engine][pert] = eresults_dict
                    self.epert_bound_dict[engine][pert] = ebound_dict
                    self.epert_free_dict[engine][pert] = efree_dict

                    pickle_loaded = True

                except:
                    logging.error(
                        f"could not load pickles for {pert} in {engine}. Was it analysed for convergence?"
                    )

                    pickle_loaded = False

                if compute_missing and not pickle_loaded:
                    path_to_dir = f"{self.output_folder}/{engine}/{pert}{name}"
                    try:
                        validate.folder_path(path_to_dir)
                    except:
                        path_to_dir = None
                        logging.error(
                            f"{engine} {pert}{name} does not exist in the searched output locations."
                        )
                        continue

                    analysed_pert = analyse(
                        path_to_dir, pert=pert, analysis_prot=self.analysis_options
                    )
                    analysed_pert._save_pickle = True

                    logging.info(
                        f"calculating convergence for {engine} {pert}{name} ..."
                    )
                    analysed_pert.calculate_convergence()

                    pickle_ext = analyse.pickle_ext(
                        self.analysis_options.dictionary(), pert, engine
                    ).split("truncate")[0]

                    with open(
                        f"{path_to_dir}/pickle/spert_results_dict_{pickle_ext}.pickle",
                        "rb",
                    ) as file:
                        sresults_dict = pickle.load(file)
                    with open(
                        f"{path_to_dir}/pickle/epert_results_dict_{pickle_ext}.pickle",
                        "rb",
                    ) as file:
                        eresults_dict = pickle.load(file)
                    with open(
                        f"{path_to_dir}/pickle/spert_bound_dict_{pickle_ext}.pickle",
                        "rb",
                    ) as file:
                        sbound_dict = pickle.load(file)
                    with open(
                        f"{path_to_dir}/pickle/epert_bound_dict_{pickle_ext}.pickle",
                        "rb",
                    ) as file:
                        ebound_dict = pickle.load(file)
                    with open(
                        f"{path_to_dir}/pickle/spert_free_dict_{pickle_ext}.pickle",
                        "rb",
                    ) as file:
                        sfree_dict = pickle.load(file)
                    with open(
                        f"{path_to_dir}/pickle/epert_free_dict_{pickle_ext}.pickle",
                        "rb",
                    ) as file:
                        efree_dict = pickle.load(file)

                    self.spert_results_dict[engine][pert] = sresults_dict
                    self.spert_bound_dict[engine][pert] = sbound_dict
                    self.spert_free_dict[engine][pert] = sfree_dict
                    self.epert_results_dict[engine][pert] = eresults_dict
                    self.epert_bound_dict[engine][pert] = ebound_dict
                    self.epert_free_dict[engine][pert] = efree_dict

    def successful_perturbations(
        self, engine: str, perts: Optional[list] = None
    ) -> tuple:
        """calculate how many successful runs

        Args:
            engine (str): the engine to calc for
            perts (list, optional): The list of perts to use. Defaults to None.

        Returns:
            tuple: (val, percen, perturbations)
        """

        res_dict = self.calc_pert_dict[engine]
        engine = validate.engine(engine)
        if perts:
            perts = validate.is_list(perts)
        else:
            perts = res_dict.keys()

        perturbations = []
        if self._is_computed_dicts:
            val = 0
            for key in res_dict.keys():
                if key in perts:
                    if not isnan(res_dict[key][0]):
                        val += 1
                        perturbations.append(key)

            percen = (val / len(perts)) * 100

            logging.info(
                f"{val} out of {len(perts)} have results, which is {percen} %."
            )
            return (val, percen, perturbations)

        else:
            logging.error("please compute results from results files first.")
            return (None, None, None)

    def failed_perturbations(self, engine: str) -> list:
        """get the failed perturbations for an engine.

        Args:
            engine (str): the MD engine.

        Returns:
            list: failed perturbations.
        """

        engine = validate.engine(engine)

        val, percen, perturbations = self.successful_perturbations(engine)

        failed_perts = []

        for pert in self.perturbations:
            if pert not in perturbations:
                failed_perts.append(pert)

        return failed_perts

    def draw_failed_perturbations(self, engine: str):
        """Draw the failed perturbations.

        Args:
            engine (str): the MD engine.
        """

        engine = validate.engine(engine)

        perturbations = self.failed_perturbations(engine)

        if perturbations:
            graph = network_graph(
                self._ligands_dict[engine],
                self._perturbations_dict[engine],
                self.calc_pert_dict[engine],
                ligands_folder=self.ligands_folder,
            )
            for pert in perturbations:
                graph.draw_perturbation(pert)

    def draw_perturbations(self, pert_list: list):
        """Draw certain perturbations.

        Args:
            engine (str): the MD engine.
        """
        graph = network_graph(
            self.ligands,
            self.perturbations,
            ligands_folder=self.ligands_folder,
        )
        for pert in pert_list:
            graph.draw_perturbation(pert)

    def disconnected_ligands(self, engine: str) -> list:
        """Get the disconnected ligands.

        Args:
            engine (str): The MD engine.

        Returns:
            list: List of disconnected ligands.
        """

        engine = validate.engine(engine)
        val, percen, perturbations = self.successful_perturbations(engine)
        graph = network_graph(
            self._ligands_dict[engine],
            perturbations,
            ligands_folder=self.ligands_folder,
        )
        ligs = graph.disconnected_ligands()

        return ligs

    def draw_disconnected_ligands(self, engine: str):
        """Draw the disconnected ligands.

        Args:
            engine (str): The MD engine.
        """

        ligands = self.disconnected_ligands(engine)
        graph = network_graph(
            self._ligands_dict[engine],
            self._perturbations_dict[engine],
            ligands_folder=self.ligands_folder,
        )
        for lig in ligands:
            graph.draw_ligand(lig)

    def get_outliers(
        self, threshold: Union[int, float] = 10, name: Optional[str] = None
    ) -> list:
        """get outliers above a certain difference to the experimental.

        Args:
            threshold (float, optional): difference threshold above which to remove. Defaults to 10.
            name (str, optional): name of the data (engine or other results). Defaults to None.
        """

        # can get from dict or dataframe
        # probably best from plotting object

        plot_obj = self._initialise_plotting_object(check=True)
        threshold = validate.is_float(threshold)

        perts = []

        if name:
            names = plot_obj._validate_in_names_list(name, make_list=True)
        else:
            names = self.other_results_names + self.engines

        for name in names:
            freenrg_df_plotting = plot_obj.freenrg_df_dict["experimental"][name][
                "pert"
            ].dropna()
            x = freenrg_df_plotting[f"freenrg_experimental"]
            y = freenrg_df_plotting["freenrg_calc"]
            # get an array of the MUE values comparing experimental and FEP values. Take the absolute values.
            mue_values = abs(x - y)

            # find the n ligand names that are outliers.
            perts = [
                key
                for pert, key in zip(
                    mue_values.gt(threshold), mue_values.gt(threshold).keys()
                )
                if pert
            ]

        return perts

    def remove_outliers(
        self, threshold: Union[int, float] = 10, name: Optional[str] = None
    ):
        """remove outliers above a certain difference to the experimental.

        Args:
            threshold (float, optional): difference threshold above which to remove. Defaults to 10.
            name (str, optional): name of the data (engine or other results). Defaults to None.
        """

        perts = self.get_outliers(threshold, name)

        for pert in perts:
            self.remove_perturbations(pert, name=name)

            logging.info(f"removed {pert} from perturbations as outlier for {name}.")

        self._compute_dicts()

        # remove plotting object as needs to be reintialised with new perturbations
        self._plotting_object = None
        self._histogram_object = None
        self._stats_object = None

    def sort_ligands_by_binding_affinity(self, engine: Optional[str] = None) -> list:
        """Ligands sorted by binding affinity.

        Args:
            engine (Optional[str], optional): _description_. Defaults to None.

        Returns:
            list: _description_
        """
        if not self._is_computed_dicts:
            self._compute_dicts()

        if engine:
            engine = validate.engine(engine)
        else:
            self.compute_consensus()
            logging.info(
                "sorting ligands based on consensus scoring, as no engine was provided."
            )
            engine = f"consensus_{'_'.join(str(engine) for engine in self.engines)}"

        df = pd.DataFrame(
            self.cinnabar_calc_val_dict[engine], index=["value", "error"]
        ).transpose()

        return df.sort_values(by="value", ascending=True)

    def sort_ligands_by_experimental_binding_affinity(self):
        df = pd.DataFrame(self.exper_val_dict, index=["value", "error"]).transpose()

        return df.sort_values(by="value", ascending=True)

    def compute_consensus(self, names: Optional[list] = None):
        if names:
            self._validate_in_names_list(names, make_list=True)
        else:
            names = self.engines

        consensus_pert_dict = {}

        for pert in self.perturbations:
            consensus_pert_dict[pert] = []

            for engine in names:
                consensus_pert_dict[pert].append(self.calc_pert_dict[engine][pert][0])

            consensus_pert_dict[pert] = (
                np.mean(consensus_pert_dict[pert]),
                sem(consensus_pert_dict[pert]),
            )

            df = pd.DataFrame.from_dict(
                consensus_pert_dict, orient="index", columns=["freenrg", "error"]
            )
            df.index.name = "perturbations"
            df = df.reset_index()
            df[["lig_0", "lig_1"]] = df["perturbations"].str.split("~", expand=True)
            df["engine"] = f"consensus_{'_'.join(str(engine) for engine in names)}"
            df["analysis"] = self.file_ext
            df["method"] = "None"
            df = df.drop(labels="perturbations", axis=1)
            df = df[
                ["lig_0", "lig_1", "freenrg", "error", "engine", "analysis", "method"]
            ]

            df.to_csv(
                f"{self.results_folder}/consensus_score_{self.engines}_{self.file_ext}.csv",
                sep=",",
                index=False,
            )

        self.compute_other_results(
            f"{self.results_folder}/consensus_score_{self.engines}_{self.file_ext}.csv",
            name="consensus",
        )

    def add_ligands_folder(self, folder: str):
        """add a ligands folder so the ligands can be visualised

        Args:
            folder (str): ligand file location folder
        """

        self.ligands_folder = validate.folder_path(folder)

    def draw_graph(
        self,
        use_cinnabar: bool = False,
        engines: Optional[list] = None,
        successful_perturbations: bool = True,
    ):
        """draw the network graph.

        Args:
            use_cinnabar (bool): whether to use the cinnabar data or the self computed data. Defaults to False.
            engines (str/list, optional): engine to draw the network for. Defaults to None, draws for each engine.
            successful_perturbations (bool): whether to only draw the successful runs. Only useable if cinnabar is set to False. Defaults to True.
        """

        if engines:
            engines = self._validate_in_names_list(engines, make_list=True)
        else:
            engines = self.engines

        if use_cinnabar:
            for engine in engines:
                file_name = f"{self.graph_dir}/cinnabar_network_{engine}_{self.file_ext}_{self.net_ext}.png"
                self._cinnabar_networks[engine].draw_graph(file_name=file_name)

        else:
            successful_perturbations = validate.boolean(successful_perturbations)

            file_dir = validate.folder_path(
                f"{self.results_folder}/network", create=True
            )

            if successful_perturbations:
                for engine in engines:
                    val, percen, perturbations = self.successful_perturbations(engine)
                    graph = network_graph(
                        self.ligands,
                        perturbations,
                        self.calc_pert_dict[engine],
                        file_dir=file_dir,
                        ligands_folder=self.ligands_folder,
                    )
                    graph.draw_graph(title=engine)
            else:
                for engine in engines:
                    graph = network_graph(
                        self._ligands_dict[engine],
                        self._perturbations_dict[engine],
                        self.calc_pert_dict[engine],
                        file_dir=file_dir,
                        ligands_folder=self.ligands_folder,
                    )
                    graph.draw_graph(title=engine)

    def _initialise_plotting_object(self, check: bool = False):
        """intialise the plotting object

        Args:
            check (bool, optional): whether to check the plotting object. Defaults to False.

        Returns:
            pipeline.analysis.plotting_engines: the plotting object.
        """

        # if not checking, always make
        if not check:
            self._plotting_object = plotting_engines(analysis_object=self)

        # if checking, first see if it exists and if not make
        elif check:
            if not self._plotting_object:
                self._plotting_object = plotting_engines(analysis_object=self)

        return self._plotting_object

    def _initialise_histogram_object(self, check: bool = False):
        """intialise the histogram plotting object

        Args:
            check (bool, optional): whether to check the plotting histogram object. Defaults to False.

        Returns:
            pipeline.analysis.plotting_engines: the plotting histogram object.
        """

        # if not checking, always make
        if not check:
            self._histogram_object = plotting_histogram(analysis_object=self)

        # if checking, first see if it exists and if not make
        elif check:
            if not self._histogram_object:
                self._histogram_object = plotting_histogram(analysis_object=self)

        return self._histogram_object

    def plot_bar_ddG(self, engines: Optional[list] = None, **kwargs):
        """plot the bar plot of the perturbations.

        Args:
            engines (str, optional): engine to plot for. Defaults to None, will use all.
        """
        if engines:
            engines = self._validate_in_names_list(engines, make_list=True)
        else:
            engines = self.engines + self.other_results_names
        engines.append("experimental")

        plot_obj = self._initialise_plotting_object(check=True)
        plot_obj.bar(pert_val="pert", names=engines, **kwargs)

    def plot_bar_dG(self, engines: Optional[list] = None, **kwargs):
        """plot the bar plot of the values per ligand.

        Args:
            engines (str, optional): engine to plot for. Defaults to None, will use all.
        """
        if engines:
            engines = self._validate_in_names_list(engines, make_list=True)
        else:
            engines = self.engines + self.other_results_names
        engines.append("experimental")

        plot_obj = self._initialise_plotting_object(check=True)
        plot_obj.bar(pert_val="val", names=engines, **kwargs)

    def plot_bar_leg(self, engine, leg="bound", **kwargs):
        engine = validate.is_list(engine, make_list=True)

        plot_obj = self._initialise_plotting_object(check=True)

        plotting_dict = {
            "title": f"{leg} for {self.file_ext.replace('_',',')}, {self.net_ext.replace('_',',')}"
        }
        for key, value in kwargs.items():
            plotting_dict[key] = value

        plot_obj.bar(pert_val=leg, names=engine, **plotting_dict)

    def plot_scatter_ddG(
        self, engines: Optional[list] = None, use_cinnabar: bool = False, **kwargs
    ):
        """plot the scatter plot of the perturbations.

        Args:
            engines (str, optional): engine to plot for. Defaults to None, will use all.
            use_cinnabar (bool, optional): whether to plot via cinnabar. Defaults to False.
        """

        if engines:
            engines = self._validate_in_names_list(engines, make_list=True)
        else:
            engines = self.engines + self.other_results_names

        if use_cinnabar:
            for engine in engines:
                plotting.plot_DDGs(
                    self._cinnabar_networks[engine].graph,
                    filename=f"{self.graph_dir}/DDGs_{engine}_{self.file_ext}_{self.net_ext}.png",
                    title=f"DDGs for {engine}, {self.net_ext}",
                    **{"figsize": 5},
                )  # with {self.file_ext}

        else:
            plot_obj = self._initialise_plotting_object(check=True)
            plot_obj.scatter(pert_val="pert", y_names=engines, **kwargs)

    def plot_scatter_dG(
        self, engines: Optional[list] = None, use_cinnabar: bool = False, **kwargs
    ):
        """plot the scatter plot of the values per ligand.

        Args:
            engines (str, optional): engine to plot for. Defaults to None, will use all.
            use_cinnabar (bool, optional): whether to plot via cinnabar. Defaults to False.
        """

        if engines:
            engines = self._validate_in_names_list(engines, make_list=True)
        else:
            engines = self.engines + self.other_results_names

        if use_cinnabar:
            for engine in engines:
                plotting.plot_DGs(
                    self._cinnabar_networks[engine].graph,
                    filename=f"{self.graph_dir}/DGs_{engine}_{self.file_ext}_{self.net_ext}.png",
                    title=f"DGs for {engine}, {self.net_ext}",
                    **{"figsize": 5},
                )

        else:
            plot_obj = self._initialise_plotting_object(check=True)
            plot_obj.scatter(pert_val="val", y_names=engines, **kwargs)

    def plot_eng_vs_eng(
        self,
        engine_a: str = None,
        engine_b: str = None,
        pert_val: str = "pert",
        **kwargs,
    ):
        """plot scatter plot of engine_a vs engine_b

        Args:
            engine_a (str): engine_a. Defaults to None.
            engine_b (str): engine_b. Defaults to None.
            pert_val (str): whether perturbations 'pert' or values per ligand 'val'. Defaults to "pert".
        """

        plot_obj = self._initialise_plotting_object(check=True)

        if pert_val == "pert":
            binding = "$\Delta\Delta$G$_{bind}$ (kcal/mol)"
        elif pert_val == "val":
            binding = "$\Delta$G$_{bind}$ (kcal/mol)"
        plotting_dict = {
            "title": f"{engine_a} vs {engine_b}\n for {self.file_ext}, {self.net_ext}",
            "y label": f"{engine_a} " + binding,
            "x label": f"{engine_b} " + binding,
            "key": False,
        }
        for key, value in kwargs.items():
            plotting_dict[key] = value

        plot_obj.scatter(
            pert_val=pert_val, y_names=engine_a, x_name=engine_b, **plotting_dict
        )

    def plot_outliers(
        self,
        engines: Optional[list] = None,
        no_outliers: int = 5,
        pert_val: str = "pert",
        **kwargs,
    ):
        """plot scatter plot with annotated outliers.

        Args:
            engine (list, optional): engine to plot for. Defaults to None.
            no_outliers (int, optional): number of outliers to annotate. Defaults to 5.
            pert_val (str, optional): whether plotting 'pert' ie perturbations or 'val' ie values (per ligand result). Defaults to None.

        """

        plot_obj = self._initialise_plotting_object(check=True)
        plot_obj.scatter(
            pert_val=pert_val,
            y_names=engines,
            outliers=True,
            no_outliers=no_outliers,
            **kwargs,
        )

    def plot_histogram_sem(
        self, engines: Optional[list] = None, pert_val: str = "pert"
    ):
        """plot histograms for the sem of the result (either pert or val).

        Args:
            engines (list, optional): engines to plot for. Defaults to None.
            pert_val (str, optional): whether perturbations 'pert' or values per ligand 'val'. Defaults to "pert".
        """

        if pert_val == "pert":
            type_error = ["SEM_pert"]
        elif pert_val == "val":
            type_error = ["per_lig"]
        else:
            raise ValueError("pert_val must be 'pert' or 'val'")

        self._plot_histogram(engines, type_error)

    def plot_histogram_legs(self, engines: Optional[list] = None):
        """plot histograms for the errors per leg.

        Args:
            engines (list, optional): engines to plot for. Defaults to None.
        """

        self._plot_histogram(engines, ["bound", "free"])

    def plot_histogram_repeats(self, engines: Optional[list] = None):
        """plot histograms for the errors per repeat.

        Args:
            engines (list, optional): engines to plot for. Defaults to None.
        """
        self._plot_histogram(engines, ["repeat"])

    def _plot_histogram(self, engines: Optional[list] = None, type_errors: str = None):
        """internal function for plotting histograms"""

        hist_obj = self._initialise_histogram_object(check=True)

        if not engines:
            engines = self.engines
        else:
            engines = validate.is_list(engines, make_list=True)

        for type_error in type_errors:
            for engine in engines:
                hist_obj.histogram(name=engine, type_error=type_error)
            hist_obj.histogram_distribution(names=engines, type_error=type_error)

    def plot_convergence(self, engines: Optional[list] = None):
        if not self.spert_results_dict:
            raise EnvironmentError(
                f"please run 'calculate_convergence' first with the main_dir set."
            )

        else:
            if not engines:
                engines = self.engines
            else:
                engines = validate.engines(engines)

            plot_obj = self._initialise_plotting_object(check=True)
            plot_obj.plot_convergence(engines=engines)

    def _initialise_stats_object(self, check: bool = False):
        """intialise the object for statistical analysis.

        Args:
            check (bool, optional): whether to check. Defaults to False.

        Returns:
            pipeline.analysis.stats_engines: statistics object
        """

        # if not checking, always make
        if not check:
            self._stats_object = stats_engines(analysis_object=self)

        # if checking, first see if it exists and if not make
        elif check:
            if not self._stats_object:
                self._stats_object = stats_engines(analysis_object=self)

        return self._stats_object

    def calc_mae_engines(self, pert_val: str = None, engines: Optional[list] = None):
        """calculate the Mean Absolute Error (MAE) vs experimental results

        Args:
            pert_val (str, optional): whether plotting 'pert' ie perturbations or 'val' ie values (per ligand result). Defaults to None.
            engines (list, optional): names of engines / other results names to calculate for.

        Returns:
            tuple: of dataframe of value and error (df, df_err)
        """

        return self._calc_stats_engines(pert_val, engines, statistic="MAE")

    def calc_mad_engines(self, pert_val: str = None, engines: Optional[list] = None):
        """calculate the Mean Absolute Deviation (MAD) for between all the engines.

        Args:
            pert_val (str, optional): whether plotting 'pert' ie perturbations or 'val' ie values (per ligand result). Defaults to None.
            engines (list, optional): names of engines / other results names to calculate for.

        Returns:
            tuple: of dataframe of value and error (df, df_err)
        """

        return self._calc_stats_engines(pert_val, engines, statistic="MAD")

    def calc_rmse_engines(self, pert_val: str = None, engines: Optional[list] = None):
        """calculate the Root Mean Squared Error (RMSE) for between all the engines and the experimental.

        Args:
            pert_val (str, optional): whether plotting 'pert' ie perturbations or 'val' ie values (per ligand result). Defaults to None.
            engines (list, optional): names of engines / other results names to calculate for.

        Returns:
            tuple: of dataframe of value and error (df, df_err)
        """

        return self._calc_stats_engines(pert_val, engines, statistic="RMSE")

    def calc_spearmans_rank_engines(
        self, pert_val: str = None, engines: Optional[list] = None
    ):
        """calculate the Spearman's Rank correlation coefficient for between all the engines and the experimental.

        Args:
            pert_val (str, optional): whether plotting 'pert' ie perturbations or 'val' ie values (per ligand result). Defaults to None.
            engines (list, optional): names of engines / other results names to calculate for.

        Returns:
            tuple: of dataframe of value and error (df, df_err)
        """

        return self._calc_stats_engines(pert_val, engines, statistic="Spearmans")

    def _calc_stats_engines(
        self,
        pert_val: str = None,
        engines: Optional[list] = None,
        statistic: str = None,
    ) -> tuple:
        """internal function to wrap around the stats object and return dataframe.

        Args:
            pert_val (str, optional): whether for pert or val. Defaults to None.
            engines (Optional[list], optional): engines to calculate for. Defaults to None.
            statistic (str, optional): The statistic to calculate for. Defaults to None.

        Returns:
            tuple: (df, df_err)
        """

        self._initialise_stats_object(check=True)

        statistic_dict = {
            "MAD": self._stats_object.compute_mue,
            "MAE": self._stats_object.compute_mue,
            "RMSE": self._stats_object.compute_rmse,
            "Spearmans": self._stats_object.compute_rho,
        }

        func = statistic_dict[statistic]

        if engines:
            engines = self._validate_in_names_list(engines, make_list=True)
        else:
            engines = self.engines + self.other_results_names

        if statistic == "MAD":
            enginesb = engines
        else:
            enginesb = ["experimental"]

        pv = validate.pert_val(pert_val)

        df = pd.DataFrame(columns=engines, index=enginesb)
        df_err = pd.DataFrame(columns=engines, index=enginesb)

        # iterate compared to experimental
        for combo in it.product(engines, enginesb):
            eng1 = combo[0]
            eng2 = combo[1]

            values = func(pv, x=eng1, y=eng2)
            mean_absolute_error = values[0]  # the computed statistic
            err = values[1]  # the stderr from bootstrapping

            # loc index, column
            df.loc[eng2, eng1] = mean_absolute_error
            df_err.loc[eng2, eng1] = err

        df.to_csv(
            f"{self.results_folder}/{statistic}_{pert_val}_{self.file_ext}.csv", sep=" "
        )
        df_err.to_csv(
            f"{self.results_folder}/{statistic}_err_{pert_val}_{self.file_ext}.csv",
            sep=" ",
        )

        return df, df_err

    # freenergworkflows stuff for comparison
    def _add_fwf_path(self, fwf_path: str):
        # using freenergworkflows
        if not fwf_path:
            raise ValueError("pls incl the path to freenergworkflows")

        fwf_path = validate.folder_path(fwf_path)

        if fwf_path not in sys.path:
            sys.path.insert(1, fwf_path)

        self._fwf_path = fwf_path

    def _get_exp_fwf(self) -> tuple:
        """get experimental values using freenergworkflows

        Returns:
            tuple: (exp_lig_dict, exp_pert_dict)
        """
        if not self._fwf_path:
            raise ValueError("need fwf path added using _add_fwf_path(fwf_path)")
        import experiments

        # first need to convert the yml file into one useable by freenergworkflows
        exp_file_dat = f"{self.exp_file.split('.')[0]}_exp_dat.dat"
        convert.yml_into_freenrgworkflows(self.exp_file, exp_file_dat)
        experiments = experiments.ExperimentalData()
        experiments.compute_affinities(
            exp_file_dat, data_type="IC50", comments="#", delimiter=","
        )
        experimental_DDGs = experiments.freeEnergiesInKcal

        exp_pert_dict, exp_lig_dict = make_dict.experimental_from_freenrgworkflows(
            experimental_DDGs, self.ligands, self.perturbations
        )
        self._fwf_experimental_DDGs = experimental_DDGs

        return exp_lig_dict, exp_pert_dict

    def _get_ana_fwf(self, engine: str = None) -> dict:
        """get experimental values using freenergworkflows

        Args:
            engine (str): name of engine. Defaults to None.

        Raises:
            ValueError: need an engine

        Returns:
            dict: freenerg dict of results for that engine
        """
        # using freenergworkflows
        if not self._fwf_path:
            raise ValueError("need fwf path added using _add_fwf_path(fwf_path)")
        import networkanalysis

        if not engine:
            raise ValueError("please incl an engine")

        # using the network analyser
        nA = networkanalysis.NetworkAnalyser()

        first_file = False
        nf = 0
        for file_name in self._results_repeat_files[engine]:
            # rewrite the file to include only lig_0, lig_1, freenrg, error, engine
            new_file_name = f"{self.results_folder}/fwf_{engine}_file_{nf}.csv"
            data = pd.read_csv(file_name, delimiter=",")
            header_data = data[["lig_0", "lig_1", "freenrg", "error", "engine"]]
            clean_data = header_data.replace("kcal/mol", "", regex=True)
            pd.DataFrame.to_csv(clean_data, new_file_name, sep=",", index=False)
            nf += 1

            if first_file is False:
                nA.read_perturbations_pandas(
                    new_file_name, comments="#", source="lig_0", target="lig_1"
                )
                first_file = True
            else:
                # add more replicates to the graph. FreeNrgWorkflows will take care of averaging
                # the free energies as well as propagating the error.
                nA.add_data_to_graph_pandas(
                    new_file_name, comments="#", source="lig_0", target="lig_1"
                )

        # set network analyser graph as graph
        self.fwf_graph = nA._graph

        computed_relative_DDGs = nA.freeEnergyInKcal

        # this is the per ligand results
        freenrg_dict = make_dict.from_freenrgworkflows_network_analyser(
            computed_relative_DDGs
        )
        self._fwf_computed_relative_DDGs.update({engine: computed_relative_DDGs})

        return freenrg_dict

    def _get_stats_fwf(self, engine: str = None) -> tuple:
        """get stats using freenergworkflows for the ligands

        Args:
            engine (str): name of engine. Defaults to None.

        Raises:
            ValueError: need an engine

        Returns:
            tuple: r_confidence, tau_confidence, mue_confidence
        """
        # using freenergworkflows
        if not self._fwf_path:
            raise ValueError("need fwf path added using _add_fwf_path(fwf_path)")
        import stats

        if not engine:
            raise ValueError("please incl an engine")

        # these are the per ligand results
        computed_relative_DDGs = self._fwf_computed_relative_DDGs[engine]
        self._get_exp_fwf()
        experimental_DDGs = self._fwf_experimental_DDGs

        _stats = stats.freeEnergyStats()
        _stats.generate_statistics(
            computed_relative_DDGs, experimental_DDGs, repeats=10000
        )
        r_confidence = _stats.R_confidence
        tau_confidence = _stats.tau_confidence
        mue_confidence = _stats.mue_confidence
        logging.info(
            "R confidence is:   %.2f < %.2f < %.2f"
            % (r_confidence[1], r_confidence[0], r_confidence[2])
        )
        logging.info(
            "MUE confidence is: %.2f < %.2f < %.2f"
            % (mue_confidence[1], mue_confidence[0], mue_confidence[2])
        )
        logging.info(
            "Tau confidence is: %.2f < %.2f < %.2f"
            % (tau_confidence[1], tau_confidence[0], tau_confidence[2])
        )

        return r_confidence, tau_confidence, mue_confidence

    def _get_mad_fwf(self, enginesa: str, enginesb: str) -> tuple:
        mad_df = pd.DataFrame(columns=enginesa, index=enginesb)
        mad_df_err = pd.DataFrame(columns=enginesa, index=enginesb)

        for combo in it.product(enginesa, enginesb):
            eng1 = combo[0]
            eng2 = combo[1]

            freenrg_dict_eng1 = self._get_ana_fwf(eng1)
            freenrg_dict_eng2 = self._get_ana_fwf(eng2)

            new_freenrg_dict = {}

            for key in freenrg_dict_eng1:
                new_freenrg_dict[key] = (
                    freenrg_dict_eng1[key][0],  # value
                    freenrg_dict_eng1[key][1],  # error
                    freenrg_dict_eng2[key][0],
                    freenrg_dict_eng2[key][1],
                )

            df = pd.DataFrame.from_dict(
                new_freenrg_dict,
                columns=["eng1_value", "eng1_err", "eng2_value", "eng2_err"],
                orient="index",
            ).dropna()

            values = stats_engines.compute_stats(
                x=df["eng1_value"],
                y=df["eng2_value"],
                xerr=df["eng1_err"],
                yerr=df["eng2_err"],
                statistic="MUE",
            )
            mean_absolute_deviation = values[0]  # the computed statitic
            mad_err = values[1]

            # loc index, column
            mad_df.loc[eng2, eng1] = mean_absolute_deviation
            mad_df_err.loc[eng2, eng1] = mad_err

        mad_df.to_csv(f"{self.results_folder}/MAD_fwf_{self.file_ext}.csv", sep=" ")
        mad_df_err.to_csv(
            f"{self.results_folder}/MAD_err_fwf_{self.file_ext}.csv", sep=" "
        )

        return mad_df, mad_df_err

    def perturbing_atoms_and_overlap(
        self, prep_dir: str = None, outputs_dir: str = None, **kwargs
    ):
        if prep_dir:
            prep_dir = validate.folder_path(prep_dir)
            calc_atom_mappings = True
        else:
            logging.info(
                "please provide the prep dir to use for calculating atom mappings. Will only look at overlap."
            )
            calc_atom_mappings = False

        if outputs_dir:
            outputs_dir = validate.folder_path(outputs_dir)
        else:
            logging.error("please provide the dir where the outputs are located")
            return

        pert_dict = self.exper_pert_dict

        with open(f"{self.results_folder}/perturbing_overlap.dat", "w") as f:
            writer = csv.writer(f)
            writer.writerow(
                [
                    "perturbation",
                    "engine",
                    "perturbing_atoms",
                    "percen_overlap_okay",
                    "too_small_avg",
                    "diff_to_exp",
                    "error",
                ]
            )
            for pert, engine in it.product(self.perturbations, self.engines):
                logging.info(f"running {pert}, {engine}....")

                if calc_atom_mappings:
                    lig_0 = pert.split("~")[0]
                    lig_1 = pert.split("~")[1]

                    # Load equilibrated inputs for both ligands
                    system0 = BSS.IO.readMolecules(
                        [
                            f"{prep_dir}/{lig_0}_lig_equil_solv.rst7",
                            f"{prep_dir}/{lig_0}_lig_equil_solv.prm7",
                        ]
                    )
                    system1 = BSS.IO.readMolecules(
                        [
                            f"{prep_dir}/{lig_1}_lig_equil_solv.rst7",
                            f"{prep_dir}/{lig_1}_lig_equil_solv.prm7",
                        ]
                    )

                    pert_atoms = pipeline.prep.merge.no_perturbing_atoms_average(
                        system0, system1, **kwargs
                    )

                else:
                    pert_atoms = None

                try:
                    ana_obj = pipeline.analysis.analyse(
                        f"{outputs_dir}/{engine}/{pert}",
                        analysis_prot=self.analysis_options,
                    )
                    avg, error, repeats_tuple_list = ana_obj.analyse_all_repeats()
                    percen_okay, too_smalls_avg = ana_obj.check_overlap()
                    diff = abs(pert_dict[pert][0] - avg.value())
                    err = error.value()
                except Exception as e:
                    logging.error(e)
                    percen_okay = None
                    too_smalls_avg = None
                    diff = None
                    err = None

                row = [pert, engine, pert_atoms, percen_okay, too_smalls_avg, diff, err]
                writer.writerow(row)
