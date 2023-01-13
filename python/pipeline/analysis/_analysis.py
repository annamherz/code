# import libraries
import warnings
import BioSimSpace as BSS
from BioSimSpace import Units as _Units
import pickle
import os as _os
import numpy as _np
import math as _math
from scipy.stats import sem
import pickle
import csv

from ..utils import *


class analyse():
    """class to analyse a work dir 
    """

    def __init__(self, work_dir):
        # instantiate the class with the work directory

        self._work_dir = validate.folder_path(work_dir)
        self._pickle_dir = validate.folder_path(f"{self._work_dir}/pickle", create=True)

        # get the perturbation name and engine from the folder path
        try:
            self.perturbation = self._work_dir.split("/")[-1]
            self.ligand_0 = self.perturbation.split("~")[0]
            self.ligand_1 = self.perturbation.split("~")[1]
            self.engine = validate.engine(
                self._work_dir.split("/")[-2].replace("_extracted", ""))
        except:
            warnings.warn(
                "was unable to get the perturbation name and engine from the file path.\n please add these to the class using self.perturbation and self.engine")

        # initialise dictionaries
        self._bound_pmf_dict = {}  # for the intial results
        self._free_pmf_dict = {}
        self._bound_matrix_dict = {}  # for the energy matrix
        self._free_matrix_dict = {}
        self._bound_val_dict = {}  # for the actual value for all the windows
        self._free_val_dict = {}
        self._bound_err_dict = {}  # for the final error defined
        self._free_err_dict = {}
        self.repeats_tuple_list = []

        # intialise other things
        self._get_repeat_folders()
        self._set_default_options()
        self._file_ext()
        self._pickle_ext()
        self.is_analysed = False


    def _set_default_options(self):

        self.estimator = "MBAR"
        self.method = "alchemlyb"
        self._check_overlap = True
        self._save_pickle = True
        self._try_pickle = True

        # for the preprocessing
        self._auto_equilibration = False
        self._statistical_inefficiency = False
        self._truncate_percentage = 0  # no truncation
        self._truncate_keep = "end"

    def _file_ext(self):

        file_ext = (f"{self.estimator}_{self.method}_"+
                    f"eq{str(self._auto_equilibration).lower()}_"+
                    f"stats{str(self._statistical_inefficiency).lower()}_"+
                    f"truncate{str(self._truncate_percentage)}{self._truncate_keep}")

        self.file_ext = file_ext

        return file_ext

    def _pickle_ext(self):

        if not self.file_ext:
            self._file_ext()

        pickle_ext = (f"{self.perturbation}_{self.engine}_"+
                      f"{self.file_ext}")

        self.pickle_ext = pickle_ext

        return pickle_ext


    def set_options(self, options_dict):

        options_dict = validate.dictionary(options_dict)

        if "estimator" in options_dict:
            self.estimator = validate.string(options_dict["estimator"])
            if self.estimator not in ['MBAR', 'TI']:
                raise ValueError("'estimator' must be either 'MBAR' or 'TI'.")
        if "check_overlap" in options_dict:
            self._check_overlap = validate.boolean(
                options_dict["check_overlap"])
        if self._check_overlap == "True" and self.estimator != "MBAR":
            self._check_overlap = False

        if "method" in options_dict:
            self.method = validate.string(options_dict["method"])
            if self.method not in ['alchemlyb', 'native']:
                raise ValueError(
                    "'estimator' must be either 'alchemlyb' or 'native'.")
        if "save_pickle" in options_dict:
            self._save_pickle = validate.boolean(options_dict["save_pickle"])
        if "try_pickle" in options_dict:
            self._try_pickle = validate.boolean(options_dict["try_pickle"])

        if "auto_equilibration" in options_dict:
            self._auto_equilibration = validate.boolean(
                options_dict["auto_equilibration"])
        if "statistical_inefficiency" in options_dict:
            self._statistical_inefficiency = validate.boolean(
                options_dict["statistical_inefficiency"])
        if "truncate_percentage" in options_dict:
            self._truncate_percentage = validate.integer(
                options_dict["truncate_percentage"])
        if "truncate_keep" in options_dict:
            self._truncate_keep = validate.string(
                options_dict["truncate_keep"])
            if self._truncate_keep not in ['start', 'end']:
                raise ValueError(
                    "'truncate_keep' must be either 'start' or 'end'.")

        # reset the file extensions
        self._file_ext()
        self._pickle_ext()


    def _get_repeat_folders(self):

        self._work_dir = self._work_dir
        # Read how many repeats are in the direct
        # ory.
        folders = (next(_os.walk(self._work_dir))[1])
        self._b_folders, self._f_folders = [], []
        for f in folders:
            if 'bound' in f:
                self._b_folders.append(f'{f}')
            elif 'free' in f:
                self._f_folders.append(f'{f}')
            else:
                continue

        # sort the folders
        self._b_folders.sort()
        self._f_folders.sort()

        if not self._b_folders:
            raise ValueError(
                "Couldn't find any folders with 'bound' in the specified directory?")
        elif not self._f_folders:
            raise ValueError(
                "Couldn't find any folders with 'free' in the specified directory?")

        no_of_b_repeats = len(self._b_folders)
        no_of_f_repeats = len(self._f_folders)
        self._b_repeats = list(range(no_of_b_repeats))
        self._f_repeats = list(range(no_of_f_repeats))

        if no_of_b_repeats != no_of_f_repeats:
            print(
                f"There are a different number of repeats for bound ({no_of_b_repeats}) and free ({no_of_f_repeats}) for {self._work_dir}."+
                f"these are {self._b_folders} and {self._f_folders}.")
        else:
            print(
                f"There are {no_of_b_repeats} repeats for each the bound and the free for {self._work_dir}."+
                f"these are {self._b_folders} and {self._f_folders}.")

        self._b_repeats = self._b_repeats
        self._f_repeats = self._f_repeats
        self._no_of_b_repeats = no_of_b_repeats
        self._no_of_f_repeats = no_of_f_repeats


    def _check_pickle(self):

        try_pickle = True
        pickle_ext = self.pickle_ext

        # try loading in if previously calculated
        if try_pickle:
            try:
                print(
                    f"trying to locate pickles in default pickle folder, {self._pickle_dir} for {pickle_ext}...")
                with open(f"{self._pickle_dir}/bound_pmf_{pickle_ext}.pickle", "rb") as file:
                    self._bound_pmf_dict = pickle.load(file)
                with open(f"{self._pickle_dir}/free_pmf_{pickle_ext}.pickle", "rb") as file:
                    self._free_pmf_dict = pickle.load(file)
                with open(f"{self._pickle_dir}/bound_matrix_{pickle_ext}.pickle", 'rb') as file:
                    self._bound_matrix_dict = pickle.load(file)
                with open(f"{self._pickle_dir}/free_matrix_{pickle_ext}.pickle", 'rb') as file:
                    self._free_matrix_dict = pickle.load(file)
                with open(f"{self._pickle_dir}/bound_val_{pickle_ext}.pickle", 'rb') as file:
                    self._bound_val_dict = pickle.load(file)
                with open(f"{self._pickle_dir}/free_val_{pickle_ext}.pickle", 'rb') as file:
                    self._free_val_dict = pickle.load(file)
                with open(f"{self._pickle_dir}/bound_err_{pickle_ext}.pickle", 'rb') as file:
                    self._bound_err_dict = pickle.load(file)
                with open(f"{self._pickle_dir}/free_err_{pickle_ext}.pickle", 'rb') as file:
                    self._free_err_dict = pickle.load(file)
                print("pickles found!")
            except:
                print("loading pickle failed. Calculating normally.")
                try_pickle = False

        return try_pickle


    def analyse_all_repeats(self):
        """Analyse all existing free-energy data from a simulation working directory.
        """

        if self._try_pickle:
            do_pickle = analyse._check_pickle(self)
        else:
            do_pickle = False
            
        if do_pickle:
            freenrg_rel, repeats_tuple_list = analyse._analyse_all_repeats_pickle(
                self)
        else:
            freenrg_rel, repeats_tuple_list = analyse._analyse_all_repeats_normal(
                self)

        # set the average and error
        self.freenrg = freenrg_rel[0]
        self.error = freenrg_rel[1]
        self.repeats_tuple_list = repeats_tuple_list
        self.is_analysed = True

        if self._check_overlap:
            analyse.check_overlap(self)
    
        if self._save_pickle:
            if do_pickle:
                print("already using pickles, will not be saving again.")
            else:
                print("saving the pmf dictionaries for bound and free as pickles.")
                analyse.save_pickle(self)

        return (freenrg_rel[0], freenrg_rel[1], repeats_tuple_list)


    def _analyse_all_repeats_normal(self):
        """Analyse all existing free-energy data from a simulation working directory.
        """

        # list for the successful calculations
        bound_calculated = []
        free_calculated = []

        process_dict = {"auto_equilibration": self._auto_equilibration,
                        "statistical_inefficiency": self._statistical_inefficiency,
                        "truncate_percentage": self._truncate_percentage,
                        "truncate_keep": self._truncate_keep}

        # Analyse the results for each leg of the transformation.
        for b in self._b_repeats:
            try:
                name = str(b) + '_bound'
                pmf_bound, overlap_matrix_bound = BSS.FreeEnergy.Relative.analyse(
                    f'{self._work_dir}/{self._b_folders[b]}', estimator=self.estimator, method=self.method, **process_dict)
                self._bound_pmf_dict.update({name: pmf_bound})
                self._bound_matrix_dict.update({name: overlap_matrix_bound})
                bound_calculated.append(name)
            except Exception as e:
                print(e)
                print(
                    f'Unable to analyse values for {name}, which is repeat {self._b_folders[b]} in {self._work_dir}.')

        for f in self._f_repeats:
            try:
                name = str(f) + '_free'
                pmf_free, overlap_matrix_free = BSS.FreeEnergy.Relative.analyse(
                    f'{self._work_dir}/{self._f_folders[f]}', estimator=self.estimator, method=self.method, **process_dict)
                self._free_pmf_dict.update({name: pmf_free})
                self._free_matrix_dict.update({name: overlap_matrix_free})
                free_calculated.append(name)
            except Exception as e:
                print(e)
                print(
                    f'Unable to analyse values for {name}, which is repeat {self._f_folders[f]} in {self._work_dir}.')

        # create a dictionary of the calculated values, for each bound and free
        for r in self._b_repeats:
            try:
                bound_name = str(r) + '_bound'
                bound_val = (self._bound_pmf_dict[bound_name])[-1][1] - \
                    (self._bound_pmf_dict[bound_name])[0][1]
                # TODO change this so as in diff?
                bound_err = (self._bound_pmf_dict[bound_name])[-1][2]
                self._bound_val_dict.update({bound_name: bound_val})
                self._bound_err_dict.update({bound_name: bound_err})
            except:
                print(f'''Unable to compute values for {bound_name} in {self._work_dir}.\
                    Check earlier error message if these values could be analysed.''')
        for r in self._f_repeats:
            try:
                free_name = str(r) + '_free'
                free_val = (self._free_pmf_dict[free_name])[-1][1] - \
                    (self._free_pmf_dict[free_name])[0][1]
                # TODO change this so as in diff?
                free_err = (self._free_pmf_dict[free_name])[-1][2]
                self._free_val_dict.update({free_name: free_val})
                self._free_err_dict.update({free_name: free_err})
            except:
                print(f'''Unable to compute values for {free_name} in {self._work_dir}.\
                    Check earlier error message if these values could be analysed.''')

        freenrg_rel, repeats_tuple_list = analyse._calculate_freenrg(
            self, free_calculated, bound_calculated)

        return (freenrg_rel, repeats_tuple_list)


    def _analyse_all_repeats_pickle(self):
        """Analyse all existing free-energy data from a simulation working directory.
        """

        # list for the successful calculations
        bound_calculated = []
        free_calculated = []

        # Analyse the results for each leg of the transformation.
        for b in self._b_repeats:
            try:
                name = str(b) + '_bound'
                if name in self._bound_pmf_dict.keys():
                    bound_calculated.append(name)
            except Exception as e:
                print(e)
                print(
                    f'Unable to analyse values for {name}, which is repeat {self._b_folders[b]} in {self._work_dir}.')

        for f in self._f_repeats:
            try:
                name = str(f) + '_free'
                if name in self._free_pmf_dict.keys():
                    free_calculated.append(name)
            except Exception as e:
                print(e)
                print(
                    f'Unable to analyse values for {name}, which is repeat {self._f_folders[f]} in {self._work_dir}.')

        freenrg_rel, repeats_tuple_list = analyse._calculate_freenrg(
            self, free_calculated, bound_calculated)

        return (freenrg_rel, repeats_tuple_list)


    def _calculate_freenrg(self, free_calculated, bound_calculated):

        # list for all the repeats
        repeats_tuple_list = []

        # get average of free energy values just calculated
        # if there is only one repeat use the BSS difference function
        if len(self._free_val_dict.values()) == 1 and len(self._bound_val_dict.values()) == 1:
            freenrg_rel = BSS.FreeEnergy.Relative.difference(
                list(self._bound_pmf_dict.items())[0][1], list(self._free_pmf_dict.items())[0][1])
            freenrg_val = freenrg_rel[0].value()
            freenrg_err = freenrg_rel[1].value()

        # otherwise, calculate the average and the SEM
        else:
            free_vals = list(
                val/_Units.Energy.kcal_per_mol for val in self._free_val_dict.values())
            free_avg = _np.mean(free_vals)
            free_sem = sem(free_vals)
            bound_vals = list(
                val/_Units.Energy.kcal_per_mol for val in self._bound_val_dict.values())
            bound_avg = _np.mean(bound_vals)
            bound_sem = sem(bound_vals)
            freenrg_val = (bound_avg-free_avg)
            freenrg_err = (_math.sqrt(
                _math.pow(bound_sem, 2)+_math.pow(free_sem, 2)))
            freenrg_rel = (freenrg_val * _Units.Energy.kcal_per_mol,
                           freenrg_err * _Units.Energy.kcal_per_mol)

        # create tuple list of each repeat that was calculated
        # first check the length of the calculated values and check if this is also the length of the folders found
        if len(bound_calculated) != self._no_of_b_repeats:
            print("the number of calculated values for bound does not match the number of bound folders.\n maybe try reanalysing / check errors?")
        if len(free_calculated) != self._no_of_f_repeats:
            print("the number of calculated values for free does not match the number of free folders.\n maybe try reanalysing / check errors?")

        # if the numebr of calculated values is the same, match these evenly
        if len(bound_calculated) == len(free_calculated):
            print(
                f"There are {len(bound_calculated)} calculated values for each the bound and the free leg for the folders in {self._work_dir}.")
            no_of_repeats = len(bound_calculated)
            repeats = list(range(no_of_repeats))
            for r in repeats:
                freenrg_rel = BSS.FreeEnergy.Relative.difference(
                    self._bound_pmf_dict[bound_calculated[r]], self._free_pmf_dict[free_calculated[r]])
                freenrg_val = freenrg_rel[0].value()
                freenrg_err = freenrg_rel[1].value()
                repeats_tuple_list.append(
                    (f"{str(r)}_repeat", freenrg_val, freenrg_err))

        elif len(bound_calculated) != len(free_calculated):
            print(f"There are {len(bound_calculated)} calculated values for the bound and \
            {len(free_calculated)} calculated values for the free leg for the folders in {self._work_dir}.")
            # use the shorter calculated values as the number of complete repeats
            if len(bound_calculated) < len(free_calculated):
                no_of_repeats = len(bound_calculated)
            else:
                no_of_repeats = len(free_calculated)
            print(
                f"The number of calculated values do not match. {no_of_repeats} repeats will be calculated.")
            r = 0
            for b, f in zip(bound_calculated, free_calculated):
                print(f"calculating repeat {r} as {b} and {f}.")
                freenrg_rel = BSS.FreeEnergy.Relative.difference(
                    self._bound_pmf_dict[b], self._free_pmf_dict[f])
                freenrg_val = freenrg_rel[0].value()
                freenrg_err = freenrg_rel[1].value()
                repeats_tuple_list.append(
                    (f"{str(r)}_repeat", freenrg_val, freenrg_err))
                r += 1

        return freenrg_rel, repeats_tuple_list


    def save_pickle(self):

        self._pickle_dir = validate.folder_path(self._pickle_dir, create=True)
        pickle_ext = self.pickle_ext

        if not self.is_analysed:
            warnings.warn(
                "can't save pickle, not all repeats have been analysed. please self.analyse_all_repeats() first!")
        else:
            try:
                # write the pmf as a pickle
                with open(f"{self._pickle_dir}/bound_pmf_{pickle_ext}.pickle", 'wb') as handle:
                    pickle.dump(self._bound_pmf_dict, handle)
                with open(f"{self._pickle_dir}/free_pmf_{pickle_ext}.pickle", 'wb') as handle:
                    pickle.dump(self._free_pmf_dict, handle)
                with open(f"{self._pickle_dir}/bound_matrix_{pickle_ext}.pickle", 'wb') as handle:
                    pickle.dump(self._bound_matrix_dict, handle)
                with open(f"{self._pickle_dir}/free_matrix_{pickle_ext}.pickle", 'wb') as handle:
                    pickle.dump(self._free_matrix_dict, handle)
                with open(f"{self._pickle_dir}/bound_val_{pickle_ext}.pickle", 'wb') as handle:
                    pickle.dump(self._bound_val_dict, handle)
                with open(f"{self._pickle_dir}/free_val_{pickle_ext}.pickle", 'wb') as handle:
                    pickle.dump(self._free_val_dict, handle)
                with open(f"{self._pickle_dir}/bound_err_{pickle_ext}.pickle", 'wb') as handle:
                    pickle.dump(self._bound_err_dict, handle)
                with open(f"{self._pickle_dir}/free_err_{pickle_ext}.pickle", 'wb') as handle:
                    pickle.dump(self._free_err_dict, handle)
            except Exception as e:
                print(e)
                print("could not save pickle :( ")


    def check_overlap(self):

        if not self.is_analysed:
            warnings.warn(
                "can't check overlap, not all repeats have been analysed. please self.analyse_all_repeats() first!")

        else:
            # check overlap matrix if okay
            for b in self._b_repeats:
                try:
                    name = str(b) + '_bound'
                    overlap = self._bound_matrix_dict[name]
                    overlap_okay = BSS.FreeEnergy.Relative.checkOverlap(
                        overlap, estimator=self.estimator)
                except Exception as e:
                    print(e)
                    print(f"could not check overlap matrix for {name}")

            for f in self._f_repeats:
                try:
                    name = str(f) + '_free'
                    overlap = self._free_matrix_dict[name]
                    overlap_okay = BSS.FreeEnergy.Relative.checkOverlap(
                        overlap, estimator=self.estimator)
                except Exception as e:
                    print(e)
                    print(f"could not check overlap matrix for {name}")


    def plot_graphs(self):

        if not self.is_analysed:
            warnings.warn(
                "can't plot, not all repeats have been analysed. please self.analyse_all_repeats() first!")

        else:
            graph_dir = validate.folder_path(self._work_dir + '/graphs')

            if self.estimator == "MBAR":

                for b in self._b_repeats:
                    try:
                        name = str(b) + '_bound'
                        overlap = self._bound_matrix_dict[name]
                        ax = BSS.FreeEnergy.Relative.plot(
                            overlap, work_dir=graph_dir, file_name=f"{name}_overlap_MBAR_{self.pickle_ext}")
                    except Exception as e:
                        print(e)
                        print(f"could not plt overlap matrix for {name}")

                for f in self._f_repeats:
                    try:
                        name = str(f) + '_free'
                        overlap = self._free_matrix_dict[name]
                        ax = BSS.FreeEnergy.Relative.plot(
                            overlap, work_dir=graph_dir, file_name=f"{name}_overlap_MBAR_{self.pickle_ext}")
                    except Exception as e:
                        print(e)
                        print(f"could not plt overlap matrix for {name}")

            elif self.estimator == "TI":

                for b in self._b_repeats:
                    name = str(b) + '_bound'
                    overlap = self._bound_matrix_dict[name]
                    try:
                        ax = BSS.FreeEnergy.Relative.plot(
                            overlap, work_dir=graph_dir, file_name=f"{name}_dHdl_TI_{self.pickle_ext}")
                    except Exception as e:
                        print(e)
                        print(f"could not plt dhdl for {name}")

                for f in self._f_repeats:
                    name = str(f) + '_free'
                    overlap = self._free_matrix_dict[name]
                    try:
                        ax = BSS.FreeEnergy.Relative.plot(
                            overlap, work_dir=graph_dir, file_name=f"{name}_dHdl_TI_{self.pickle_ext}")
                    except Exception as e:
                        print(e)
                        print(f"could not plt dhdl for {name}")
