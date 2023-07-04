# import libraries
import BioSimSpace as BSS

import yaml
import csv
import numpy as np

from ..utils import *
from ._dictionaries import *

# conversion of constants
from scipy.constants import R, calorie

kJ2kcal = 1 / calorie
R_kJmol = R / 1000
R_kcalmol = R_kJmol * kJ2kcal


class convert:
    """class of static methods for converting data"""

    def __init__(self):
        pass

    # TODO more robust yml file conversion
    @staticmethod
    def yml_into_freenrgworkflows(exp_file, exp_file_dat):
        """convert yml format into one suitable for freenergworkflows

        Args:
            exp_file (str): yml file path of experimental data
            exp_file_dat (str): new file to write experimental data to (fwf format)
        """
        # get the experimental data into a useable format (from yml to csv)
        # for freenergworkflows, want to save as lig, Ki
        # experimental values (e.g. ic50/ki) for all ligands in our set.

        exp_file = validate.file_path(exp_file)
        exp_file_dat = validate.string(exp_file_dat)

        with open(exp_file, "r") as file:
            data = yaml.safe_load(file)  # loads as dictionary

        with open(exp_file_dat, "w") as file:
            writer = csv.writer(file, delimiter=",")
            writer.writerow(["ligand", "value", "error"])

            # the data needs to be IC50, uM
            # am assuming that ki and IC50 are the same

            for key in data.keys():  # write for each ligand that was in yaml file
                if data[key]["measurement"]["unit"] == "uM":
                    writer.writerow(
                        [
                            key,
                            data[key]["measurement"]["value"],
                            data[key]["measurement"]["error"],
                        ]
                    )
                elif data[key]["measurement"]["unit"] == "nM":
                    writer.writerow(
                        [
                            key,
                            "{:.4f}".format(data[key]["measurement"]["value"] / 1000),
                            data[key]["measurement"]["error"] / 1000,
                        ]
                    )

    @staticmethod
    def convert_M_kcal(value, magnitude="uM", temp=300):
        """convert value into kcal/mol

        Args:
            value (float): value
            magnitude (str, optional): unit/magnitude of the value. uM or nM . Defaults to "uM".
            temp (int, optional): temperature of the simulation. Defaults to 300.

        Returns:
            string: value in kcal/mol
        """
        if magnitude == "uM":
            power = 10**-6
        if magnitude == "nM":
            power = 10**-9
        # gas constant in kcal per Kelvin per mol, exp val converted into M
        kcal_val = R_kcalmol * temp * np.log(value * (power))

        return kcal_val

    @staticmethod
    def yml_into_exper_dict(exp_file, temp=300):
        """convert yml file into experimental dictionary of values.

        Args:
            exp_file (str): yml file
            temp (int, optional): Temperature to use during the conversion. Defaults to 300.

        Returns:
            dict: kcal/mol value for each ligand
        """

        exp_file = validate.file_path(exp_file)
        temp = validate.is_float(temp)

        # TODO different data formats
        with open(exp_file, "r") as file:
            data = yaml.safe_load(file)  # loads as dictionary

        exper_raw_dict = {}
        for key in data.keys():  # write for each ligand that was in yaml file
            if data[key]["measurement"]["unit"] == "uM":
                exper_raw_dict[key] = (
                    data[key]["measurement"]["value"],
                    data[key]["measurement"]["error"],
                )
            elif data[key]["measurement"]["unit"] == "nM":
                exper_raw_dict[key] = (
                    "{:.4f}".format(data[key]["measurement"]["value"] / 1000),
                    data[key]["measurement"]["error"] / 1000,
                )

        exper_val_dict = {}

        for key in exper_raw_dict.keys():
            lig = str(key)

            exp_val = float(exper_raw_dict[key][0])
            # convert into kcal mol
            exp_kcal = convert.convert_M_kcal(exp_val, temp=temp)

            # convert both upper and lower error bounds for this too
            # get average and keep this as the error
            err = float(exper_raw_dict[key][1])
            exp_upper = exp_val + err
            exp_lower = exp_val - err
            exp_upper_kcal = convert.convert_M_kcal(exp_upper, temp=temp)
            exp_lower_kcal = convert.convert_M_kcal(exp_lower, temp=temp)
            err_kcal = abs(exp_upper_kcal - exp_lower_kcal) / 2

            # add to dict
            exper_val_dict.update({lig: (exp_kcal, err_kcal)})

        return exper_val_dict

    @staticmethod
    def cinnabar_file(
        results_files, exper_val, output_file, perturbations=None, method=None
    ):
        """convert results files into format needed for cinnabar. If multiple results files, uses the average of a perturbation.

        Args:
            results_files (list): list of results files
            exper_val (dict or str): dict of experimental values or yml file of experimental values.
            output_file (str): output file path
            perturbations (list, optional): list of perturbations to include. Defaults to None.
            method (str, optional): name of the method to consider.
        """
        # files is a list of files
        results_files = validate.is_list(results_files, make_list=True)
        # output file

        if exper_val:
            # first check if the experimental values are a dict or a file
            try:
                exper_val_dict = validate.dictionary(exper_val)
            except:
                validate.file_path(exper_val)
                print("input is a file, will convert this into a dict...")
                print("please check that the conversion of values is okay.")
                exper_val_dict = convert.yml_into_exper_dict(exper_val)

            add_exper_values = True

        else:
            add_exper_values = False

        # write to a csv file
        with open(f"{output_file}.csv", "w") as cinnabar_data_file:
            writer = csv.writer(cinnabar_data_file, delimiter=",")

            if add_exper_values:
                # first, write the experimental data
                writer.writerow(["# Experimental block"])
                writer.writerow(["# Ligand", "expt_DDG", "expt_dDDG"])

                # TODO calc exp from yml, follwed by conversion into dict
                for lig in exper_val_dict.keys():
                    writer.writerow(
                        [lig, f"{exper_val_dict[lig][0]}", f"{exper_val_dict[lig][1]}"]
                    )

            # second write the perturbation data
            writer.writerow([" "])
            writer.writerow(["# Calculated block"])
            writer.writerow(
                [
                    "# Ligand1",
                    "Ligand2",
                    "calc_DDG",
                    "calc_dDDG(MBAR)",
                    "calc_dDDG(additional)",
                ]
            )

            # need to write the average of the data, otherwise cinnabar just uses the last entry
            comp_diff_dict = make_dict.comp_results(
                results_files, perturbations=perturbations, method=method,
            )

            # write to file
            for key in comp_diff_dict:
                lig_0 = key.split("~")[0]
                lig_1 = key.split("~")[1]
                comp_ddG = comp_diff_dict[key][0]
                comp_err = comp_diff_dict[key][1]

                if not comp_ddG:
                    pass
                else:
                    if perturbations:
                        pert = f"{lig_0}~{lig_1}"
                        if pert in perturbations:
                            writer.writerow([lig_0, lig_1, comp_ddG, comp_err, "0.0"])
                        else:
                            pass

                    else:
                        writer.writerow([lig_0, lig_1, comp_ddG, comp_err, "0.0"])
