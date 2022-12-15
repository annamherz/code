import csv
import os

from ._validate import *

def write_analysis_file(analysis, results_dir):

    analysis = validate.analysis(analysis, analysed=True)
    results_dir = validate.folder_path(results_dir, create = True)

    # data point for average
    data_point_avg = [analysis.ligand_0, analysis.ligand_0,
                      analysis.freenrg, analysis.error,
                      analysis.engine, analysis.pickle_ext]

    # use csv to open the results file.
    with open(f"{results_dir}/final_summary_{analysis.pickle_ext}.csv", "a") as freenrg_writefile:
        writer = csv.writer(freenrg_writefile)

        # first, write a header if the file is created for the first time.
        if os.path.getsize(f"{results_dir}/final_summary_{analysis.pickle_ext}.csv") == 0:
            print(f"Starting {results_dir}/final_summary_{analysis.pickle_ext}.csv file.")
            writer.writerow(["lig_1", "lig_2", "freenrg",
                            "error", "engine", "method"])


    with open(f"{results_dir}/final_summary_{analysis.pickle_ext}.csv", "r") as freenrg_readfile:
        # then, grab all of the data that is already in the file.
        reader = csv.reader(freenrg_readfile)
        data_entries = [row for row in reader]

    # check if our data entry is not already in the results file. Raise an error if is.
    if data_point_avg in data_entries:
        warnings.warn(
            f"Results for in {analysis.perturbation}, {analysis.engine} "+
            f"are already in {results_dir}/final_summary_{analysis.pickle_ext}.csv .")

    else:
        # at this point we know that we are writing a new entry in the results file. Append the line to the file.
        # use csv to open the results file.
        with open(f"{results_dir}/final_summary_{analysis.pickle_ext}.csv", "a") as freenrg_writefile:
            writer = csv.writer(freenrg_writefile)
            print(
                f"Writing results. Average free energy of binding is {analysis.freenrg}"+
                f"and the error is {analysis.error} for {analysis.perturbation}, {analysis.engine}.")
            writer.writerow(data_point_avg)


    # write results for each repeat also
    no_repeats = list(range(len(analysis.repeats_tuple_list)))
    # use csv to open the results file.
    for r in no_repeats:
        data_point = [analysis.ligand_0,
                      analysis.ligand_1,
                      analysis.repeats_tuple_list[r][1],
                      analysis.repeats_tuple_list[r][2],
                      analysis.engine,
                      analysis.pickle_ext
                      ]
        results_file_path = f"{results_dir}/repeat_{no_repeats.index(r)}_{analysis.pickle_ext}.csv"
        with open(results_file_path, "a") as freenrg_writefile:
            writer = csv.writer(freenrg_writefile)

            # first, write a header if the file is created for the first time.
            if os.path.getsize(results_file_path) == 0:
                print(f"Starting {results_file_path} file.")
                writer.writerow(["lig_1", "lig_2", "freenrg",
                                "error", "engine", "estimator", "method"])


        with open(results_file_path, "r") as freenrg_readfile:
            # then, grab all of the data that is already in the file.
            reader = csv.reader(freenrg_readfile)
            data_entries = [row for row in reader]

        # check if our data entry is not already in the results file. Raise an error if is.
        if data_point in data_entries:
            warnings.warn(
                f"Results for in {analysis.perturbation}, {analysis.engine} are already in {results_file_path}.")

        else:
            # at this point we know that we are writing a new entry in the results file. Append the line to the file.
            # use csv to open the results file.
            with open(results_file_path, "a") as freenrg_writefile:
                writer = csv.writer(freenrg_writefile)
                print(
                    f"Writing results. For repeat {r}, free energy of binding is "+
                    f"{analysis.repeats_tuple_list[r][1]} and the error is {analysis.repeats_tuple_list[r][2]} "+
                    f"for {analysis.perturbation}, {analysis.engine}.")
                writer.writerow(data_point)
