import csv
import os

from ._validate import *

def write_analysis_file(analysis, results_dir):

    analysis = validate.analysis(analysis, analysed=True)
    results_dir = validate.folder_path(results_dir, create=True)

    # data point for average
    data_point_avg = [analysis.ligand_0, analysis.ligand_1,
                      analysis.freenrg, analysis.error,
                      analysis.engine, analysis.file_ext]

    # use csv to open the results file.
    final_summary_file = f"{results_dir}/final_summary_{analysis.engine.upper()}_{analysis.file_ext}.csv"
    with open(final_summary_file, "a") as freenrg_writefile:
        writer = csv.writer(freenrg_writefile)

        # first, write a header if the file is created for the first time.
        if os.path.getsize(final_summary_file) == 0:
            print(f"Starting {final_summary_file} file.")
            writer.writerow(["lig_0", "lig_1", "freenrg",
                            "error", "engine", "method"])


    with open(final_summary_file, "r") as freenrg_readfile:
        # then, grab all of the data that is already in the file.
        reader = csv.reader(freenrg_readfile)
        data_entries = [row for row in reader]

    # check if our data entry is not already in the results file. Raise an error if is.
    if data_point_avg in data_entries:
        warnings.warn(
            f"Results for in {analysis.perturbation}, {analysis.engine} "+
            f"are already in {final_summary_file} .")

    else:
        # at this point we know that we are writing a new entry in the results file. Append the line to the file.
        # use csv to open the results file.
        with open(final_summary_file, "a") as freenrg_writefile:
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
                      analysis.file_ext
                      ]
        results_file_path = f"{results_dir}/freenrg_repeat_{no_repeats.index(r)}_{analysis.engine.upper()}_{analysis.file_ext}.csv"
        with open(results_file_path, "a") as freenrg_writefile:
            writer = csv.writer(freenrg_writefile)

            # first, write a header if the file is created for the first time.
            if os.path.getsize(results_file_path) == 0:
                print(f"Starting {results_file_path} file.")
                writer.writerow(["lig_0", "lig_1", "freenrg",
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
    
    # write results for each the bound and the free too
    for bf in ["bound", "free"]:

        if bf == "free":
            val_dict = analysis._free_val_dict
            err_dict = analysis._free_err_dict
            name = bf
        elif bf == "bound":
            val_dict = analysis._bound_val_dict
            err_dict = analysis._bound_err_dict
            name = bf

        no_repeats = list(range(len(analysis.repeats_tuple_list)))
        # use csv to open the results file.
        for r in no_repeats:
            data_point = [analysis.ligand_0,
                        analysis.ligand_1,
                        val_dict[f"{r}_{name}"],
                        err_dict[f"{r}_{name}"],
                        analysis.engine,
                        analysis.file_ext
                        ]
            results_file_path = f"{results_dir}/{name}_repeat_{no_repeats.index(r)}_{analysis.engine.upper()}_{analysis.file_ext}.csv"
            with open(results_file_path, "a") as freenrg_writefile:
                writer = csv.writer(freenrg_writefile)

                # first, write a header if the file is created for the first time.
                if os.path.getsize(results_file_path) == 0:
                    print(f"Starting {results_file_path} file.")
                    writer.writerow(["lig_0", "lig_1", "freenrg",
                                    "error", "engine", "estimator", "method"])


            with open(results_file_path, "r") as freenrg_readfile:
                # then, grab all of the data that is already in the file.
                reader = csv.reader(freenrg_readfile)
                data_entries = [row for row in reader]

            # check if our data entry is not already in the results file. Raise an error if is.
            if data_point in data_entries:
                warnings.warn(
                    f"Results for {name} in {analysis.perturbation}, {analysis.engine} are already in {results_file_path}.")

            else:
                # at this point we know that we are writing a new entry in the results file. Append the line to the file.
                # use csv to open the results file.
                with open(results_file_path, "a") as freenrg_writefile:
                    writer = csv.writer(freenrg_writefile)
                    print(
                        f"Writing results. For repeat {name} {r}, free energy of binding is "+
                        f"{val_dict[f'{r}_{name}']} and the error is {err_dict[f'{r}_{name}']} "+
                        f"for {analysis.perturbation}, {analysis.engine}.")
                    writer.writerow(data_point)


def write_atom_mappings(lig_0, lig_1, ligand_0, ligand_1, mapping, output_file):

    # data point for average
    data_point = [lig_0, lig_1, ligand_0, ligand_1, mapping]

    # use csv to open the results file.
    with open(output_file, "a") as writefile:
        writer = csv.writer(writefile, delimiter=";")

        # first, write a header if the file is created for the first time.
        if os.path.getsize(output_file) == 0:
            print(f"Starting {output_file} file.")
            writer.writerow(["lig_0", "lig_1", "lig_0_atoms", "lig_1_atoms", "mapping"])


    with open(output_file, "r") as readfile:
        # then, grab all of the data that is already in the file.
        reader = csv.reader(readfile)
        data_entries = [row for row in reader]

    # check if our data entry is not already in the results file. Raise an error if is.
    if data_point in data_entries:
        warnings.warn(
            f"this atom mapping ahs already been written.")

    else:
        # at this point we know that we are writing a new entry in the results file. Append the line to the file.
        # use csv to open the results file.
        with open(output_file, "a") as writefile:
            writer = csv.writer(writefile, delimiter=";")
            print(
                f"Writing results.")
            writer.writerow(data_point)


def write_modified_results_files(results_files, perturbations, output_folder=None, extra_options=None):

    len_results_files = 0
    for file in results_files:
        validate.file_path(file)
        len_results_files += 1
    # print(f"there are : {len_results_files} results files.")

    if not output_folder:
        # will write as folder of first results file
        output_folder = validate.folder_path(results_files[0].replace(results_files[0].split("/")[-1], "")[:-1])
        print(f"using {output_folder} to write the results as none specified...")
    
    # set extra_options variables as defaults
    engines = [eng.upper() for eng in BSS.FreeEnergy.engines()] # use all

    if extra_options:
        extra_options = validate.dictionary(extra_options)

        if "engine" in extra_options.keys():
            engines = [validate.engine(extra_options["engine"])]
        if "engines" in extra_options.keys():
            engines = validate.is_list(extra_options["engines"])
            for engine in engines:
                engine_val = validate.engine(engine)
                engines = [engine_val if i == engine else i for i in engines]
       
    mod_results_files = []

    for file in results_files:
        new_file_name = f"{output_folder}/results_{results_files.index(file)}_{'_'.join(engines)}.csv"
        with open(new_file_name, "w") as result_file:

            writer = csv.writer(result_file, delimiter=",")
            writer.writerow(["lig_1","lig_2","freenrg","error","engine"])

            for row, index in pd.read_csv(file).iterrows():
                pert = f"{index['lig_1']}~{index['lig_2']}"
                if pert in perturbations and index['engine'].strip() in engines:
                        writer.writerow([index['lig_1'], index['lig_2'], index['freenrg'], index['error'], index['engine']])    

            mod_results_files.append(new_file_name)