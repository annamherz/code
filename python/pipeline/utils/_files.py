import csv
import os
import numpy as np
import pandas as pd

from ._validate import *

def write_analysis_file(analysis, results_dir, method=None):
    """write the analysis file for the analysis object

    Args:
        analysis (pipeline.analysis.analyse): the analysed object from a AFE run with the pipeline
        results_dir (str): folder path to the results directory
        method (str, optional): method used, for description in the written file. Defaults to None.
    """
    analysis = validate.analysis(analysis, analysed=True)
    results_dir = validate.folder_path(results_dir, create=True)

    if not method:
        method = "None"

    # data point for average
    data_point_avg = [analysis.ligand_0,
                      analysis.ligand_1,
                      str(analysis.freenrg), # need as string so can compare to existing entries sometimes
                      str(analysis.error),
                      analysis.engine,
                      analysis.file_ext,
                      method
                      ]

    # use csv to open the results file.
    final_summary_file = f"{results_dir}/final_summary_{analysis.engine.upper()}_{analysis.file_ext}.csv"
    with open(final_summary_file, "a") as freenrg_writefile:
        writer = csv.writer(freenrg_writefile)

        # first, write a header if the file is created for the first time.
        if os.path.getsize(final_summary_file) == 0:
            print(f"Starting {final_summary_file} file.")
            writer.writerow(["lig_0", "lig_1", "freenrg",
                            "error", "engine", "analysis", "method"])


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
                f"Writing results. Average free energy of binding is {str(analysis.freenrg)} "+
                f"and the error is {str(analysis.error)} for {analysis.perturbation}, {analysis.engine}.")
            writer.writerow(data_point_avg)


    # write results for each calculated repeat
    no_repeats = list(range(len(analysis.repeats_tuple_list)))
    # use csv to open the results file.
    for r in no_repeats:
        data_point = [analysis.ligand_0,
                      analysis.ligand_1,
                      str(analysis.repeats_tuple_list[r][1]),
                      str(analysis.repeats_tuple_list[r][2]),
                      analysis.engine,
                      analysis.file_ext,
                      method
                      ]
        results_file_path = f"{results_dir}/freenrg_repeat_{no_repeats.index(r)}_{analysis.engine.upper()}_{analysis.file_ext}.csv"
        with open(results_file_path, "a") as freenrg_writefile:
            writer = csv.writer(freenrg_writefile)

            # first, write a header if the file is created for the first time.
            if os.path.getsize(results_file_path) == 0:
                print(f"Starting {results_file_path} file.")
                writer.writerow(["lig_0", "lig_1", "freenrg",
                                "error", "engine", "analysis", "method"])


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
                    f"{str(analysis.repeats_tuple_list[r][1])} and the error is {str(analysis.repeats_tuple_list[r][2])} "+
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

        # use csv to open the results file.
        for key in val_dict.keys():

            r= key.split('_')[0]
            
            try:
                data_point = [analysis.ligand_0,
                            analysis.ligand_1,
                            str(val_dict[f"{key}"]),
                            str(err_dict[f"{key}"]),
                            analysis.engine,
                            analysis.file_ext,
                            method
                            ]
            except:
                # if that repeat does not have a value, write none in the output file.
                data_point = [analysis.ligand_0,
                            analysis.ligand_1,
                            str(np.nan),
                            str(np.nan),
                            analysis.engine,
                            analysis.file_ext,
                            method
                            ]                
            results_file_path = f"{results_dir}/{name}_repeat_{r}_{analysis.engine.upper()}_{analysis.file_ext}.csv"
            with open(results_file_path, "a") as freenrg_writefile:
                writer = csv.writer(freenrg_writefile)

                # first, write a header if the file is created for the first time.
                if os.path.getsize(results_file_path) == 0:
                    print(f"Starting {results_file_path} file.")
                    writer.writerow(["lig_0", "lig_1", "freenrg",
                                    "error", "engine", "analysis", "method"])


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
                        f"{str(val_dict[f'{r}_{name}'])} and the error is {str(err_dict[f'{r}_{name}'])} "+
                        f"for {analysis.perturbation}, {analysis.engine}.")
                    writer.writerow(data_point)


def write_atom_mappings(lig_0, lig_1, ligand_0, ligand_1, mapping, output_file):
    """write the atom mappings previously found using merge.atom_mappings() to a file.

    Args:
        lig_0 (str): name of ligand 0
        lig_1 (str): name of ligand 1
        ligand_0 (): result of ligand_0.getAtoms()
        ligand_1 (): result of ligand_0.getAtoms()
        mapping (dict): mapping of the atoms
        output_file (str): output file path
    """
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
            f"this atom mapping has already been written.")

    else:
        # at this point we know that we are writing a new entry in the results file. Append the line to the file.
        # use csv to open the results file.
        with open(output_file, "a") as writefile:
            writer = csv.writer(writefile, delimiter=";")
            print(
                f"Writing results.")
            writer.writerow(data_point)


def write_modified_results_files(results_files, perturbations, output_folder=None, extra_options=None):
    """write modified results files that contain only specific perturbations.
    In the extra options, can specify 'engine' or 'engines'.

    Args:
        results_files (list): list of results files that are to be modified
        perturbations (list): list of perturbations to include
        output_folder (str, optional): output folder for new files. Defaults to None. Uses folder of passed results files.
        extra_options (dict, optional): dictionary of extra options, such as 'engine' or 'engines'. Defaults to None.

    Returns:
        list: list of the modified results files
    """

    results_files = validate.is_list(results_files, make_list=True)

    # find the length of the passed results files and validate this
    len_results_files = 0
    for file in results_files:
        validate.file_path(file)
        len_results_files += 1

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

    # create a list for the modified results files
    mod_results_files = []

    # write the new files
    for file in results_files:
        new_file_name = f"{output_folder}/results_{results_files.index(file)}_{'_'.join(engines)}.csv"
        with open(new_file_name, "w") as result_file:

            writer = csv.writer(result_file, delimiter=",")
            writer.writerow(["lig_0","lig_1","freenrg","error","engine","analysis","method"])

            # read the file as a csv, should have the same headings as well
            for row, index in pd.read_csv(file).iterrows():
                pert = f"{index['lig_0']}~{index['lig_1']}"
                if pert in perturbations and index['engine'].strip() in engines:
                        
                        # check if have analysis and method in it 
                        try:
                            ana_str = index['analysis']
                        except:
                            ana_str = "not specified"
                        
                        try:
                            method_str = index['method']
                        except:
                            method_str = "None"
                        
                        # write the row
                        writer.writerow([index['lig_0'],
                                         index['lig_1'], 
                                         index['freenrg'],
                                         index['error'],
                                         index['engine'],
                                         ana_str,
                                         method_str])    

            mod_results_files.append(new_file_name)
    
    return mod_results_files


def write_protocol(query_dict, file_path):
    """write the protocol dictionary as a file for the pipeline scripts.

    Args:
        query_dict (dict): dictionary of protocol options
        file_path (str): file path to where the file should be written.
    """

    file = validate.string(file_path)
    query_dict = validate.dictionary(query_dict)

    # write in the style needed for the protocol
    with open(file, "w") as protocol_file:
        writer = csv.writer(protocol_file)
        for query in query_dict.keys():
            writer.writerow([f"{query} = {query_dict[query]}"])


def write_ligands(ligand_names, file_path):
    """write the ligands to a file for the pipeline scripts.

    Args:
        ligand_names (list): list of ligand names
        file_path (str): file path to where the file should be written.
    """

    ligand_names = validate.is_list(ligand_names)
    file = validate.string(file_path)

    with open(file, "w") as ligands_file:
        writer = csv.writer(ligands_file)
        for lig in ligand_names:
            writer.writerow([lig])               

def write_lomap_scores(pert_network_dict, file_path):
    """write the lomap scores to a file.

    Args:
        pert_network_dict (dict): dictionary of perturbations with their LOMAP scores.
        file_path (str): file path to where the file should be written.
    """

    file = validate.string(file_path)
    pert_network_dict = validate.dictionary(pert_network_dict)

    with open(file, "w") as scores_file:
        writer = csv.writer(scores_file)

        for transf in sorted(pert_network_dict.keys()):
            score = pert_network_dict[transf]
            writer.writerow([transf[0], transf[1], score])

def write_network(pert_network_dict, protocol, file_path):
    """write the network to a file using the protocol for the pipeline.

    Args:
        pert_network_dict (dict): dictionary of perturbations with their LOMAP scores.
        protocol (pipeline.prep.pipeline_protocol): pipeline protocol. Needed for lambda windows and engine.
        file_path (str): file path to where the file should be written.
    """

    # validate inputs
    file = validate.string(file_path)
    pert_network_dict = validate.dictionary(pert_network_dict)
    protocol = validate.pipeline_protocol(protocol)

    # write perts file. Base the lambda schedule on the file generated in the previous cell.
    np.set_printoptions(formatter={'float': '{: .4f}'.format})

    with open(file, "w") as network_file:

        writer = csv.writer(network_file, delimiter=" ")
        
        for pert, lomap_score in pert_network_dict.items():
            # # based on the provided (at top of notebook) lambda allocations and LOMAP threshold, decide allocation.
            # if lomap_score == None or lomap_score < float(node.getInput("LOMAP Threshold")):
            #     num_lambda = node.getInput("DiffLambdaWindows")
            # else:
            #     num_lambda = node.getInput("LambdaWindows")
            
            num_lambda = protocol.num_lambda()# same no lamdda windows for all
        
            # given the number of allocated lambda windows, generate an array for parsing downstream.
            lam_array_np = np.around(np.linspace(0, 1, int(num_lambda)), decimals=5)

            # make the array into a format readable by bash.
            lam_array = str(lam_array_np).replace("[ ", "").replace("]", "").replace("  ", ",").replace('\n', '')

            # write out both directions for this perturbation.
            for eng in protocol.engines():
                writer.writerow([pert[0], pert[1], len(lam_array_np), lam_array, eng])
                # writer.writerow([pert[1], pert[0], len(lam_array_np), lam_array, engine])   