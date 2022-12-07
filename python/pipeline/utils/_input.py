import sys
import os
import BioSimSpace as BSS

# TODO import validate
from ._validate import *

def engine_network(engine, file_path):
    """Generate a network file for only the engine specified.

    Args:
        engine (str): engine to generate the network for
        file_path (str): network.dat file containing multiple engines
    """
    validate_query.engine(engine)
    validate_query.path(file_path)

    try:
        print("using the folder of the file as the location to write the output file...")
        output_file = f"{file_path.rsplit('.',1)[0]}_{engine.lower()}.dat"

    except:
        print("assuming the file path is just network_combined.dat as not provided...")
        file_path = "network_combined.dat"
        if not os.path.exists(file_path):
            raise ValueError(f"{file_path} does not exist in the current folder.")
        output_file = f"network_combined_{engine.lower()}.dat"

    with open(file_path, "r") as file:
        with open(output_file, "w") as f:
                for line in file:
                        if engine in line:
                                f.write(f"{line}")


def read_protocol(prot_file):

    # Read the protocol.dat to get all parameter infos.
    search_dict = {"ligand forcefield":None, "protein forcefield":None,
                    "solvent":None, "box edges":None, "box type": None,
                    "protocol":None, "sampling":None, "HMR":None,
                    "repeats":None, "keep trajectories":None}
    
    # get value regardless of the order of the protocol.dat
    for search in search_dict.keys():
        with open(f"{prot_file}", "r") as file:
            for line in file:
                if search.casefold() in line.casefold():
                    search_dict[search] = line.strip()
                    break
    
    # get all the queries
    ligff_query = (search_dict["ligand forcefield"].rstrip().replace(" ", "").split("=")[-1]).lower() # ligand ff
    protff_query = search_dict["protein forcefield"].rstrip().replace(" ", "").split("=")[-1] # protein ff
    solvent_query = search_dict["solvent"].rstrip().replace(" ", "").split("=")[-1] # solvent ff
    boxsize_query = search_dict["box edges"].rstrip().replace(" ", "").split("=")[-1] # box size and type
    box_axis_length = boxsize_query.split("*")[0]
    box_axis_unit_query = boxsize_query.split("*")[1]
    boxtype_query = (search_dict["box type"].rstrip().replace(" ", "").split("=")[-1]).lower()
    sampling_query = search_dict["sampling"].rstrip().replace(" ", "").split("=")[-1].split("*")[0]
    sampling_unit_query = search_dict["sampling"].rstrip().replace(" ", "").split("=")[-1].split("*")[1]
    hmr_query = search_dict["HMR"].rstrip().replace(" ", "").split("=")[-1]
    repeats_query = search_dict["repeats"].rstrip().replace(" ", "").split("=")[-1]

    # validate if the queries are permitted
    query_dict = { "ligand forcefield":ligff_query, "protein forcefield":protff_query,
                    "solvent":solvent_query, "box length":box_axis_length, "box unit":box_axis_unit_query,
                    "boxtype":boxtype_query, "sampling":sampling_query, 
                    "sampling unit":sampling_unit_query, "HMR":hmr_query, "repeats":repeats_query}

    query = validate_query(query_dict)
    try:
        query.validate_all() # validate all inputs provided in the protocol file
    except Exception as e:
        print(f"There is a problem with the input provided in {prot_file}.\n Exception is:\n {e}")

    # convert the queries into the correct format for using for the rest of the pipeline

    parameters_dict = { "ligand forcefield":ligff, "protein forcefield":protff,
                    "solvent":solvent, "box length":box_length, "box unit":box_unit,
                    "boxtype":boxtype, "sampling":sampling, 
                    "sampling unit":sampling_unit, "HMR":hmr, "repeats":repeats}

    return query_dict