import sys
import os
import BioSimSpace as BSS

from ._validate import *

def engine_network(engine, file_path):
    """Generate a network file for only the engine specified.

    Args:
        engine (str): engine to generate the network for
        file_path (str): network.dat file containing multiple engines
    """
    engine = validate.engine(engine)
    file_path = validate.file_path(file_path)

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



class check_protocol():

    def __init__(self, file):
        # instantiate the class with the protocol file
        self._prot_file = validate.file_path(file)
        print(self._prot_file)
        self._query_dict = check_protocol._read_protocol(self)
        # update query dict if anything is missing
        self._query_dict = check_protocol._check_query(self)

    def _read_protocol(self):
        """reads the protocol file into a dictionary

        Raises:
            ValueError: if the protocol file does not contain "="

        Returns:
            dict: dictionary of the entries in the protocol file,
            seperated into key:value pairs by the "=" and units by "*".
        """

        query_dict = {}

        # get value regardless of the order of the protocol.dat
        with open(f"{self._prot_file}", "r") as file:
            for line in file:
                
                if "=" not in line:
                    raise ValueError("protocol file may only contain lines with '=', eg 'ligand forcefield = GAFF2'")
                elif "*" in line: # if unit, get as seperate dict entry
                    query_dict[f"{line.split('=')[0].strip().lower()}"] = line.split('=')[-1].strip().split("*")[0]
                    query_dict[f"{line.split('=')[0].strip().lower()} unit"] = line.split('=')[-1].strip().split("*")[-1]
                else:
                    query_dict[f"{line.split('=')[0].strip().lower()}"] = line.split('=')[-1].strip()

        return query_dict

    def _check_query(self):
        """fills in any gaps of the dict from the protocol file
        """
        
        query_dict = self._query_dict

        default_dict = {'ligand forcefield': 'gaff2',
                    'solvent': 'TIP3P',
                    'box edges': '30',
                    'box edges unit': 'angstrom',
                    'box type': 'cubic',
                    'sampling': '2',
                    'sampling unit': 'ns',
                    'hmr': 'False',
                    'repeats': '1',
                    'trajectories': "None",
                    'protein forcefield': 'ff14SB'
                    }
        
        for query in default_dict.keys():
            if query not in query_dict.keys():
                query_dict[query] = default_dict[query]
                print(f"{query} not found in protocol. {default_dict[query]} will be used.")
            
        return query_dict

    def validate(self):
        """validates all the input from the dict from the protocol file.
        """
        
        query_dict = self._query_dict

        try:
            self.ligand_forcefield = validate.lig_ff(query_dict['ligand forcefield'])
            self.protein_forcefield = validate.prot_ff(query_dict['protein forcefield'])
            self.solvent = validate.solvent_ff(query_dict['solvent'])
            self.box_edges = validate.integer(query_dict['box edges'])
            self.box_edges_unit = validate.box_edges_unit(query_dict['box edges unit'])
            self.box_type = validate.box_type(query_dict['box type'])
            self.sampling = validate.integer(query_dict['sampling'])
            self.sampling_unit = validate.sampling_unit(query_dict['sampling unit'])
            self.hmr = validate.boolean(query_dict['hmr'])
            self.repeats = validate.integer(query_dict['repeats'])
            self.trajectories = validate.trajectories(query_dict['trajectories'])

        except Exception as e:
            print(f"There is a problem with the input provided in {self._prot_file}.\n Exception is:\n {e}")
