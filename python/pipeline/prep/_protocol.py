import sys
import os
import csv
import BioSimSpace as BSS

from ..utils._validate import *


class pipeline_protocol():

    def __init__(self, file, auto_validate=False):
        """class for storing and validating protocol options from files or dictionary.

        Args:
            file (file path (or dictionary)): file with the options for the protocol. Alternatively a dictionary.
            auto_validate (bool, optional): Whether to automatically validate the given input. Defaults to False.
        """

        try:
            # instantiate the class with the protocol file
            self._prot_file = validate.file_path(file)
            # turn into a query dict as needed
            self._query_dict = self._read_protocol()
            is_file = True        
        except Exception as e:
            print(e)
            print("file was not recognised, trying to read as dictionary...")
            is_file = False
        
        if not is_file:
            try:
                self._query_dict = validate.dictionary(file)
                self._prot_file = None
            except Exception as e:
                print(e)
                raise TypeError("dictionary wasn't recognised either.")
        
        # update query dict if anything is missing
        self._query_dict = self._check_query(self.default_dict())

        # validated if needed
        auto_validate = validate.boolean(auto_validate)
        if auto_validate:
            self.validate()
            self._is_validated = True
        else:
            self._is_validated = False

    def default_dict(self):

        default_dict = {'ligand forcefield': 'gaff2',
                    'solvent': 'TIP3P',
                    'box edges': '30',
                    'box edges unit': 'angstrom',
                    'box type': 'cubic',
                    'sampling': '2',
                    'sampling unit': 'ns',
                    'hmr': 'False',
                    'hmr factor': 'auto',
                    'timestep overwrite': 'True',
                    'timestep' : "2",
                    "timestep unit" : "fs",
                    'repeats': '1',
                    'trajectories': "None",
                    'protein forcefield': 'ff14SB',
                    'start temperature': '0',
                    'end temperature': '300',
                    'temperature': '300',
                    'temperature unit': 'kelvin',
                    'pressure': '1',
                    'pressure unit': 'bar',
                    'minimisation steps': '10000',
                    'equilibrium runtime': '100',
                    'equilibrium runtime unit': 'ps'
                    }
        
        return default_dict

    def dictionary(self):
        """calls the dictionary made from the file, validated or not.

        Returns:
            dict: key value pairs for the given protocol
        """
        return self._query_dict


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


    def _check_query(self, default_dict):
        """fills in any gaps of the dict from the protocol file
        """
        
        query_dict = self._query_dict

        # remove any unrecognised dictionary entries
        for query in list(query_dict):
            if query not in default_dict.keys():
                del query_dict[query]
                print(f"{query} removed from the protocol as not recognised.\n please use only:\n {default_dict.keys()}")

        # add any missing protocol entries
        for query in default_dict.keys():
            if query not in query_dict.keys():
                query_dict[query] = default_dict[query]
                print(f"{query} not found in protocol. {default_dict[query]} will be used.")
        
            
        return query_dict

    
    def rewrite_protocol(self, file_path=None):
        """ rewrite the protocol to the same file after validation
        """
        if file_path:
            new_file = validate.string(file_path)
        else:
            if self._prot_file:
                new_file = self._prot_file
            else:
                raise ValueError("there is no file to overwrite, please provide a 'file_path'")

        query_dict = self._query_dict

        with open(new_file, "w") as protocol_file:
            writer = csv.writer(protocol_file)
            for query in query_dict.keys():
                writer.writerow([f"{query} = {query_dict[query]}"])


    def validate(self):
        """validates all the input from the dict from the protocol file.
        """
        
        query_dict = self._query_dict

        try:
            # validate all the input dict and replace in the query dict
            self.ligand_forcefield = validate.lig_ff(query_dict['ligand forcefield'])
            self.protein_forcefield = validate.prot_ff(query_dict['protein forcefield'])
            self.solvent = validate.solvent_ff(query_dict['solvent'])
            self.box_edges = validate.integer(query_dict['box edges'])
            self.box_edges_unit = validate.box_edges_unit(query_dict['box edges unit'])
            self.box_type = validate.box_type(query_dict['box type'])
            self.sampling = validate.integer(query_dict['sampling'])
            self.sampling_unit = validate.time_unit(query_dict['sampling unit'])
            self.repeats = validate.integer(query_dict['repeats'])
            self.trajectories = validate.trajectories(query_dict['trajectories'])
            self.start_temperature = validate.integer(query_dict['start temperature'])
            self.end_temperature = validate.integer(query_dict['end temperature'])
            self.temperature = validate.integer(query_dict['temperature'])
            self.temperature_unit = validate.temperature_unit(query_dict['temperature unit'])
            self.pressure = validate.integer(query_dict['pressure'])
            self.pressure_unit = validate.pressure_unit(query_dict['pressure unit'])
            self.min_steps = validate.integer(query_dict['minimisation steps'])
            self.eq_runtime = validate.integer(query_dict['equilibrium runtime'])
            self.eq_runtime_unit = validate.time_unit(query_dict['equilibrium runtime unit'])
            
            # choose timestep based on whether HMR is applied or not
            # this is important as BSS hmr mixin considers the timestep for the default auto
            self.hmr = validate.boolean(query_dict['hmr'])
            if self.hmr:
                if query_dict['hmr factor'].lower() == "auto":
                    self.hmr_factor = "auto"
                else:
                    self.hmr_factor = validate.is_float(query_dict['hmr factor'])
            self.timestep_overwrite = validate.boolean(query_dict['timestep overwrite'])
            self.timestep_unit = validate.time_unit(query_dict['timestep unit'])
            self.timestep = validate.integer(query_dict["timestep"])

            # will overwrite the provided timestep based on default values
            if self.timestep_overwrite:
                self.timestep_unit = validate.time_unit("fs")
                if self.hmr:
                    self.timestep = validate.integer(4)
                else:
                    self.timestep = validate.integer(2)

            # is now validated
            self._is_validated = True

        except ValueError as e:
            print(f"There is a problem with the input provided in {self._prot_file}.\n Error is:\n {e}")

    
    def print_protocol(self):

        if not self._is_validated:
            print("please validate the protocol first.")
        
        else:
            query_dict = self._query_dict

            for query in query_dict.keys():
                print(f"{query} : {query_dict[query]}")


class analysis_protocol(pipeline_protocol):

    def __init__(self, file, auto_validate=False):
        # inherit the init from other protocol too
        super().__init__(file, auto_validate)

    # read protocol inherited
    # rewrite protocol also inherited
    # print protocol also inherited
    # check query also inherited but with new default dict

    def default_dict(self):

        default_dict = {'estimator': "MBAR",
                        "method":"alchemlyb",
                        "check overlap":True,
                        "try pickle":True,
                        'save pickle':True,
                        "auto equilibration": False,
                        "statistical inefficiency": False,
                        "truncate percentage": 0,
                        "truncate keep":"end",
                        "mbar method": None # robust or default
                        }
        
        return default_dict

    def validate(self):
        """validates all the input from the dict from the protocol file.
        """
        
        query_dict = self._query_dict

        try:
            # validate all the input dict
            self.estimator = validate.estimator(query_dict['estimator'])
            self.method = validate.analysis_method(query_dict['method'])
            self.check_overlap = validate.boolean(query_dict['check overlap'])
            self.try_pickle = validate.boolean(query_dict['try pickle'])
            self.save_pickle = validate.boolean(query_dict['save pickle'])
            self.auto_equililbration = validate.boolean(query_dict['auto equilibration']) 
            self.statistical_inefficiency = validate.boolean(query_dict['statistical inefficiency'])
            self.truncate_percentage = validate.integer(query_dict["truncate percentage"])
            self.truncate_keep = validate.truncate_keep(query_dict["truncate keep"])
            self.mbar_method = validate.mbar_method(query_dict["mbar method"])
            
            # is now validated
            self._is_validated = True

        except ValueError as e:
            print(f"There is a problem with the input provided in {self._prot_file}.\n Error is:\n {e}")

