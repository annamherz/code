import sys
import os
import csv
import BioSimSpace as BSS

from ..utils._validate import *
from ..utils._files import *

class pipeline_protocol():

    def __init__(self, file=None, auto_validate=False, verbose=True):
        """class for storing and validating protocol options from files or dictionary.

        Args:
            file (file path (or dictionary)): file with the options for the protocol. Alternatively a dictionary.
            auto_validate (bool, optional): Whether to automatically validate the given input. Defaults to False.
        """

        self.verbose = validate.boolean(verbose)

        if file:
            try:
                # instantiate the class with the protocol file
                self._prot_file = validate.file_path(file)
                # turn into a query dict as needed
                self._query_dict = self._read_protocol()
                is_file = True        
            except Exception as e:
                print(e)
                print("Not recognised as file, trying to read as dictionary...")
                is_file = False
            
            if not is_file:
                try:
                    self._query_dict = validate.dictionary(file)
                    self._prot_file = None
                except Exception as e:
                    print(e)
                    raise TypeError("dictionary wasn't recognised either.")
        else:
            print("no file or dict passed, using entirely default values...")
            # set the query dict as the default dict
            self._query_dict = self.default_dict()
        
        # update query dict if anything is missing
        self._query_dict = self._check_query()

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
                    'equilibrium runtime unit': 'ps',
                    'engine':"ALL",
                    "fepprep":"start"
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


    def _check_query(self):
        """fills in any gaps of the dict from the protocol file
        """
        
        query_dict = self._query_dict
        default_dict = self.default_dict()

        # remove any unrecognised dictionary entries
        for query in list(query_dict):
            if query not in default_dict.keys():
                del query_dict[query]
                if self.verbose:
                    print(f"{query} removed from the protocol as not recognised.\n please use only:\n {default_dict.keys()}")

        # add any missing protocol entries
        for query in default_dict.keys():
            if query not in query_dict.keys():
                query_dict[query] = default_dict[query]
                if self.verbose:
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

        write_protocol(self._query_dict, new_file)


    def validate(self):
        """validates all the input from the current query dict (taken from the protocol file or dict).
        """
        
        query_dict = self._query_dict

        try:
            # validate all the input dict and replace in the query dict
            self.ligand_forcefield(query_dict['ligand forcefield'])
            self.protein_forcefield(query_dict['protein forcefield'])
            self.solvent(query_dict['solvent'])
            self.box_edges(query_dict['box edges'])
            self.box_edges_unit(query_dict['box edges unit'])
            self.box_type(query_dict['box type'])
            self.sampling(query_dict['sampling'])
            self.sampling_unit(query_dict['sampling unit'])
            self.repeats(query_dict['repeats'])
            self.trajectories(query_dict['trajectories'])
            self.start_temperature(query_dict['start temperature'])
            self.end_temperature(query_dict['end temperature'])
            self.temperature(query_dict['temperature'])
            self.temperature_unit(query_dict['temperature unit'])
            self.pressure(query_dict['pressure'])
            self.pressure_unit(query_dict['pressure unit'])
            self.min_steps(query_dict['minimisation steps'])
            self.eq_runtime(query_dict['equilibrium runtime'])
            self.eq_runtime_unit(query_dict['equilibrium runtime unit'])
            self.engine(query_dict['engine'])
            self.fepprep(query_dict["fepprep"])
            
            # choose timestep based on whether HMR is applied or not
            # this is important as BSS hmr mixin considers the timestep for the default auto
            self.hmr(query_dict['hmr'])
            if self._hmr:
                if query_dict['hmr factor'].lower() == "auto":
                    self._hmr_factor = "auto"
                else:
                    self.hmr_factor(query_dict['hmr factor'])

            self.timestep(query_dict["timestep"])
            self.timestep_unit(query_dict['timestep unit'])
            self.timestep_overwrite(query_dict['timestep overwrite'])

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


    # changing any protocol settings
    # validate and update the internal query dictionary

    # this is not part of the default protocol and is found in the network file
    # it needs to be allocated before fepprep
    def fepprep(self, value=None):

        if value:
            options_list = ["start","middle","both"]
            if value not in options_list:
                raise ValueError(f"fepprep option must be in {options_list}")
            self._query_dict["fepprep"] = value
            self._fepprep = value
        else:
            value = self._fepprep

        return value

    def num_lambda(self, value=None):

        if value:
            value = validate.integer(value)
            self._query_dict["num lambda"] = value
            self._num_lambda = value
        else:
            value = self._num_lambda

        return value

    def engine(self, value=None):

        if value:
            value = validate.engines(value)
            self._query_dict["engine"] = value
            self._engine = value
        else:
            value = self._engine

        return value

    def ligand_forcefield(self, value=None):

        if value:
            value = validate.lig_ff(value)
            self._query_dict["ligand forcefield"] = value
            self._ligand_forcefield = value
        else:
            value = self._ligand_forcefield

        return value

    def protein_forcefield(self, value=None):

        if value:
            value = validate.prot_ff(value)
            self._query_dict["protein forcefield"] = value
            self._protein_forcefield = value
        else:
            value = self._protein_forcefield

        return value

    def solvent(self, value=None):

        if value:
            value = validate.solvent_ff(value)
            self._query_dict["solvent"] = value
            self._solvent = value
        else:
            value = self._solvent

        return value

    def box_edges(self, value=None):

        if value:
            value = validate.integer(value)
            self._query_dict["box edges"] = value
            self._box_edges = value
        else:
            self._box_edges

        return value

    def box_edges_unit(self, value=None):

        if value:
            value = validate.box_edges_unit(value)
            self._query_dict["box edges unit"] = value
            self._box_edges_unit = value
        else:
            value = self._box_edges_unit

        return value

    def box_type(self, value=None):

        if value:
            value = validate.box_type(value)
            self._query_dict["box type"] = value
            self._box_type = value
        else:
            value = self._box_type

        return value

    def sampling(self, value=None):

        if value:
            value = validate.integer(value)
            self._query_dict["sampling"] = value
            self._sampling = value
        else:
            value = self._sampling

        return value

    def sampling_unit(self, value=None):

        if value:
            value = validate.time_unit(value)
            self._query_dict["sampling unit"] = value
            self._sampling_unit = value
        else:
            value = self._sampling_unit

        return value    
    
    def repeats(self, value=None):

        if value:
            value = validate.integer(value)
            self._query_dict["repeats"] = value
            self._repeats = value
        else:
            value = self._repeats

        return value  
    
    def trajectories(self, value=None):

        if value:
            value = validate.trajectories(value)
            self._query_dict["trajectories"] = value
            self._trajectories = value
        else:
            value = self._trajectories

        return value  

    def start_temperature(self, value=None):

        if value:
            value = validate.integer(value)
            self._query_dict["start temperature"] = value
            self._start_temperature = value
        else:
            value = self._start_temperature

        return value  

    def end_temperature(self, value=None):

        if value:
            value = validate.integer(value)
            self._query_dict["end temperature"] = value
            self._end_temperature = value
        else:
            value = self._end_temperature

        return value  

    def temperature(self, value=None):

        if value:
            value = validate.integer(value)
            self._query_dict["temperature"] = value
            self._temperature = value
        else:
            value = self._temperature

        return value  

    def temperature_unit(self, value=None):

        if value:
            value = validate.temperature_unit(value)
            self._query_dict["temperature unit"] = value
            self._temperature_unit = value
        else:
            value = self._temperature_unit

        return value  

    def pressure(self, value=None):

        if value:
            value = validate.integer(value)
            self._query_dict["pressure"] = value
            self._pressure = value
        else:
            value = self._pressure

        return value  

    def pressure_unit(self, value=None):

        if value:
            value = validate.pressure_unit(value)
            self._query_dict["pressure unit"] = value
            self._pressure_unit = value
        else:
            value = self._pressure_unit

        return value  

    def min_steps(self, value=None):

        if value:
            value = validate.integer(value)
            self._query_dict["minimisation steps"] = value
            self._min_steps = value
        else:
            value = self._min_steps

        return value  

    def eq_runtime(self, value=None):

        if value:
            value = validate.integer(value)
            self._query_dict["equilibrium runtime"] = value
            self._eq_runtime = value
        else:
            value = self._eq_runtime

        return value  

    def eq_runtime_unit(self, value=None):

        if value:
            value = validate.time_unit(value)
            self._query_dict["equilibrium runtime unit"] = value
            self._eq_runtime_unit = value
        else:
            value = self._eq_runtime_unit

        return value  

    def hmr(self, value=None):

        if value:
            value = validate.boolean(value)
            self._query_dict["hmr"] = value
            self._hmr = value
        else:
            value = self._hmr

        return value  

    def hmr_factor(self, value=None):

        if value:
            if self._hmr:
                if value.lower() == "auto":
                    self._hmr_factor = "auto"
                else:
                    try:
                        self._hmr_factor = validate.is_float(value)
                    except:
                        raise ValueError("hmr_factor must be 'auto' or a integer/float")
                
                self._query_dict['hmr factor'] = value

                return value  
            
            else:
                print("'hmr' must be set to True for a hmr factor to be applied")
                return
        
        else:
            value = self._hmr_factor
            return value


    def timestep_overwrite(self, value=None):

        if value:
            value = validate.boolean(value)
            self._query_dict["timestep overwrite"] = value
            self._timestep_overwrite = value
        else:
            value = self._timestep_overwrite

        # will overwrite the provided timestep based on default values
        if self._timestep_overwrite:
            self.timestep_unit("fs")
            if self._hmr:
                self.timestep(4)
            else:
                self.timestep(2)

        return value

    def timestep_unit(self, value=None):

        if value:
            value = validate.time_unit(value)
            self._query_dict["timestep unit"] = value
            self._timestep_unit = value
        else:
            value = self._timestep_unit

        return value

    def timestep(self, value=None):

        if value:
            value = validate.integer(value)
            self._query_dict["timestep"] = value
            self._timestep = value
        else:
            value = self._timestep

        return value        

class analysis_protocol(pipeline_protocol):

    def __init__(self, file=None, auto_validate=False, verbose=False):
        # inherit the init from other protocol too
        super().__init__(file, auto_validate, verbose)

    # read protocol inherited
    # rewrite protocol also inherited
    # print protocol also inherited
    # check query also inherited but with new default dict

    def default_dict(self):

        default_dict = {'estimator': "MBAR",
                        "method":"alchemlyb",
                        "check overlap":"True",
                        "try pickle":"True",
                        'save pickle':"True",
                        "auto equilibration": "False",
                        "statistical inefficiency": "False",
                        "truncate percentage": "0",
                        "truncate keep":"end",
                        "mbar method": "None" # robust or default
                        }
        
        return default_dict

    def validate(self):
        """validates all the input from the dict from the protocol file.
        """
        
        query_dict = self._query_dict

        try:
            # validate all the input dict
            self.estimator(query_dict['estimator'])
            self.analysis_method(query_dict['method'])
            self.check_overlap(query_dict['check overlap'])
            self.try_pickle(query_dict['try pickle'])
            self.save_pickle(query_dict['save pickle'])
            self.auto_equilibration(query_dict['auto equilibration']) 
            self.statistical_inefficiency(query_dict['statistical inefficiency'])
            self.truncate_percentage(query_dict["truncate percentage"])
            self.truncate_keep(query_dict["truncate keep"])
            self.mbar_method(query_dict["mbar method"])
            
            # is now validated
            self._is_validated = True

        except ValueError as e:
            print(f"There is a problem with the input provided in {self._prot_file}.\n Error is:\n {e}")


    def estimator(self, value=None):

        if value:
            value = validate.estimator(value)
            self._query_dict["estimator"] = value
            self._estimator = value
        else:
            value = self._estimator

        return value
    
    def analysis_method(self, value=None):

        if value:
            value = validate.analysis_method(value)
            self._query_dict["method"] = value
            self._method = value
        else:
            value = self._method

        return value

    def check_overlap(self, value=None):

        if value:
            value = validate.boolean(value)
            self._query_dict["check overlap"] = value
            self._check_overlap = value
        else:
            value = self._check_overlap

    def try_pickle(self, value=None):

        if value:
            value = validate.boolean(value)
            self._query_dict["try pickle"] = value
            self._try_pickle = value
        else:
            value = self._try_pickle

        return value

    def save_pickle(self, value=None):

        if value:
            value = validate.boolean(value)
            self._query_dict["save pickle"] = value
            self._save_pickle = value
        else:
            value = self._save_pickle

        return value
    
    def auto_equilibration(self, value=None):

        if value:
            value = validate.boolean(value)
            self._query_dict["auto equilibration"] = value
            self._auto_equilibration = value
        else:
            value = self._auto_equilibration

        return value

    def statistical_inefficiency(self, value=None):

        if value:
            value = validate.boolean(value)
            self._query_dict["statistical inefficiency"] = value
            self._statistical_inefficiency = value
        else:
            value = self._statistical_inefficiency

        return value

    def truncate_percentage(self, value=None):

        if value:
            value = validate.integer(value)
            self._query_dict["truncate percentage"] = value
            self._truncate_percentage = value
        else:
            value = self._truncate_percentage

        return value

    def truncate_keep(self, value=None):

        if value:
            value = validate.truncate_keep(value)
            self._query_dict["truncate keep"] = value
            self._truncate_keep = value
        else:
            value = self._truncate_keep

        return value

    def mbar_method(self, value=None):

        if value:
            value = validate.mbar_method(value)
            self._query_dict["mbar method"] = value
            self._mbar_method = value
        else:
            value = self._mbar_method

        return value
