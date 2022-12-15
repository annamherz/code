import sys
import os
import csv
import BioSimSpace as BSS

from ..utils._validate import *


class pipeline_protocol():

    def __init__(self, file):
        # instantiate the class with the protocol file
        self._prot_file = validate.file_path(file)
        print(self._prot_file)
        self._query_dict = pipeline_protocol._read_protocol(self)
        # update query dict if anything is missing
        self._query_dict = pipeline_protocol._check_query(self)


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

    
    def rewrite_protocol(self):
        """ rewrite the protocol to the same file after validation
        """
        new_file = self._prot_file
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
            # validate all the input dict
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
            self.timestep_unit = validate.time_unit("fs")
            if self.hmr:
                self.timestep = validate.integer(4)
            else:
                self.timestep = validate.integer(2)


        except ValueError as e:
            print(f"There is a problem with the input provided in {self._prot_file}.\n Error is:\n {e}")
