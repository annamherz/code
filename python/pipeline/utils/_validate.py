# validate various inputs

import BioSimSpace as BSS
from BioSimSpace.Units.Length import angstrom as _angstrom
import os

class validate_query():

    def __init__(self, sampling=None):
        self._sampling = sampling
  

    def path(self, file_path):
        if not isinstance(file_path, str):
            raise TypeError("'file_path' must be of type 'str'.")

        if not os.path.exists(file_path):
            raise ValueError(f"{file_path} does not exist!")


    def engine(engine):
        if not isinstance(engine, str):
            raise TypeError("'engine' must be of type 'str'.")

        if engine.upper() not in BSS.FreeEnergy.engines():
            raise ValueError(f"'engine' must be one of {BSS.FreeEnergy.engines()}.")
    

    def lig_ff(lig_ff):

        if not isinstance(lig_ff, str):
            raise TypeError("'lig_ff' must be of type 'str'.")

        lig_ff_list = ["sage","parsely","gaff","gaff2"]
        if lig_ff not in lig_ff_list:
            raise ValueError(f"'lig_ff' must be one of {lig_ff_list}.")


    def prot_ff(prot_ff):

        if not isinstance(prot_ff, str):
            raise TypeError("'prot_ff' must be of type 'str'.")

        prot_ff_list = ['ff03','ff14SB','ff99','ff99SB','ff99SBildn']
        if prot_ff not in prot_ff_list:
            raise ValueError(f"'prot_ff' must be one of {prot_ff_list}.")


    def solvent_ff(solvent_ff):

        if not isinstance(solvent_ff, str):
            raise TypeError("'solvent_ff' must be of type 'str'.")

        if solvent_ff not in BSS.Solvent.waterModels():
            raise ValueError(f"'prot_ff' must be one of {BSS.Solvent.waterModels()}.")


    def sampling(self, sampling=None):
        if not sampling:
            sampling = self._sampling
        if not isinstance(sampling, int):
            if not isinstance(sampling, str):
                raise ValueError("sampling must be of type 'int' or 'str'")
    

    def sampling_unit(self, sampling_unit):

        if not isinstance(sampling_unit, str):
            raise TypeError("'sampling_unit' must be of type 'str'.")

        sampling_unit_list = ["ns", "ps"]
        if sampling_unit not in sampling_unit_list:
            raise ValueError(f"'prot_ff' must be one of {sampling_unit_list}.")


    def validate_all(self):
        validate_query.sampling(self)




    # query_dict = { "ligand forcefield":ligff_query, "protein forcefield":protff_query,
    #                 "solvent":solvent_query, "box length":box_axis_length, "box unit":box_axis_unit_query,
    #                 "boxtype":boxtype_query, "sampling":sampling_query, 
    #                 "sampling unit":sampling_unit, "HMR":hmr, "repeats":repeats}

class validate(validate_query):

    def __init__(self, file_path=None):
        self._file_path = file_path

    def sampling(self, sampling):
        super().sampling(sampling)
        if isinstance(sampling,str):
            sampling_int = int(sampling)
        
        return sampling_int
    
    def path(self, file_path):
        super().path(file_path)
        return file_path


    # if sampling_unit_query == "ns":
    #     sampling_unit = BSS.Units.Time.nanosecond
    # elif sampling_unit_query == "ps":
    #     sampling_unit = BSS.Units.Time.picosecond

    # Check if HMR
    
    # if hmr == "True":
    #     hmr = True
    # elif hmr == "False":
    #     hmr = False
    # if not isinstance(hmr, bool):
    #     raise NameError("Hydrogen Mass Repartitioning must be of type bool."
    #                     + "on the eighth line of protocol.dat in the shape of (e.g.):\nHMR = True")