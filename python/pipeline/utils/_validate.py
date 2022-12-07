# validate various inputs

import BioSimSpace as BSS
from BioSimSpace.Units.Length import angstrom as _angstrom
import os


class validate():

    def __init__(self):
        pass
  
    @staticmethod
    def file_path(file_path):
        if not isinstance(file_path, str):
            raise TypeError("'file_path' must be of type 'str'.")

        if not os.path.exists(file_path):
            raise ValueError(f"{file_path} does not exist!")
        
        return file_path


    @staticmethod
    def engine(engine):
        if not isinstance(engine, str):
            raise TypeError("'engine' must be of type 'str'.")

        if engine.upper() not in BSS.FreeEnergy.engines():
            raise ValueError(f"'engine' must be one of {BSS.FreeEnergy.engines()}.")
        
        return engine.upper()
    

    @staticmethod
    def lig_ff(lig_ff):

        if not isinstance(lig_ff, str):
            raise TypeError("'lig_ff' must be of type 'str'.")
        else:
            lig_ff = lig_ff.lower()

        lig_ff_list = ["sage","parsely","gaff","gaff2", "openff_unconstrained-2.0.0", "openff_unconstrained-1.3.0", "openff"]
        if lig_ff not in lig_ff_list:
            raise ValueError(f"'lig_ff' must be one of {lig_ff_list}.")
        
        if lig_ff == "sage":
            lig_ff = "openff_unconstrained-2.0.0"
        elif lig_ff == "parsely":
            lig_ff = "openff_unconstrained-1.3.0"
        elif lig_ff == "openff":
            lig_ff = "openff_unconstrained-2.0.0"
        
        return lig_ff


    @staticmethod
    def prot_ff(prot_ff):

        if not isinstance(prot_ff, str):
            raise TypeError("'prot_ff' must be of type 'str'.")

        prot_ff_list = ['ff03','ff14SB','ff99','ff99SB','ff99SBildn']
        if prot_ff not in prot_ff_list:
            raise ValueError(f"'prot_ff' must be one of {prot_ff_list}.")

        return prot_ff


    @staticmethod
    def solvent_ff(solvent_ff):

        if not isinstance(solvent_ff, str):
            raise TypeError("'solvent_ff' must be of type 'str'.")

        if solvent_ff.lower() not in BSS.Solvent.waterModels():
            raise ValueError(f"'solvent_ff' must be one of {BSS.Solvent.waterModels()}.")

        return solvent_ff.lower()

    @staticmethod
    def sampling(sampling):

        if not isinstance(sampling, int):
            if not isinstance(sampling, str):
                raise ValueError("sampling must be of type 'int' or 'str'")

        return int(sampling)


    @staticmethod
    def sampling_unit(sampling_unit):

        if not isinstance(sampling_unit, str):
            raise TypeError("'sampling_unit' must be of type 'str'.")

        sampling_unit_list = ["ns", "ps"]
        if sampling_unit not in sampling_unit_list:
            raise ValueError(f"'sampling_unit' must be one of {sampling_unit_list}.")
        
        if sampling_unit == "ns":
            sampling_unit = BSS.Units.Time.nanosecond
        elif sampling_unit == "ps":
            sampling_unit = BSS.Units.Time.picosecond        

        return sampling_unit

        
    @staticmethod
    def box_edges(box_edges):

        if not isinstance(box_edges, int):
            if not isinstance(box_edges, str):
                raise TypeError("'box_edges' must be of type 'str' or 'int.")

        return int(box_edges)
        

    @staticmethod
    def box_edges_unit(box_edges_unit):

        if not isinstance(box_edges_unit, str):
            raise TypeError("'box_edges_unit' must be of type 'str'.")

        box_edges_unit_list = ["angstrom", "nm"]
        if box_edges_unit not in box_edges_unit_list:
            raise ValueError(f"'box_edges_unit' must be one of {box_edges_unit_list}.")

        if box_edges_unit == "angstrom":
            box_edges_unit = BSS.Units.Length.angstrom
        elif box_edges_unit == "nm":
            box_edges_unit = BSS.Units.Length.nanometer   

        return box_edges_unit


    @staticmethod
    def box_type(box_type):

        if not isinstance(box_type, str):
            raise TypeError("'box_type' must be of type 'str'.")

        box_type_list = ["cubic", "truncatedOctahedron"]
        if box_type not in box_type_list:
            raise ValueError(f"'box_type' must be one of {box_type_list}.")

        return box_type


    @staticmethod
    def hmr(hmr):

        if not isinstance(hmr, bool):
            if not isinstance(hmr, str):
                raise TypeError("'hmr' must be of type 'str' (True or False) or 'bool'.")
            else:
                hmr_list = ["True", "False"]
                if hmr not in hmr_list:
                    raise ValueError(f"'hmr' must be one of {hmr_list}.")
                
                if hmr == "True":
                    hmr = True
                elif hmr == "False":
                    hmr = False
        
        return hmr


    @staticmethod
    def repeats(repeats):

        if not isinstance(repeats, str):
            if not isinstance(repeats, int):
                raise TypeError("'repeats' must be of type 'str' or 'int'.")

        return int(repeats)


    @staticmethod
    def trajectories(trajectories):

        if not isinstance(trajectories, str):
            raise TypeError("'trajectories' must be of type 'str'.")

        trajectories_list = ["None", "0,0.5,1", "0,1", "All"]
        if trajectories not in trajectories_list:
            raise ValueError(f"'trajectories' must be one of {trajectories_list}.")

        return trajectories

    @staticmethod
    def num_lambda(num_lambda):

        if not isinstance(num_lambda, str):
            if not isinstance(num_lambda, int):
                raise TypeError("'num_lambda' must be of type 'str' or 'int'.")

        num_lambda_list = [11]
        if int(num_lambda) not in num_lambda_list:
            raise ValueError(f"'num_lambda' must be one of {num_lambda}.")

        return int(num_lambda)