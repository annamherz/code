import BioSimSpace as BSS
from BioSimSpace._SireWrappers import System as _System
import os
import warnings
import pipeline


class validate():
    """class of staticmethods to return validated input
    """

    def __init__(self):
        pass


    @staticmethod
    def string(string):
        """validates the string

        Args:
            string (str): text

        Raises:
            TypeError: must be of type 'str'

        Returns:
            str: text
        """
        if not isinstance(string, str):
            raise TypeError(f"{string} / 'string' must be of type 'str'.")
        
        return string


    @staticmethod
    def file_path(file_path, create=False):
        """validates the provided file path

        Args:
            file_path (str): path to file

        Raises:
            TypeError: must be of type 'str'
            ValueError: path must exist!

        Returns:
            str: file path
        """
        if not isinstance(file_path, str):
            raise TypeError("'file_path' must be of type 'str'.")

        if not os.path.exists(file_path):
            raise ValueError(f"{file_path} does not exist!")
        
        return file_path


    @staticmethod
    def folder_path(folder_path, create=False):
        """validates the provided file path

        Args:
            folder_path (str): path to file

        Raises:
            TypeError: must be of type 'str'
            ValueError: path must exist!

        Returns:
            str: file path
        """
        if not isinstance(folder_path, str):
            raise TypeError("'folder_path' must be of type 'str'.")

        if create:
            if not os.path.exists(folder_path):
                warnings.warn(f"{folder_path} does not exist!")
                os.makedirs(folder_path)
            else:
                pass   
        else:
            if not os.path.exists(folder_path):
                raise ValueError(f"{folder_path} does not exist!")
        
        return folder_path


    @staticmethod
    def engine(engine):
        """validates the provided engine

        Args:
            engine (str): engine to be used for the MD / FEP run

        Raises:
            TypeError: must be of type 'str'
            ValueError: must be one of BSS.FreeEnergy.engines().

        Returns:
            str: engine.upper()
        """
        if not isinstance(engine, str):
            raise TypeError("'engine' must be of type 'str'.")

        if engine.upper() not in BSS.FreeEnergy.engines():
            raise ValueError(f"'engine' must be one of {BSS.FreeEnergy.engines()}.")
        
        return engine.upper()
    

    @staticmethod
    def lig_ff(lig_ff):
        """validates the provided ligand forcefield

        Args:
            lig_ff (str): name of the ligand forcefield

        Raises:
            TypeError: must be of type 'str'
            ValueError: must be either sage, parsely, gaff or gaff2.

        Returns:
            str: forcefield that can be used with BSS.Parameters.parameterise()
        """

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
        """validates the provided protein forcefield

        Args:
            prot_ff (str): name of the protein forcefield

        Raises:
            TypeError: must be of type 'str'
            ValueError: must be one of the BSS.Parameters.forcefields() protein forcefields.

        Returns:
            str: forcefield that can be used with BSS.Parameters.parameterise()
        """

        if not isinstance(prot_ff, str):
            raise TypeError("'prot_ff' must be of type 'str'.")

        prot_ff_list = ['ff03','ff14SB','ff99','ff99SB','ff99SBildn']
        if prot_ff not in prot_ff_list:
            raise ValueError(f"'prot_ff' must be one of {prot_ff_list}.")

        return prot_ff


    @staticmethod
    def solvent_ff(solvent_ff):
        """validates the provided solvent forcefield / water model

        Args:
            solvent_ff (str): name of the water model

        Raises:
            TypeError: must be of type 'str'
            ValueError: must be in BSS.Solvent.waterModels()

        Returns:
            str: water model that can be used with BSS.Solvent.solvate()
        """

        if not isinstance(solvent_ff, str):
            raise TypeError("'solvent_ff' must be of type 'str'.")

        if solvent_ff.lower() not in BSS.Solvent.waterModels():
            raise ValueError(f"'solvent_ff' must be one of {BSS.Solvent.waterModels()}.")

        return solvent_ff.lower()

    @staticmethod
    def integer(integer):
        """validates the provided integer so that it is an integer

        Args:
            integer (str or int or float): an integer / number to be converted into an integer

        Raises:
            TypeError: must be of type 'str' or 'int' or 'float'
            ValueError: 'str' could not be converted into an integer

        Returns:
            int: integer in integer format
        """

        if not isinstance(integer, int):
            if not isinstance(integer, float):
                if not isinstance(integer, str):
                    raise TypeError(f"{integer} must be of type 'int', 'float' or 'str'")
                else:
                    try:
                        if "." in integer:
                            integer = float(integer)
                        else:
                            integer = int(integer)
                    except:
                        raise ValueError(f"{integer} could not be converted into an integer")
        
        if isinstance(integer, float):
            print(f"{integer} will be turned into {int(integer)}")
            integer = round(integer)

        return integer


    @staticmethod
    def time_unit(time_unit):
        """validates the sampling unit

        Args:
            time_unit (str or BSS.Types.Time): unit of sampling

        Raises:
            TypeError: must be of type 'str' or 'BSS.Types.Time'
            ValueError: if not of BSS.Types.Time, must be 'ns' or 'ps'.

        Returns:
            BSS.Types.Time: in BSS Time
        """

        if not isinstance(time_unit, str):
            if not isinstance(time_unit, BSS.Types.Time):
                raise TypeError("'time_unit' must be of type 'str' or 'BSS.Types.Time'.")

        else:
            time_unit_list = ["ns", "ps", "fs"]
            if time_unit not in time_unit_list:
                raise ValueError(f"'time_unit' must be one of {time_unit_list}.")
            
            if time_unit == "ns":
                time_unit = BSS.Units.Time.nanosecond
            elif time_unit == "ps":
                time_unit = BSS.Units.Time.picosecond      
            elif time_unit == "fs":
                time_unit = BSS.Units.Time.femtosecond      

        return time_unit
       

    @staticmethod
    def box_edges_unit(box_edges_unit):
        """validates the box edges unit

        Args:
            box_edges_unit (str or BSS.Types.Length): unit of length

        Raises:
            TypeError: must be of type 'str' or 'BSS.Types.Length'
            ValueError: if not of BSS.Types.Length, must be 'nm' or 'angstrom'.

        Returns:
            BSS.Types.Length: in BSS Length
        """

        if not isinstance(box_edges_unit, str):
            if not isinstance(box_edges_unit, BSS.Types.Length):
                raise TypeError("'box_edges_unit' must be of type 'str' or 'BSS.Types.Length'.")

        else:
            box_edges_unit = box_edges_unit.lower()
            box_edges_unit_list = ["angstrom", "a", "nanometer", "nm"]
            if box_edges_unit not in box_edges_unit_list:
                raise ValueError(f"'box_edges_unit' must be one of {box_edges_unit_list}.")

            if box_edges_unit == "angstrom" or box_edges_unit == "a":
                box_edges_unit = BSS.Units.Length.angstrom
            elif box_edges_unit == "nm" or box_edges_unit == "nanometer":
                box_edges_unit = BSS.Units.Length.nanometer   

        return box_edges_unit


    @staticmethod
    def temperature_unit(temperature_unit):
        """validates the temperature unit

        Args:
            temperature_unit (str or BSS.Types.Temperature): unit of temperature

        Raises:
            TypeError: must be of type 'str' or 'BSS.Types.Temperature'
            ValueError: if not of BSS.Types.Temperature
        Returns:
            BSS.Types.Temperature: in BSS Temperature
        """

        if not isinstance(temperature_unit, str):
            if not isinstance(temperature_unit, BSS.Types.Temperature):
                raise TypeError("'temperature_unit' must be of type 'str' or 'BSS.Types.Temperature'.")

        else:
            temperature_unit = temperature_unit.lower()
            temperature_unit_list = ["kelvin", "k", "celsius", "c"]
            if temperature_unit not in temperature_unit_list:
                raise ValueError(f"'temperature_unit' must be one of {temperature_unit_list}.")

            if temperature_unit == "k" or temperature_unit == "kelvin":
                temperature_unit = BSS.Units.Temperature.kelvin
            elif temperature_unit == "c" or temperature_unit == "celsius":
                temperature_unit = BSS.Units.Temperature.celsius   

        return temperature_unit


    @staticmethod
    def pressure_unit(pressure_unit):
        """validates the pressure unit

        Args:
            pressure_unit (str or BSS.Types.Pressure): unit of pressure

        Raises:
            TypeError: must be of type 'str' or 'BSS.Types.Pressure'
            ValueError: if not of BSS.Types.Pressure
        Returns:
            BSS.Types.Pressure: in BSS Pressure
        """

        if not isinstance(pressure_unit, str):
            if not isinstance(pressure_unit, BSS.Types.Pressure):
                raise TypeError("'pressure_unit' must be of type 'str' or 'BSS.Types.Pressure'.")

        else:
            pressure_unit = pressure_unit.lower()
            pressure_unit_list = ["atm", "atmosphere", "bar"]
            if pressure_unit not in pressure_unit_list:
                raise ValueError(f"'pressure_unit' must be one of {pressure_unit_list}.")

            if pressure_unit == "atm" or pressure_unit == "atmosphere":
                pressure_unit = BSS.Units.Pressure.atm
            elif pressure_unit == "bar":
                pressure_unit = BSS.Units.Pressure.bar  

        return pressure_unit
        

    @staticmethod
    def box_type(box_type):
        """validates the box type

        Args:
            box_type (str): type of box to be used for solvation.

        Raises:
            TypeError: must be of type 'str'
            ValueError: must be one of the box types accepted by minimum_solvation()

        Returns:
            str: box type
        """

        if not isinstance(box_type, str):
            raise TypeError("'box_type' must be of type 'str'.")

        box_type_list = ["cubic", "truncatedOctahedron", "octahedral"]
        if box_type not in box_type_list:
            raise ValueError(f"'box_type' must be one of {box_type_list}.")

        return box_type


    @staticmethod
    def boolean(boolean):
        """validates boolean input

        Args:
            boolean (str or bool): boolean

        Raises:
            TypeError: must be of type 'str' or 'bool'
            ValueError: if 'str' , must be 'True'/'1' or 'False'/'0'

        Returns:
            bool: boolean
        """

        if boolean != 1 and boolean != 0:
            if not isinstance(boolean, bool):
                if not isinstance(boolean, str):
                    raise TypeError(f"{boolean} must be of type 'str' (True or False) or 'bool'.")
                else:
                    boolean_list = ["True", "False","1","0"]
                    if boolean not in boolean_list:
                        raise ValueError(f"{boolean} must be one of {boolean_list}.")
        else:
            boolean = str(boolean)


        if boolean == "True" or boolean == "1":
            boolean = True
        elif boolean == "False" or boolean == "0":
            boolean = False
        
        return boolean


    @staticmethod
    def trajectories(trajectories):
        """validate input for trajectories being saved or not

        Args:
            trajectories (str): whether to save the trajectory or not

        Raises:
            TypeError: must be of type 'str'
            ValueError: must be in trajectories_list so can be recognised by future scripts

        Returns:
           str: how many of the trajectories will be saved
        """

        if not isinstance(trajectories, str):
            raise TypeError("'trajectories' must be of type 'str'.")

        trajectories_list = ["None", "0,0.5,1", "0,1", "All"]
        if trajectories not in trajectories_list:
            raise ValueError(f"'trajectories' must be one of {trajectories_list}.")

        return trajectories


    @staticmethod
    def num_lambda(num_lambda):
        """validate number of lambdas to be run

        Args:
            num_lambda (str or int): number of lambdas

        Raises:
            TypeError: must be of type 'str' or 'int'
            ValueError: must be accepted number of lambda windows.

        Returns:
            int: number of lambdas
        """

        validate.integer(num_lambda)

        num_lambda_list = [11]
        if num_lambda not in num_lambda_list:
            raise ValueError(f"'num_lambda' must be one of {num_lambda}.")

        return num_lambda


    @staticmethod
    def pipeline_protocol(protocol, fepprep=False):
        """check if the passed protocol is a correct pipeline_protocol

        Args:
            protocol (pipeline.prep.pipeline_protocol): a previously
            read and validated pipeline protocol

            fepprep (bool) : 

        Returns:
            pipeline.prep.pipeline_protocol: the protocol if it is okay
        """

        if not isinstance(protocol, pipeline.prep.pipeline_protocol):
            raise TypeError("'protocol' must be of type 'pipeline.prep.pipeline_protocol'.")

        if not isinstance(fepprep, bool):
            raise TypeError("'fepprep' must be of type 'bool'.")

        # validate incase it wasnt
        protocol.validate()

        # if fepprep, check that an engine and no of lam was specified
        if fepprep:
            if not hasattr(protocol, "num_lambda"):
                warnings.warn("the provided protocol does not have attribute 'num_lambda'.\n 11 lambda windows will be used...")
                protocol.num_lambda = 11
            if not hasattr(protocol, "engine"):
                raise TypeError("protocol must have an engine to be used for fep.\n please set an engine using protocol.engine = 'ENGINE' ")
            else:
                protocol.engine = validate.engine(protocol.engine)

        return protocol


    @staticmethod
    def system(system):
        """checks if it is a BSS system

        Args:
            system (BioSimSpace._SireWrappers._system.System): the system

        Raises:
            TypeError: if not BioSimSpace._SireWrappers._system.System

        Returns:
            BioSimSpace._SireWrappers._system.System: the passed system
        """


        if not isinstance(system, _System):
            raise TypeError("'system' must be a BSS system!.")

        return system


    @staticmethod
    def analysis(analysis, analysed=True):
        """checks if it is an analysed pipeline.analysis.analyse 

        Args:
            system (pipeline.analysis.analyse): an analysed pipeline analysis

        Raises:
            TypeError: if not pipeline.analysis.analyse

        Returns:
            pipeline.analysis.analyse: the passed system
        """


        if not isinstance(analysis, pipeline.analysis.analyse):
            raise TypeError("'analysis' must be a pipeline.analysis.analyse!.")
        
        if analysed == True:
            if not analysis.is_analysed:
                raise ValueError("'analysis' must have already been analysed!.")

        return analysis