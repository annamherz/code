import BioSimSpace as BSS
import sys

try:
    import pipeline
except:
    print("adding code to the pythonpath...")
    code = '/home/anna/Documents/code/python'
    if code not in sys.path:
        sys.path.insert(1, code)
    import pipeline

from pipeline import *

file = "/home/anna/Documents/benchmark/tyk2_benchmark/execution_model_rbfenn_test/protocol.dat"

protocol = utils.check_protocol(file) # instantiate the protocol as an object
protocol.validate() # validate all the input

protocol.ligand_forcefield