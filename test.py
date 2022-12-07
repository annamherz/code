import sys

try:
    import pipeline
except:
    print("adding code to the pythonpath...")
    code = '/home/anna/Documents/code/python'
    if code not in sys.path:
        sys.path.insert(1, code)
    import pipeline

from pipeline.utils._validate import *

prot_dict = {"sampling":2}
try:
    validate_query.validate_protocol_dict(prot_dict)
except Exception as e:
    print(f"There is a problem with the input provided in file\n Exception is:\n {e}")