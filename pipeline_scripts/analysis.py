# analysis script for a single pert
# anna

import BioSimSpace as BSS
import sys
import os as _os

BSS.setVerbose = True

try:
    import pipeline
except:
    print("adding code to the pythonpath...")
    code = '/home/anna/Documents/code/python'
    if code not in sys.path:
        sys.path.insert(1, code)
    import pipeline

from pipeline.analysis import *
from pipeline.utils import write_analysis_file
from pipeline.prep import *

# someway that this is for the transformations listed in the csv file made at the start for the setup
# or alternatively do this in bash
trans = sys.argv[1].rstrip()
engine = sys.argv[2].rstrip()

ana_file = os.environ["ana_file"]

# options
analysis_options = analysis_protocol(ana_file, auto_validate=True)
analysis_options.rewrite_protocol()

main_dir = _os.environ["MAINDIRECTORY"]

# TODO final results folder has a method for the analysis option used?

# find correct path, use extracted if it exists
if _os.path.exists(f"{main_dir}/outputs_extracted/{engine}/{trans}"):
    path_to_dir = f"{main_dir}/outputs_extracted/{engine}/{trans}"
    final_results_folder = f"{main_dir}/outputs_extracted/results"
elif _os.path.exists(f"{main_dir}/outputs/{engine}_extracted/{trans}"):
    path_to_dir = f"{main_dir}/outputs/{engine}_extracted/{trans}"
    final_results_folder = f"{main_dir}/outputs/results"
else:
    path_to_dir = f"{main_dir}/outputs/{engine}/{trans}"
    final_results_folder = f"{main_dir}/outputs/results"

if not _os.path.exists(path_to_dir):
    raise OSError(f"{path_to_dir} does not exist.")

print(f'analysing results for {path_to_dir}')
print(f"using {analysis_options} for analysis")

# using the pipeline module for analysis
analysed_pert = analyse(path_to_dir)
analysed_pert.set_options(analysis_options)
avg, error, repeats_tuple_list = analysed_pert.analyse_all_repeats()
analysed_pert.plot_graphs()

# write the final result
write_analysis_file(analysed_pert, final_results_folder)
# TODO change so also class method