# analysis script for plotting and stats for all the diff engines
# anna

from inspect import BoundArguments
import warnings
import sys
import os

if '/home/anna/Documents/cinnabar' not in sys.path:
    sys.path.insert(1, '/home/anna/Documents/cinnabar')
import cinnabar

try:
    import pipeline
except:
    print("adding code to the pythonpath...")
    code = '/home/anna/Documents/code/python'
    if code not in sys.path:
        sys.path.insert(1, code)
    import pipeline

from pipeline.analysis import *


# location of repeat files
# results folder for output
# main_dir = os.environ["MAINDIRECTORY"]
# net_file = os.environ["net_file"]
main_dir = "/backup/anna/benchmark/tyk2"
net_file = "/home/anna/Documents/benchmark/tyk2_benchmark/execution_model/network_combined.dat"

# choose location for the experimental file
exp_file = "/home/anna/Documents/benchmark/inputs/experimental/tyk2.yml"

if os.path.exists(f"{main_dir}/outputs_extracted/results"):
    results_folder = f"{main_dir}/outputs_extracted/results"
elif os.path.exists(f"{main_dir}/outputs/results"):
    results_folder = f"{main_dir}/outputs/results"
else:
    raise ValueError(f"results directory not found in the {main_dir}. please make sure results were written using the analysis script previously in the pipeline")

output_folder = f"{main_dir}/analysis"
# TODO format options for final analysis
# options that should be done in dict form
# options
# analysis_options = {"engines": "all",
#                     "bar":True,
#                     "scatter":True,
#                     "cycle_closure": True,
#                     "convergence": True,
#                     "outliers": 3,
#                     "XXX":"XXX",
#                     # "try_pickle":True,'save_pickle':True
#                     }

# TODO set these somewhere else externally so will use in both scripts?
# This should be the analysis options chosen during the analysis
analysis_options = {'estimator': "MBAR", "method":"alchemlyb",
                    "check_overlap":True,
                    "try_pickle":True, 'save_pickle':True,
                    "auto_equilibration": False,
                    "statistical_inefficiency": False,
                    "truncate_percentage": 0,
                    "truncate_keep":"start"}

file_ext = str(f"{analysis_options['estimator']}_{analysis_options['method']}_"+
            f"eq{str(analysis_options['auto_equilibration']).lower()}_"+
            f"stats{str(analysis_options['statistical_inefficiency']).lower()}_"+
            f"truncate{str(analysis_options['truncate_percentage'])}{analysis_options['truncate_keep']}")

all_analysis_object = analysis_network(results_folder,
                                       exp_file=None,
                                       net_file=net_file,
                                       output_folder=output_folder,
                                       file_ext=file_ext
                                        )
# all_analysis_object.set_options(analysis_options)
all_analysis_object.get_experimental(exp_file) # add experimental files

print("computing results...")
all_analysis_object.compute()

# can add any other results files
# all_analysis_object.compute_other_results(file_name=None, name=None)

# plotting
print("plotting results...")
# plot bar graphs for all engines
all_analysis_object.plot_bar_pert()
all_analysis_object.plot_bar_lig()

# plot scatter plots
all_analysis_object.plot_scatter_lig(use_cinnabar=True)
all_analysis_object.plot_scatter_pert(use_cinnabar=True)

# all_analysis_object.cycle_closure()

# put out statistics