# analysis script for a single pert
# anna

import BioSimSpace as BSS
import sys
import os as _os

BSS.setVerbose = True

if '/home/anna/Documents/cinnabar' not in sys.path:
    sys.path.insert(1, '/home/anna/Documents/cinnabar')
import cinnabar

print("adding code to the pythonpath...")
code = '/home/anna/Documents/code/python'
if code not in sys.path:
    sys.path.insert(1, code)
import pipeline

print(pipeline.__file__)

from pipeline import *
from pipeline.analysis import *
from pipeline.prep import *
from pipeline.utils import *

bench_folder = f"/home/anna/Documents/benchmark"
protein = "tyk2"
main_dir = f"{bench_folder}/extracted/{protein}"

# choose location for the files
net_file = f"{main_dir}/execution_model/network_lomap.dat"
ana_file = f"{main_dir}/execution_model/analysis_protocol.dat"
exp_file = f"{bench_folder}/inputs/experimental/{protein}.yml"

if os.path.exists(f"{main_dir}/outputs_extracted/results"):
    results_folder = f"{main_dir}/outputs_extracted/results"
elif os.path.exists(f"{main_dir}/outputs/results"):
    results_folder = f"{main_dir}/outputs/results"
else:
    raise ValueError(f"results directory not found in the {main_dir}. please make sure results were written using the analysis script previously in the pipeline")

output_folder = validate.folder_path(f"{main_dir}/analysis", create=True)

all_analysis_object = analysis_network(results_folder,
                                       exp_file=exp_file,
                                       net_file=net_file,
                                       output_folder=output_folder,
                                       analysis_ext=ana_file
                                        )

all_analysis_object.compute(cycle_closure=True)

# # fwf exp data
# print("fwf")
# exp_dicts = res_obj._get_exp_fwf(fwf_path='/home/anna/Documents/september_2022_workshops/freenrgworkflows/networkanalysis/')
# for key in exp_dicts[0]:
#     print(f"{key} : {exp_dicts[0][key][0]}")
#     # dG_kj_mol = dG_kj_mol - np.mean(dG_kj_mol) # nomralises the free energy by subtracting the mean from them
#     # in fwf

# print("self normalised")
# # print(res_obj.exper_val_dict)
# values_exp = []
# for val in res_obj.exper_val_dict.values():
#     values_exp.append(float(val[0]))
# avg_exp = np.mean(values_exp)
# for key in res_obj.exper_val_dict:
#     print(f"{key} : {float(res_obj.exper_val_dict[key][0]) - avg_exp}")

# x_data = np.asarray([node[1]["exp_DG"] for node in res_obj._cinnabar_networks["SOMD"].graph.nodes(data=True)])


# for node in res_obj._cinnabar_networks["SOMD"].graph.nodes(data=True):
#     print(node[1]["name"], node[1]["exp_DG"], ( float(node[1]["exp_DG"]) - np.mean(x_data)) )

# x_data = x_data - np.mean(x_data)
# print(x_data)

# print(res_obj.normalised_exper_val_dict["SOMD"])
