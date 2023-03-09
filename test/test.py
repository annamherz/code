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

print(pipeline.__file__)

from pipeline import *
from pipeline.analysis import *
from pipeline.prep import *
from pipeline.utils import *

transf = ["lig_ejm31~lig_ejm45", "lig_ejm44~lig_ejm45","lig_ejm44~lig_ejm49"]
engine = "GROMACS"
main_dir = "/home/anna/Documents/benchmark/tyk2_benchmark"
# methods = ["1fs", "2fs_HMR4", "4fs_HMR4", "4fs_HMR3", "2fs_HMR3", "2fs"]
methods = ["2fs"]

# options
ana_file = f"{main_dir}/execution_model/analysis_protocol.dat"
analysis_options = analysis_protocol(ana_file, auto_validate=True)
# analysis_options.rewrite_protocol()

for method in methods:
    for trans in transf:
        path_to_dir = f"{main_dir}/outputs_extracted/{engine}/{method}/{trans}"
        print(path_to_dir)
        final_results_folder = f"{main_dir}/outputs_extracted/results"

        # try:
        # using the pipeline module for analysis
        analysed_pert = analyse(path_to_dir)
        analysed_pert.set_options(analysis_options)
        # analysed_pert.set_options({"try pickle":False})
        avg, error, repeats_tuple_list = analysed_pert.analyse_all_repeats()
        # analysed_pert.plot_graphs()
        data_point_avg = [analysed_pert.ligand_0, analysed_pert.ligand_1,
                    analysed_pert.freenrg, analysed_pert.error,
                    analysed_pert.engine, analysed_pert.file_ext, method]
        print(data_point_avg)
        # write_analysis_file(analysed_pert, final_results_folder, method=method)
        # except Exception as e:
        #     print(e)
        #     print(f"could not analyse {path_to_dir}")


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
