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

ligands_folder = "/home/anna/Documents/benchmark/inputs/mcl1/ligands"
lig_name = "lig_60"

protocol = pipeline_protocol()
protocol.ligand_forcefield("sage")
print(protocol.ligand_forcefield())
# load ligand, these should already be in the correct position
try: # sdf first
    ligand = BSS.IO.readMolecules(f"{ligands_folder}/{lig_name}.sdf")[0]
except: # mol2 if sdf is not available
    ligand = BSS.IO.readMolecules(f"{ligands_folder}/{lig_name}.mol2")[0]

# paramaterise the ligand
print(f"Parameterising {lig_name}...")
lig_p = ligprep.lig_paramaterise(ligand, protocol.ligand_forcefield())

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

work_dir = "aa"
kwargs = {}
# for presentation
work_dir = "path/to/perturbation"
BSS.Relative.FreeEnergy.analyse(work_dir, estimator="MBAR", method="alchemlyb", **kwargs)

