# import libraries
import warnings
import BioSimSpace as BSS
from BioSimSpace import Units as _Units
import pickle
import os
import numpy as _np
import math as _math
from scipy.stats import sem
import pickle
import csv

from ..utils import *



# class plot_convergence():
#     """class to plot convergence
#     """

#     def __init__(self, work_dir):
#         # instantiate the class with the work directory

#         # need the output directory

#         self._work_dir = validate.folder_path(work_dir)
#         self._pickle_dir = validate.folder_path(f"{self._work_dir}/pickle", create=True)

#         # get the perturbation name and engine from the folder path
#         try:
#             self.perturbation = self._work_dir.split("/")[-1]
#             self.ligand_0 = self.perturbation.split("~")[0]
#             self.ligand_1 = self.perturbation.split("~")[1]
#             self.engine = validate.engine(
#                 self._work_dir.split("/")[-2].replace("_extracted", ""))
#         except:
#             warnings.warn(
#                 "was unable to get the perturbation name and engine from the file path.\n please add these to the class using self.perturbation and self.engine")



#     def plot_convergence_single(self, perturbation, engines=None):

#         if engines:
#             engines = self._plotting_engines(engines)
#         # if no engines provided, use the defaults that were set based on the analysis object
#         else:
#             engines = self.engines

#         tra = perturbation

#         # plot the convergence w time
#         prot = "tyk2"

#         col = self.colours[eng]

#         for leg in [ 'free','bound']:
#             plt.figure()
#             lines = []
#             for eng,col in zip(engine,colour):
#                 # with open (f'/home/anna/Documents/benchmark/{prot}/outputs/{eng}/{tra}/{leg}_pmf_{tra}_{eng}.pickle', 'rb') as handle:
#                 with open (f'/home/anna/Documents/benchmark/{prot}/pickles/{leg}_pmf_{tra}_{eng}.pickle', 'rb') as handle:
#                     pmf_dict = pickle.load(handle)
#                 lines += plt.plot(0,0,c=col, label=eng)
#                 for repeat in pmf_dict:
#                     pmf = pmf_dict[repeat]
#                     x =[]
#                     y=[]
#                     yerr = []
#                     for p in pmf:
#                         x.append(p[0])
#                         y.append(p[1]*(1/BSS.Units.Energy.kcal_per_mol))
#                         yerr.append(p[2]*(1/BSS.Units.Energy.kcal_per_mol))
#                     plt.errorbar(x,y,yerr=yerr,color=col, ecolor='black')
#             plt.xlim(xmin=0,xmax=1)
#             plt.ylabel("Computed $\Delta$G$_{transformation}$ / kcal$\cdot$mol$^{-1}$")
#             plt.xlabel("Lambda")
#             labels = [l.get_label() for l in lines]
#             plt.legend(lines, labels)
#             plt.title(f"Convergence, {leg} for {tra}")
#             plt.savefig(f'/home/anna/Documents/benchmark/{prot}/outputs/{eng}/{tra}/convergence_{leg}.png')

#         # plotting delta delta G

#         plt.figure()
#         lines = []
#         for eng in engines:
#             # with open (f'/home/anna/Documents/benchmark/{prot}/outputs/{eng}/{tra}/bound_pmf_{tra}_{eng}.pickle', 'rb') as handle:
#             with open (f'/home/anna/Documents/benchmark/{prot}/pickles/bound_pmf_{tra}_{eng}.pickle', 'rb') as handle:
#                 bound_pmf_dict = pickle.load(handle)
#             # with open (f'/home/anna/Documents/benchmark/{prot}/outputs/{eng}/{tra}/free_pmf_{tra}_{eng}.pickle', 'rb') as handle:
#             with open (f'/home/anna/Documents/benchmark/{prot}/pickles/free_pmf_{tra}_{eng}.pickle', 'rb') as handle:
#                 free_pmf_dict = pickle.load(handle)
#             lines += plt.plot(0,0,c=colour_dict[eng], label=eng)
#             for repf,repb in zip(free_pmf_dict,bound_pmf_dict):
#                 bound_pmf = bound_pmf_dict[repb]
#                 free_pmf = free_pmf_dict[repf]
#                 x = []
#                 y = []
#                 yerr = []
#                 for pb,pf in zip(bound_pmf,free_pmf):
#                     x.append(pb[0])
#                     y.append((pb[1]*(1/BSS.Units.Energy.kcal_per_mol))-(pf[1]*(1/BSS.Units.Energy.kcal_per_mol)))
#                     yerr.append((pb[2]*(1/BSS.Units.Energy.kcal_per_mol))+(pf[2]*(1/BSS.Units.Energy.kcal_per_mol)))
#                 plt.errorbar(x,y,yerr=yerr,color=colour_dict[eng])
#         plt.xlim(xmin=0,xmax=1)
#         plt.ylabel("Computed $\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
#         plt.xlabel("Lambda")
#         labels = [l.get_label() for l in lines]
#         plt.legend(lines, labels)
#         plt.title(f"Convergence for {tra}")
#         plt.savefig(f'/home/anna/Documents/benchmark/{prot}/outputs/{eng}/{tra}/convergence_deltadeltaG.png')

