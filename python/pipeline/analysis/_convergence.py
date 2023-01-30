# import libraries
import warnings
import BioSimSpace as BSS
from BioSimSpace import Units as _Units
import pickle
import os
import numpy as _np
import math as _math
import matplotlib.pyplot as plt
from scipy.stats import sem
import pickle
import csv

from ..utils import *
from . import plotting_engines

class plot_convergence():
    """class to plot convergence
    """

    def __init__(self, outputs_dir, perturbations=None, engines=None, file_ext=None, res_folder=None):

        # need the outputs directory
        self.outputs_dir = validate.folder_path(outputs_dir)
        self.engines = validate.engines(engines)

        self.file_ext = validate.string(file_ext)

        if not perturbations:
            raise ValueError("please include perturbations")
        else:
            self.perturbations = validate.is_list(perturbations)

        if not res_folder:
            self.res_folder = validate.folder_path(f"{outputs_dir}/convergence", create=True)
        else:
            self.res_folder = validate.folder_path(res_folder, create=True)

        self.set_colours()

    def set_colours(self, colour_dict=None):
        
        set_colour_dict = plotting_engines._set_colours(colour_dict)
        self.colours = set_colour_dict

        return set_colour_dict

    def plot_convergence_all(self):

        for pert in self.perturbations:
            lig_0 = pert.split("~")[0]
            lig_1 = pert.split("~")[1]

            self.plot_convergence_single(pert, engines=self.engines)


    def plot_convergence_single(self, perturbation, engines=None):

        for leg in [ 'free','bound']:
            
            plt.figure()
            lines = []
            
            for eng in engines:
                col = self.colours[eng]
                with open (f'{self.outputs_dir}/{eng}/{perturbation}/pickle/{leg}_pmf_{perturbation}_{eng}_{self.file_ext}.pickle', 'rb') as handle:
                    pmf_dict = pickle.load(handle)
                lines += plt.plot(0,0,c=col, label=eng)
                for repeat in pmf_dict:
                    pmf = pmf_dict[repeat]
                    x =[]
                    y=[]
                    yerr = []
                    for p in pmf:
                        x.append(p[0])
                        y.append(p[1]*(1/BSS.Units.Energy.kcal_per_mol))
                        yerr.append(p[2]*(1/BSS.Units.Energy.kcal_per_mol))
                    plt.errorbar(x,y,yerr=yerr,color=col)
            plt.xlim(xmin=0,xmax=1)
            plt.ylabel("Computed $\Delta$G$_{transformation}$ / kcal$\cdot$mol$^{-1}$")
            plt.xlabel("Lambda")
            labels = [l.get_label() for l in lines]
            plt.legend(lines, labels)
            plt.title(f"Convergence, {leg} for {perturbation}")
            plt.savefig(f'{self.res_folder}/{perturbation}_convergence_{leg}.png')

        # plotting delta delta G

        plt.figure()
        lines = []
        for eng in engines:
            # with open (f'/home/anna/Documents/benchmark/{prot}/outputs/{eng}/{perturbation}/bound_pmf_{perturbation}_{eng}.pickle', 'rb') as handle:
            with open (f'{self.outputs_dir}/{eng}/{perturbation}/pickle/bound_pmf_{perturbation}_{eng}_{self.file_ext}.pickle', 'rb') as handle:
                bound_pmf_dict = pickle.load(handle)
            # with open (f'/home/anna/Documents/benchmark/{prot}/outputs/{eng}/{perturbation}/free_pmf_{perturbation}_{eng}.pickle', 'rb') as handle:
            with open (f'{self.outputs_dir}/{eng}/{perturbation}/pickle/free_pmf_{perturbation}_{eng}_{self.file_ext}.pickle', 'rb') as handle:
                free_pmf_dict = pickle.load(handle)
            lines += plt.plot(0,0,c=self.colours[eng], label=eng)
            for repf,repb in zip(free_pmf_dict,bound_pmf_dict):
                bound_pmf = bound_pmf_dict[repb]
                free_pmf = free_pmf_dict[repf]
                x = []
                y = []
                yerr = []
                for pb,pf in zip(bound_pmf,free_pmf):
                    x.append(pb[0])
                    y.append((pb[1]*(1/BSS.Units.Energy.kcal_per_mol))-(pf[1]*(1/BSS.Units.Energy.kcal_per_mol)))
                    yerr.append((pb[2]*(1/BSS.Units.Energy.kcal_per_mol))+(pf[2]*(1/BSS.Units.Energy.kcal_per_mol)))
                plt.errorbar(x,y,yerr=yerr,color=self.colours[eng])
        plt.xlim(xmin=0,xmax=1)
        plt.ylabel("Computed $\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
        plt.xlabel("Lambda")
        labels = [l.get_label() for l in lines]
        plt.legend(lines, labels)
        plt.title(f"Convergence for {perturbation}")
        plt.savefig(f'{self.res_folder}/{perturbation}_convergence_deltadeltaG.png')


        # TODO plot average of each window

        plt.figure()
        lines = []
        for eng in engines:
            # with open (f'/home/anna/Documents/benchmark/{prot}/outputs/{eng}/{perturbation}/bound_pmf_{perturbation}_{eng}.pickle', 'rb') as handle:
            with open (f'{self.outputs_dir}/{eng}/{perturbation}/pickle/bound_pmf_{perturbation}_{eng}_{self.file_ext}.pickle', 'rb') as handle:
                bound_pmf_dict = pickle.load(handle)
            # with open (f'/home/anna/Documents/benchmark/{prot}/outputs/{eng}/{perturbation}/free_pmf_{perturbation}_{eng}.pickle', 'rb') as handle:
            with open (f'{self.outputs_dir}/{eng}/{perturbation}/pickle/free_pmf_{perturbation}_{eng}_{self.file_ext}.pickle', 'rb') as handle:
                free_pmf_dict = pickle.load(handle)
            lines += plt.plot(0,0,c=self.colours[eng], label=eng)

            # get the x (lambda) windows
            x = []
            for pf in free_pmf_dict[0]:
                x.append(pf[0])

            y_avg = []
            y_err = []

            for val in x:

                y = []
                yerr = []

            for repf,repb in zip(free_pmf_dict,bound_pmf_dict):
                bound_pmf = bound_pmf_dict[repb]
                free_pmf = free_pmf_dict[repf]
                
                for pb,pf in zip(bound_pmf,free_pmf):
                    y.append((pb[1]*(1/BSS.Units.Energy.kcal_per_mol))-(pf[1]*(1/BSS.Units.Energy.kcal_per_mol)))
                    yerr.append((pb[2]*(1/BSS.Units.Energy.kcal_per_mol))+(pf[2]*(1/BSS.Units.Energy.kcal_per_mol)))
                


                
            plt.errorbar(x,y_avg,yerr=yerr,color=self.colours[eng])
        
        plt.xlim(xmin=0,xmax=1)
        plt.ylabel("Computed $\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
        plt.xlabel("Lambda")
        labels = [l.get_label() for l in lines]
        plt.legend(lines, labels)
        plt.title(f"Convergence for {perturbation}")
        plt.savefig(f'{self.res_folder}/{perturbation}_convergence_deltadeltaG.png')