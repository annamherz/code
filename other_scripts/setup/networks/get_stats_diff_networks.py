#!/usr/bin/python

# libraries

# import libraries
import BioSimSpace as BSS
import os
import sys
import glob
import csv
import numpy as np
import networkx as nx
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from sklearn.preprocessing import minmax_scale, MinMaxScaler
import itertools

import warnings

warnings.filterwarnings("ignore")

print("adding code to the pythonpath...")
code = "/home/anna/Documents/code/python"
if code not in sys.path:
    sys.path.insert(1, code)
import pipeline

from pipeline.prep import *
from pipeline.utils import *
from pipeline.setup import initialise_pipeline


import random
import math
import pandas as pd
import subprocess

fwf_path = (
    "/home/anna/Documents/september_2022_workshops/freenrgworkflows/networkanalysis"
)
if fwf_path not in sys.path:
    sys.path.insert(1, fwf_path)

import networkanalysis

from sklearn.metrics import mean_absolute_error as mae
from sklearn.metrics import mean_squared_error as mse

from scipy.stats import sem as sem
from scipy.stats import bootstrap, norm
from scipy.stats import spearmanr


# dictionary
best_fit_dict = {}

protein = sys.argv[1]
lr = sys.argv[2]  # lomap or rbfenn

exec_folder = f"/home/anna/Documents/benchmark/{protein}_benchmark/execution_model"

exper_val_dict = None

repeat_dict = {}
repeat = 0

if lr == "lomap":
    if protein == "tyk2":
        no_perts = 24
    elif protein == "mcl1":
        no_perts = 48
    elif protein == "p38":
        no_perts = 48
elif lr == "rbfenn":
    if protein == "tyk2":
        no_perts = 30
    elif protein == "mcl1":
        no_perts = 45
    elif protein == "p38":
        no_perts = 51

for file_ext in ["", "-a-optimal", "-d-optimal"]:
    repeat_dict[f"{lr}{file_ext}"] = {}

while repeat < 5001:
    command = (
        "/usr/bin/Rscript /home/anna/Documents/other_workflows/yang2020_optimal_designs/me/optimal_designs_prenorm.R %s %s %s"
        % (protein, lr, no_perts)
    )

    try:
        result = subprocess.run(command, shell=True, capture_output=True)
    except:
        print("failed to run R script.")
        continue

    for file_ext in ["-a-optimal", "-d-optimal"]:
        network_name = f"{lr}{file_ext}"

        print(f"{protein}, {network_name}, {repeat}")

        # rename files
        file = f"/home/anna/Documents/benchmark/{protein}_benchmark/execution_model/network_{network_name}.dat"

        commands = [
            "sed -i 's/LIGAND_1/lig0/g' %s" % (file),
            "sed -i 's/LIGAND_2/lig1/g' %s" % (file),
        ]

        for command in commands:
            try:
                result = subprocess.run(command, shell=True, capture_output=True)
            except:
                continue

        perts, ligs = pipeline.utils.get_info_network(f"{file}")

        if not exper_val_dict:
            exper_val_dict = pipeline.analysis.convert.yml_into_exper_dict(
                f"/home/anna/Documents/benchmark/inputs/experimental/{protein}.yml",
                temp=300,
            )
            normalised_exper_val_dict = pipeline.analysis.make_dict.exper_from_ligands(
                exper_val_dict, sorted(ligs), normalise=True
            )
        else:
            pass

        pert_dict = pipeline.analysis.make_dict.exper_from_perturbations(
            exper_val_dict, perts
        )

        # create fep pairwise ddG by randomly adding error to true value, which here is taken to be the experimental
        # sigmna2 fep is 1.0
        variance = 1
        pert_dict_fep = {}
        for pert in pert_dict:
            true_val = pert_dict[pert][0]
            # make value centred on experimental with variance sigma2fep
            rand_val = random.normalvariate(mu=true_val, sigma=math.sqrt(variance))
            pert_dict_fep[pert] = rand_val

        df = pd.DataFrame.from_dict(pert_dict_fep, orient="index").reset_index()
        df[["lig_0", "lig_1"]] = df["index"].str.split("~", expand=True)
        df = df.drop(columns="index")
        df = df.rename(columns={0: "freenrg"})
        df["error"] = 0.5
        df["engine"] = "SOMD"
        df = df[["lig_0", "lig_1", "freenrg", "error", "engine"]]

        # write into file for network analysis
        new_file_name = f"{protein}_{network_name}_fwf_file.csv"
        pd.DataFrame.to_csv(df, new_file_name, sep=",", index=False)

        nA = networkanalysis.NetworkAnalyser()

        nA.read_perturbations_pandas(
            new_file_name, comments="#", source="lig_0", target="lig_1"
        )

        computed_relative_DDGs = nA.freeEnergyInKcal
        freenrg_dict = (
            pipeline.analysis.make_dict.from_freenrgworkflows_network_analyser(
                computed_relative_DDGs
            )
        )

        x = [x[0] for x in normalised_exper_val_dict.values()]
        y = [y[0] for y in freenrg_dict.values()]
        # pipeline.analysis.stats_engines.compute_stats(x=x,y=y, statistic="MUE")
        dg_error = mae(x, y)

        x = [x[0] for x in pert_dict.values()]
        y = [y for y in pert_dict_fep.values()]
        # pipeline.analysis.stats_engines.compute_stats(x=x,y=y, statistic="MUE")
        ddg_error = mae(x, y)

        x = [x[0] for x in normalised_exper_val_dict.values()]
        y = [y[0] for y in freenrg_dict.values()]
        # pipeline.analysis.stats_engines.compute_stats(x=x,y=y, statistic="MUE")
        dg_error_mse = mse(x, y)

        x = [x[0] for x in pert_dict.values()]
        y = [y for y in pert_dict_fep.values()]
        # pipeline.analysis.stats_engines.compute_stats(x=x,y=y, statistic="MUE")
        ddg_error_mse = mse(x, y)

        x = [x[0] for x in normalised_exper_val_dict.values()]
        y = [y[0] for y in freenrg_dict.values()]
        # pipeline.analysis.stats_engines.compute_stats(x=x,y=y, statistic="MUE")
        coef, p = spearmanr(x, y)

        repeat_dict[network_name][repeat] = (
            dg_error,
            ddg_error,
            dg_error_mse,
            ddg_error_mse,
            coef,
            p,
        )

    repeat += 1

for file_ext in ["-a-optimal", "-d-optimal"]:  #
    # make a df

    df = pd.DataFrame.from_dict(
        repeat_dict[f"{lr}{file_ext}"],
        columns=[
            "mae_dG",
            "mae_ddG",
            "mse_dG",
            "mse_ddG",
            "spearman_coeff",
            "spearman_p",
        ],
        orient="index",
    )
    df.index.name = "repeat"
    df.to_csv(f"{exec_folder}/network_{network_name}_stats.csv")
