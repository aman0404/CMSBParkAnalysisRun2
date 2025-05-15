import numpy as np
import dask.dataframe as dd
from dask.distributed import Client, LocalCluster
import glob
import math
from itertools import repeat
import mplhep as hep
import time
import matplotlib.pyplot as plt
import copy
import random as rand
import pandas as pd
from functools import reduce

from ROOT import TCanvas, gStyle, TH1F,TH2F, TLegend, TFile, RooRealVar, RooCBShape, gROOT, RooFit, RooAddPdf, RooDataHist, RooDataSet
from ROOT import RooArgList

from array import array
import sys
import argparse


cluster = LocalCluster()
cluster.scale(20) # create 4 local workers
client = Client(cluster)

load_fields = [
        "Sphericity",
        "DisplacedTrackN",
        "SumPt",
    ]

paths = "/cms/kaur/output/*.parquet"

data_files = glob.glob(paths)

df_temp = dd.read_parquet(data_files)

df_data   = df_temp[load_fields]

print("computation complete")

binning = (
    [j for j in range(0, 25, 1)]
    + [25]
)

df_data = df_data[(df_data["Sphericity"] > 0.8) & (df_data["DisplacedTrackN"] >= 9)]

data_sumpT = df_data["SumPt"].compute().values

h_data = TH1F("h_data", "h_data", len(binning)-1, array('d', binning))

for i in range(len(data_sumpT)):
    h_data.Fill(data_sumpT[i])


file2 = TFile("out_plots.root","RECREATE")
file2.cd()
h_data.Write()
file2.Close()

