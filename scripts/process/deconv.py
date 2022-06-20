"""
Script to deconvolute visium slides into cell types using a pre-trained regression model.
"""

import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

# this line forces theano to use the GPU and should go before importing cell2location
os.environ["THEANO_FLAGS"] = 'device=cuda0,floatX=float32,force_device=True'

import cell2location
import scvi


import argparse


# Read command line and set args
parser = argparse.ArgumentParser(prog='qc', description='Generates regression model from raw single nuc data')
parser.add_argument('-s', '--slide_path', help='Path to slide adata', required=True)
parser.add_argument('-r', '--reg_path', help='Path to regression model', required=True)
parser.add_argument('-n', '--n_cells_spot', help='Number of cells per spot', required=False, type=float, default=30)
parser.add_argument('-a', '--d_apha', help='Regularization parammeter for technical effects', required=False, type=float, default=20)
parser.add_argument('-o', '--path_output', help='Path were to save deconvolution', required=True)
args = vars(parser.parse_args())

slide_path = args['slide_path']
reg_path = args['reg_path']
n_cells_spot = int(args['n_cells_spot'])
d_alpha = int(args['d_alpha'])
path_output = args['path_output']

# Read inputs
adata_vis = sc.read_visium(slide_path)
inf_aver = pd.read_csv(reg_path, index_col=0)

# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")

# create and train the model
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=n_cells_spot,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=d_alpha
)
mod.view_anndata_setup()

# Train
mod.train(max_epochs=30000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
          use_gpu=True)

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
)




