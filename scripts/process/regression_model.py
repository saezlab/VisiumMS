"""
Script to generate a regression model for cell2location from raw sn-seq data
"""

import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
from numpy.random import default_rng
import os

# this line forces theano to use the GPU and should go before importing cell2location
os.environ["THEANO_FLAGS"] = 'device=cuda0,floatX=float32,force_device=True'

import cell2location
import scvi

from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel

import argparse


# Read command line and set args
parser = argparse.ArgumentParser(prog='qc', description='Generates regression model from raw single nuc data')
parser.add_argument('-r', '--path_raw', help='Path to raw single nuclei adata', required=True)
parser.add_argument('-n', '--label_name', help='Label of cell type from raw data', required=True)
parser.add_argument('-s', '--sample_id', help='Label of sample id from raw data', required=True)
#parser.add_argument('-g', '--path_genes', help='Path to genes universe csv', required=True)
parser.add_argument('-p', '--perc_cells', help='Percentage of cells to subsample for speed and memory usage improvements', required=False, type=float, default=0.2)
parser.add_argument('-o', '--path_output', help='Path were to save model', required=True)
args = vars(parser.parse_args())

path_raw = args['path_raw']
label_name = args['label_name']
sample_id = args['sample_id']
#gene_uni = args['path_genes']
perc_cells = args['perc_cells']
path_output = args['path_output']

#Load
adata_raw = anndata.read_h5ad(path_raw)
adata_raw = adata_raw[~adata_raw.obs[label_name].isna(), :]


"""
Estimating expression signatures
"""
# Subsample cells from atlas for fast performance
rng = default_rng(seed=420)

t_cell_ids = []

# Iterate each cell type
for cell_type in adata_raw.obs[label_name].unique():
    
    # Select cells from a cell type
    msk = adata_raw.obs[label_name] == cell_type
    cell_ids = adata_raw.obs.index[msk]
    
    n_cells = int(np.ceil(perc_cells * len(cell_ids)))
    
    cell_ids = rng.choice(cell_ids, size=n_cells, replace=False)
    t_cell_ids.extend(cell_ids)
    
adata_raw = adata_raw[t_cell_ids]


"""
Basic QC
"""
# Filter by gene universe
#g_uni = pd.read_csv(gene_uni, index_col=0).index 
#inter = np.intersect1d(adata_raw.var.index, g_uni)
#adata_raw = adata_raw[:, inter].copy()

# Filter by cell2loc thresholds
selected = filter_genes(adata_raw, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
adata_raw = adata_raw[:, selected].copy()

"""
Prepare anndata for the regression model 
"""

RegressionModel.setup_anndata(adata=adata_raw,
                              # 10X reaction / sample / batch
                              batch_key=sample_id,
                              # cell type, covariate used for constructing signatures
                              labels_key=label_name
)

# Run regression model
mod = RegressionModel(adata_raw)

# Training 
mod.train(max_epochs=250, batch_size=2500, train_size=1, lr=0.002, use_gpu=True)


# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_raw = mod.export_posterior(
    adata_raw, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_raw.varm.keys():
    inf_aver = adata_raw.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_raw.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_raw.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_raw.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_raw.uns['mod']['factor_names']

inf_aver.to_csv(os.path.join(path_output,'inf_aver.csv'))

# Plot results

