import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
from numpy.random import default_rng
import os
import matplotlib.pyplot as plt
# this line forces theano to use the GPU and should go before importing cell2location
os.environ["THEANO_FLAGS"] = 'device=cuda0,floatX=float32,force_device=True'
import cell2location
import scvi
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--inp_path', required=True)
parser.add_argument('-e','--max_epochs', required=True)
parser.add_argument('-b','--batch_size', required=True)
parser.add_argument('-l','--lr', required=True)
parser.add_argument('-n','--num_samples', required=True)
parser.add_argument('-p','--plot_path', required=True)
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

inp_path = args['inp_path']
max_epochs = int(args['max_epochs'])
batch_size = int(args['batch_size'])
lr = float(args['lr'])
num_samples = int(args['num_samples'])
plot_path = args['plot_path']
out_path = args['out_path']

# Read annot sc data
print('Read adata')
adata = sc.read_h5ad(inp_path)
del adata.X
adata.X = adata.layers['counts'].copy()
del adata.layers['counts']

# Set up object
cell2location.models.RegressionModel.setup_anndata(
    adata=adata,
    batch_key='Sample id',
    labels_key='leiden'
)

# Run regression model
mod = RegressionModel(adata)
mod.view_anndata_setup()

# Training
print('Training:')
mod.train(
    max_epochs=max_epochs,
    batch_size=batch_size,
    train_size=1,
    lr=lr,
    use_gpu=True
)

# Generate gene signature
print('Generate signature')
adata = mod.export_posterior(
    adata,
    sample_kwargs={
        'num_samples': num_samples,
        'batch_size': batch_size,
        'use_gpu': True
    }
)

# Export estimated expression in each cluster
print('Export results')
if 'means_per_cluster_mu_fg' in adata.varm.keys():
    inf_aver = adata.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' for i in adata.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata.var[[f'means_per_cluster_mu_fg_{i}' for i in adata.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata.uns['mod']['factor_names']

# Write
print('Write')
inf_aver.to_csv(out_path)

# Plot
print('Plot')
fig, ax = plt.subplots(1, 1, figsize=(5, 3), dpi=150)
mod.plot_history(iter_start=20, ax=ax)
fig.suptitle('Regression model')
fig.savefig(plot_path, bbox_inches='tight')
