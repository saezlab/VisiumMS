import scanpy as sc
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from plotting import plot_ngene_diff, plot_hvg_nbatches, plot_ngenes_vs_counts, plot_sorted_rank
import argparse
import os

"""
Script to plot different QC metrics after merging the data.
"""

# Read command line and set args
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_path', help='Input path to merged object', required=True)
parser.add_argument('-v', '--version', help='Version of the merging', required=False)
parser.add_argument('-o', '--output_dir', help='Output directory', required=True)
args = vars(parser.parse_args())

input_path = args['input_path']
version = args['version']
output_path = args['output_dir']
###############################

# Read merged object
adata = sc.read_h5ad(input_path)

# Plot HVG filtering QC plots
fig = plt.figure(figsize=(12,9), dpi=150, tight_layout=True, facecolor='white')
fig.suptitle('HVG filtering QC plots', fontsize=11)
gs = fig.add_gridspec(3, 3)

ax = fig.add_subplot(gs[0,0])
sc.pl.umap(adata, color='sample_id', ax=ax, frameon=False, return_fig=False, show=False)

ax = fig.add_subplot(gs[0,1])
sc.pl.umap(adata, color='lesion_type', ax=ax, frameon=False, return_fig=False, show=False)

ax = fig.add_subplot(gs[1,0])
sc.pl.umap(adata, color='doublet_score', ax=ax, frameon=False, return_fig=False, show=False)

ax = fig.add_subplot(gs[1,1])
sc.pl.umap(adata, color='diss_score', ax=ax, frameon=False, return_fig=False, show=False)

ax = fig.add_subplot(gs[0,2])
plot_hvg_nbatches(adata.var, ax)

ax = fig.add_subplot(gs[1,2])
plot_ngenes_vs_counts(adata.obs, ax, gene_thr=np.nan)

ax = fig.add_subplot(gs[2,:])
lst_samples = adata.obs['sample_id'].cat.categories
sc.pl.violin(adata, 'n_genes_by_counts', ax=ax, rotation=45, 
             groupby='sample_id', stripplot=False, show=False, order=lst_samples)

# Save
if version is not None:
    fname = 'merged_summary_{0}.png'.format(version)
else:
    fname = 'merged_summary.png'
fig.savefig(os.path.join(output_path, fname))