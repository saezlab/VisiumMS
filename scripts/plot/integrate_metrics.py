import scanpy as sc
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import argparse
import os

"""
Script to plot different metrics after integration.
"""

# Read command line and set args
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_path', help='Input path to object', required=True)
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
fig = plt.figure(figsize=(12,6), dpi=150, tight_layout=True, facecolor='white')
fig.suptitle('HVG filtering QC plots', fontsize=11)
gs = fig.add_gridspec(2, 3)

ax = fig.add_subplot(gs[0,0])
sc.pl.umap(adata, color='sample_id', ax=ax, frameon=False, return_fig=False, show=False)

ax = fig.add_subplot(gs[0,1])
sc.pl.umap(adata, color='lesion_type', ax=ax, frameon=False, return_fig=False, show=False)

ax = fig.add_subplot(gs[0,2])
sc.pl.umap(adata, color='doublet_score', ax=ax, frameon=False, return_fig=False, show=False)

ax = fig.add_subplot(gs[1,0])
sc.pl.umap(adata, color='diss_score', ax=ax, frameon=False, return_fig=False, show=False)

ax = fig.add_subplot(gs[1,1])
sc.pl.umap(adata, color='n_genes_by_counts', ax=ax, frameon=False, return_fig=False, show=False)

ax = fig.add_subplot(gs[1,2])
sc.pl.umap(adata, color='pct_counts_mt', ax=ax, frameon=False, return_fig=False, show=False)

# Save
if version is not None:
    fname = 'integrated_summary_{0}.png'.format(version)
else:
    fname = 'integrated_summary.png'
fig.savefig(os.path.join(output_path, fname))