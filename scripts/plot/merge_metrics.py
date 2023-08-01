
# python scripts/plot/merge_metrics.py --output cellbender
# python scripts/plot/merge_metrics.py --output cellranger

import scanpy as sc
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from plotting import plot_ngene_diff, plot_hvg_nbatches, plot_ngenes_vs_counts, plot_sorted_rank
import argparse
from pathlib import Path
import os

"""
Script to plot different QC metrics after merging the data.
"""

# add command line flag arguments to specify either "cellbender" or "cellranger" output
parser = argparse.ArgumentParser()
parser.add_argument("--output", type=str, required=True)
args = parser.parse_args()

# set up relative paths within the project
current_folder = Path(__file__).parent
output_dir = current_folder / ".." / ".." / "data" / "prc" / "sc"
if args.output == "cellbender":
    input_path = current_folder / ".." / ".." / "data" / "prc" / "sc" / "cellbender_merged.h5ad"
    output_dir = current_folder / ".." / ".." / "out" / "cellbender_merge"
elif args.output == "cellranger":
    input_path = current_folder / ".." / ".." / "data" / "prc" / "sc" / "cellranger_merged.h5ad"
    output_dir = current_folder / ".." / ".." / "out" / "cellranger_merge"
else:
    raise ValueError("output must be either 'cellbender' or 'cellranger'")
output_dir.mkdir(parents=True, exist_ok=True)

# verbose
print("input_path: ", input_path)
print("output_dir: ", output_dir)

version = "0"

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
fig.savefig(os.path.join(output_dir, fname))