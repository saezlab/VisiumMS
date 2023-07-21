import scanpy as sc
import scanpy.external as sce
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--input_path', required=True)
parser.add_argument('-p','--plot_path', required=True)
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

input_path = args['input_path']
plot_path = args['plot_path']
out_path = args['out_path']

# Read
adata = sc.read_h5ad(input_path)

# Run harmony
sce.pp.harmony_integrate(adata, key='Sample id', adjusted_basis='X_pca', max_iter_harmony=30)

# Run umap with updated connectivity
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# Plot
fig, axes = plt.subplots(2, 3, figsize=(12, 6), tight_layout=True, dpi=150)
axes = axes.ravel()

ax = axes[0]
sc.pl.umap(adata, color='doublet_score', ax=ax, return_fig=False, show=False, vmin=0, vmax=0.5)
ax.set_title('Doublet score')

ax = axes[1]
sc.pl.umap(adata, color='pct_counts_mt', ax=ax, return_fig=False, show=False, vmin=0, vmax=10)
ax.set_title('Fraction MT counts')

ax = axes[2]
sc.pl.umap(adata, color='Sample id', ax=ax, return_fig=False, show=False)

ax = axes[3]
adata.obs['log_n_genes_by_counts'] = np.log10(adata.obs['n_genes_by_counts'])
sc.pl.umap(adata, color='log_n_genes_by_counts', ax=ax, return_fig=False, show=False, vmin=0)
ax.set_title('Number of genes expr (log10)')
del adata.obs['log_n_genes_by_counts']

ax = axes[4]
adata.obs['log_total_counts'] = np.log10(adata.obs['total_counts'])
sc.pl.umap(adata, color='log_total_counts', ax=ax, return_fig=False, show=False, vmin=0)
ax.set_title('Total counts (log10)')
del adata.obs['log_total_counts']

ax = axes[5]
sc.pl.umap(adata, color='Lesion type', ax=ax, return_fig=False, show=False)

fig.savefig(plot_path, bbox_inches='tight')

# Save
adata.write(out_path)

