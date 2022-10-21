import pandas as pd
import numpy as np
import scanpy as sc
import os
import seaborn as sns
import matplotlib.pyplot as plt


# Defina path
fig_path = 'figures/manuscript/fig2/'
fig_name = 'umap_niches.pdf'
plt.rcParams['font.sans-serif'] = 'Arial'

# Read
adata = sc.read_h5ad('data/prc/visium/props_niches.h5ad')
adata = adata[adata.obs['leiden'] != '7']

# Plot
fig, ax = plt.subplots(1,1, figsize=(6, 6), dpi=150, facecolor='white')
sc.pl.umap(adata, color='leiden', ax=ax, size=50, frameon=False, return_fig=False, show=False)
ax.set_title('')
os.makedirs(fig_path, exist_ok=True)
fig.savefig(os.path.join(fig_path, fig_name), bbox_inches='tight')
