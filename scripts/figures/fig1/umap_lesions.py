import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

# Defina path
fig_path = 'figures/manuscript/fig1/'
fig_name = 'umap_lesions.pdf'

# Read data
adata = sc.read_h5ad('data/prc/sc/annotated.h5ad')
plt.rcParams['font.sans-serif'] = 'Arial'

# Open palette
palette = pd.read_csv('data/cond_palette.csv')
palette['rgb'] = [(r, g, b) for r, g, b in zip(palette['r'], palette['g'], palette['b'])]
palette = {k: v for k,v in zip(palette['cond'], palette['rgb'])}
adata.obs['lesion_type'] = pd.Categorical(adata.obs['lesion_type'].values, categories=palette.keys(), ordered=True)

# Plot
fig, ax = plt.subplots(1,1, figsize=(8, 8), dpi=150, facecolor='white')
sc.pl.umap(adata, color='lesion_type', ax=ax, size=25, frameon=False, return_fig=False, show=False, palette=palette)
ax.set_title('')
os.makedirs(fig_path, exist_ok=True)
fig.savefig(os.path.join(fig_path, fig_name), bbox_inches='tight')