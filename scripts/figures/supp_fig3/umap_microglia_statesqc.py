import pandas as pd
import numpy as np
import scanpy as sc
import os
import seaborn as sns
import matplotlib.pyplot as plt
import decoupler as dc


# Defina path
fig_path = 'figures/manuscript/supp_fig3/'
fig_name = 'umap_microglia_statesqc.pdf'
plt.rcParams['font.sans-serif'] = 'Arial'

# Read
adata = sc.read_h5ad('data/prc/sc/microglia.h5ad')
adata = adata[~np.isin(adata.obs['leiden'], ['5'])]

# Add palette
palette = pd.read_csv('data/cond_palette.csv')
palette['rgb'] = [(r, g, b) for r, g, b in zip(palette['r'], palette['g'], palette['b'])]
palette = {k: v for k,v in zip(palette['cond'], palette['rgb'])}

# Plot
fig, axes = plt.subplots(1,2, figsize=(10, 4), facecolor='white', tight_layout=True, dpi=150)
sc.pl.umap(adata, color='lesion_type', frameon=False, ax=axes[0], return_fig=False, show=False, size=50, palette=palette)
sc.pl.umap(adata, color='sample_id', frameon=False, ax=axes[1], return_fig=False, show=False, size=50)
axes[0].set_title('')
axes[1].set_title('')

# Write
os.makedirs(fig_path, exist_ok=True)
fig.savefig(os.path.join(fig_path, fig_name), bbox_inches='tight')
