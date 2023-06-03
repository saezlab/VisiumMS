import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Define path
fig_path = 'figures/manuscript/fig3/'
fig_name = 'umap_states.pdf'
plt.rcParams['font.sans-serif'] = 'Arial'

adata = sc.read_h5ad('data/prc/sc/microglia.h5ad')
adata = adata[~np.isin(adata.obs['leiden'], ['5'])]
adata.uns['log1p']["base"] = None

fig, axes = plt.subplots(1,1, figsize=(5, 4), facecolor='white', tight_layout=True, dpi=150)
sc.pl.umap(adata, color='leiden', frameon=False, ax=axes, return_fig=False, show=False, size=50)
axes.set_title('')

# Write
os.makedirs(fig_path, exist_ok=True)
fig.savefig(os.path.join(fig_path, fig_name), bbox_inches='tight')
