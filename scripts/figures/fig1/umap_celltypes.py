import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

# Defina path
fig_path = 'figures/manuscript/fig1/'
fig_name = 'umap_celltypes.pdf'

# Read data
adata = sc.read_h5ad('data/prc/sc/annotated.h5ad')
plt.rcParams['font.sans-serif'] = 'Arial'

# Plot
fig, ax = plt.subplots(1,1, figsize=(8, 8), dpi=150, facecolor='white')
sc.pl.umap(adata, color='leiden', ax=ax, size=25, frameon=False, legend_loc='on data', return_fig=False, show=False)
ax.set_title('')
os.makedirs(fig_path, exist_ok=True)
fig.savefig(os.path.join(fig_path, fig_name), bbox_inches='tight')
