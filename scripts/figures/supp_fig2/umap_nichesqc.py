import pandas as pd
import numpy as np
import scanpy as sc
import os
import seaborn as sns
import matplotlib.pyplot as plt
import decoupler as dc


# Defina path
fig_path = 'figures/manuscript/supp_fig2/'
fig_name = 'umap_nichesqc.pdf'
plt.rcParams['font.sans-serif'] = 'Arial'

# Read
adata = sc.read_h5ad('data/prc/visium/props_niches.h5ad')
adata = adata[adata.obs['leiden'] != '7']

# Add palette
palette = pd.read_csv('data/cond_palette.csv')
palette['rgb'] = [(r, g, b) for r, g, b in zip(palette['r'], palette['g'], palette['b'])]
palette = {k: v for k,v in zip(palette['cond'], palette['rgb'])}

# Plot
sc.set_figure_params(dpi=150, fontsize=14, figsize=(7, 7), facecolor='white')
acts = dc.get_acts(adata, 'props')
fig = sc.pl.umap(acts, color=list(acts.var_names) + ['lesion_type'], size=50, frameon=False, ncols=6, return_fig=True, palette=palette)
os.makedirs(fig_path, exist_ok=True)
fig.savefig(os.path.join(fig_path, fig_name), bbox_inches='tight')
