import scanpy as sc
import numpy as np
import pandas as pd
import decoupler as dc
from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt
import os


# Define path
fig_path = 'figures/manuscript/fig2/'
fig_name = 'dendo_vm.pdf'
plt.rcParams['font.sans-serif'] = 'Arial'

# Read
adata = sc.read_h5ad('data/prc/visium/props_niches.h5ad')
adata = adata[adata.obs['leiden'] != '7']

# Get pseudo-bulk profile
padata = dc.get_pseudobulk(adata, sample_col='sample_id', layer='counts', groups_col=None, min_prop=0.2, min_smpls=2)

# Normalize
sc.pp.normalize_total(padata, target_sum=1e4)
sc.pp.log1p(padata)

# Run dendogram
linked = linkage(padata.X, 'single', metric='cosine')

# Add palette
meta = pd.read_csv(os.path.join('data','metadata.csv')).set_index('sample_id')
palette = pd.read_csv('data/cond_palette.csv')
palette['rgb'] = [(r, g, b) for r, g, b in zip(palette['r'], palette['g'], palette['b'])]
palette = {k: v for k,v in zip(palette['cond'], palette['rgb'])}
s_palette = {k: palette[meta.loc[k, 'lesion_type']] for k in meta.index}

# Plot
fig, ax = plt.subplots(1,1, figsize=(4,4), facecolor='white', tight_layout=True, dpi=150)
labelList = padata.obs['sample_id'].values
dendrogram(linked, orientation='left', labels=labelList, distance_sort='descending', show_leaf_counts=True, ax=ax, link_color_func=lambda k: 'grey')
ylbls = ax.get_ymajorticklabels()
for lbl in ylbls:
    lbl.set_color(s_palette[lbl.get_text()])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)

# Write
os.makedirs(fig_path, exist_ok=True)
fig.savefig(os.path.join(fig_path, fig_name), bbox_inches='tight')
