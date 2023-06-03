import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import os

# Defina path
fig_path = 'figures/manuscript/fig4/'
fig_name = 'dotplot_markers.pdf'
plt.rcParams['font.sans-serif'] = 'Arial'

# Read
adata = sc.read_h5ad('data/prc/sc/astros.h5ad')
adata = adata[~np.isin(adata.obs['leiden'], ['3', '6', '7'])]
adata.uns['log1p']["base"] = None

# Compute marker genes
sc.tl.rank_genes_groups(adata, 'leiden')

# Retrieve markers
df = sc.get.rank_genes_groups_df(adata, group=None)
df = df[(df['logfoldchanges'] > 1.) & (df['pvals'] < 0.05)].sort_values(['group', 'pvals', 'logfoldchanges'], ascending=[True, True, False])

# Define markers
markers = df.groupby('group').head(7).groupby('group')['names'].apply(lambda x: list(x))
markers.index = [idx.split('.')[0] for idx in markers.index]
markers = dict(markers)

# Actualize markers
markers['5'][1] = 'SPAG17'

# Plot
fig = sc.pl.dotplot(adata, markers, groupby='leiden', standard_scale='var', return_fig=True)

os.makedirs(fig_path, exist_ok=True)
fig.savefig(os.path.join(fig_path, fig_name), bbox_inches='tight')
df.to_csv(os.path.join(fig_path, 'markers.csv'), index=False)
