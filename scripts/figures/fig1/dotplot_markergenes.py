import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt


# Defina path
fig_path = 'figures/manuscript/fig1/'
fig_name = 'dotplot_markergenes.pdf'
markers_path = 'data/markers.csv'

# Read data
adata = sc.read_h5ad('data/prc/sc/annotated.h5ad')
plt.rcParams['font.sans-serif'] = 'Arial'

# Read marker genes
markers = pd.read_csv(markers_path)
markers = dict(markers.groupby('cell_type')['gene'].apply(list))
del markers['Astros_gm']
del markers['Neuron']

# Plot
fig, ax = plt.subplots(1,1, figsize=(9, 3), dpi=150, facecolor='white')
sc.pl.dotplot(adata, markers, 'leiden', dendrogram=False, ax=ax, return_fig=False, show=False, standard_scale='var')
os.makedirs(fig_path, exist_ok=True)
fig.savefig(os.path.join(fig_path, fig_name), bbox_inches='tight')
