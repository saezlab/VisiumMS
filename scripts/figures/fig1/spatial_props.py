import scanpy as sc
import pandas as pd
import numpy as np
import decoupler as dc
import os
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


def read_slide(sample_id):

    # Read rna-seq
    slide = sc.read_visium('data/raw/visium/{0}/outs/'.format(sample_id))
    slide.var_names_make_unique()
    sc.pp.filter_genes(slide, min_cells=3)
    sc.pp.filter_cells(slide, min_genes=200)

    # Store raw counts
    slide.layers['counts'] = slide.X

    # Normalize
    sc.pp.normalize_total(slide, target_sum=1e4)
    sc.pp.log1p(slide)

    # Read props
    props = pd.read_csv('data/prc/visium/{0}/cell_props.csv'.format(sample_id), index_col=0)
    inter = slide.obs.index.intersection(props.index)
    slide.obsm['props'] = props.loc[inter]

    return slide

# Defina path
fig_path = 'figures/manuscript/fig1/'
fig_name = 'spatial_props.pdf'

# Plot
plt.rcParams['font.sans-serif'] = 'Arial'
sample_ids = ['CO40', 'MS371', 'MS377T', 'MS497I']
cell_types = [None, 'Astros', 'Oligos', 'Microglia', 'Macrophages_f']
lesions = ['Control', 'Acute', 'Chronic Active', 'Chronic Inactive']

fig, axes = plt.subplots(4, 5, facecolor='white', tight_layout=True, figsize=(15, 10))
for i, sample_id in enumerate(sample_ids):
    ax = axes[i]
    slide = read_slide(sample_id)
    acts = dc.get_acts(slide, 'props')
    lesion = lesions[i]
    for j, cell_type in enumerate(cell_types):
        sc.pl.spatial(acts, color=cell_type, size=1.5, return_fig=False, show=False, ax=ax[j],
                      frameon=False, vmin=0, vmax=0.8, cmap='magma')
        if i != 0:
            ax[j].set_title('')

os.makedirs(fig_path, exist_ok=True)
fig.savefig(os.path.join(fig_path, fig_name), bbox_inches='tight')
