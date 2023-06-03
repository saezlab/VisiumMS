import pandas as pd
import numpy as np
import scanpy as sc
import os
import matplotlib.pyplot as plt
import decoupler as dc


# Define path
fig_path = 'figures/manuscript/fig2/'
fig_name = 'spatial_eg_props.pdf'
plt.rcParams['font.sans-serif'] = 'Arial'

def read_slide(sample_id):

    # Read rna-seq
    slide = sc.read_visium('data/raw/visium/{0}/outs/'.format(sample_id))
    slide.var_names_make_unique()
    sc.pp.filter_genes(slide, min_cells=3)
    sc.pp.filter_cells(slide, min_genes=200)

    # Store raw counts
    slide.layers['counts'] = slide.X.copy()

    # Normalize
    sc.pp.normalize_total(slide, target_sum=1e4)
    sc.pp.log1p(slide)

    # Read props
    props = pd.read_csv('data/prc/visium/{0}/cell_props.csv'.format(sample_id), index_col=0)
    inter = slide.obs.index.intersection(props.index)
    slide.obsm['props'] = props.loc[inter]

    return slide

meta = pd.read_csv(os.path.join('data','metadata.csv')).set_index('sample_id')

#### Define samples to plot
samples = ['CO74', 'MS497I'] # MS377I

# Plot
fig, axes = plt.subplots(2, 4, figsize=(12, 6), dpi=150, facecolor='white')
axes = axes.ravel()

slide = read_slide('CO74')
slide = dc.get_acts(slide, 'props')
ax = axes[0]
sc.pl.spatial(slide, color='Macrophages_f', size=1.6, ax=ax, return_fig=False, show=False, frameon=False, legend_loc=False, cmap='magma')
ax = axes[1]
sc.pl.spatial(slide, color='Oligos', size=1.6, ax=ax, return_fig=False, show=False, frameon=False, legend_loc=False, cmap='magma')

slide = read_slide('CO40')
slide = dc.get_acts(slide, 'props')
ax = axes[2]
sc.pl.spatial(slide, color='Stroma', size=1.6, ax=ax, return_fig=False, show=False, frameon=False, legend_loc=False, cmap='magma')
ax = axes[3]
sc.pl.spatial(slide, color='T-cells', size=1.6, ax=ax, return_fig=False, show=False, frameon=False, legend_loc=False, cmap='magma')


slide = read_slide('MS197U')
slide = dc.get_acts(slide, 'props')
ax = axes[4]
sc.pl.spatial(slide, color='Macrophages_f', size=1.6, ax=ax, return_fig=False, show=False, frameon=False, legend_loc=False, cmap='magma')
ax = axes[5]
sc.pl.spatial(slide, color='Oligos', size=1.6, ax=ax, return_fig=False, show=False, frameon=False, legend_loc=False, cmap='magma')

slide = read_slide('MS497I')
slide = dc.get_acts(slide, 'props')
ax = axes[6]
sc.pl.spatial(slide, color='Stroma', size=1.6, ax=ax, return_fig=False, show=False, frameon=False, legend_loc=False, cmap='magma')
ax = axes[7]
sc.pl.spatial(slide, color='T-cells', size=1.6, ax=ax, return_fig=False, show=False, frameon=False, legend_loc=False, cmap='magma')

#for i, sample in enumerate(samples):
#    i_ax = axes[i]
#    lesion = meta.loc[sample, 'lesion_type']
#    slide = read_slide(sample)
#    slide = dc.get_acts(slide, 'props')
#    sc.pl.spatial(slide, color='Macrophages_f', size=1.6, ax=i_ax[0], return_fig=False, show=False, frameon=False, legend_loc=False, cmap='magma')
#    sc.pl.spatial(slide, color='Oligos', size=1.6, ax=i_ax[1], return_fig=False, show=False, frameon=False, legend_loc=False, cmap='magma')
#    sc.pl.spatial(slide, color='Stroma', size=1.6, ax=i_ax[2], return_fig=False, show=False, frameon=False, legend_loc=False, cmap='magma')
#    sc.pl.spatial(slide, color='T-cells', size=1.6, ax=i_ax[3], return_fig=False, show=False, frameon=False, legend_loc=False, cmap='magma')

# Write
os.makedirs(fig_path, exist_ok=True)
fig.savefig(os.path.join(fig_path, fig_name), bbox_inches='tight')