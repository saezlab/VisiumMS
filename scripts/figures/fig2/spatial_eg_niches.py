import pandas as pd
import numpy as np
import scanpy as sc
import os
import matplotlib.pyplot as plt
import decoupler as dc


# Define path
fig_path = 'figures/manuscript/fig2/'
fig_name = 'spatial_eg_niches.pdf'
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


def get_ann_slide(sample_id, adata):
    msk = adata.obs['sample_id'].values == sample_id
    obs = adata[msk].obs.copy()
    progeny = adata.obsm['progeny'].loc[msk].copy()
    obs.index = [''.join([idx.split('-')[0], '-1']) for idx in obs.index]
    progeny.index = [''.join([idx.split('-')[0], '-1']) for idx in progeny.index]
    slide = read_slide(sample_id)
    inter = obs.index.intersection(slide.obs.index)
    obs = obs.loc[inter]
    progeny = progeny.loc[inter]
    slide = slide[inter]
    slide.obs = obs
    slide.obsm['progeny'] = progeny
    return slide



# Read
adata = sc.read_h5ad('data/prc/visium/props_niches.h5ad')
adata = adata[adata.obs['leiden'] != '7']
meta = pd.read_csv(os.path.join('data','metadata.csv')).set_index('sample_id')

# Plot
fig, axes = plt.subplots(2, 4, figsize=(10, 5), dpi=150, facecolor='white')

slide = get_ann_slide('CO74', adata)
ax = axes[0, 0]
sc.pl.spatial(slide, color=None, size=1.6, ax=ax, return_fig=False, show=False, frameon=False, legend_loc=False, cmap='magma')
ax = axes[1, 0]
sc.pl.spatial(slide, color='leiden', size=1.6, ax=ax, return_fig=False, show=False, frameon=False, legend_loc=False, cmap='magma', title='')

slide = get_ann_slide('MS371', adata)
ax = axes[0, 1]
sc.pl.spatial(slide, color=None, size=1.6, ax=ax, return_fig=False, show=False, frameon=False, legend_loc=False, cmap='magma')
ax = axes[1, 1]
sc.pl.spatial(slide, color='leiden', size=1.6, ax=ax, return_fig=False, show=False, frameon=False, legend_loc=False, cmap='magma', title='')

slide = get_ann_slide('MS377T', adata)
ax = axes[0, 2]
sc.pl.spatial(slide, color=None, size=1.6, ax=ax, return_fig=False, show=False, frameon=False, legend_loc=False, cmap='magma')
ax = axes[1, 2]
sc.pl.spatial(slide, color='leiden', size=1.6, ax=ax, return_fig=False, show=False, frameon=False, legend_loc=False, cmap='magma', title='')

slide = get_ann_slide('MS549T', adata)
ax = axes[0, 3]
sc.pl.spatial(slide, color=None, size=1.6, ax=ax, return_fig=False, show=False, frameon=False, legend_loc=False, cmap='magma')
ax = axes[1, 3]
sc.pl.spatial(slide, color='leiden', size=1.6, ax=ax, return_fig=False, show=False, frameon=False, legend_loc=False, cmap='magma', title='')

os.makedirs(fig_path, exist_ok=True)
fig.savefig(os.path.join(fig_path, fig_name), bbox_inches='tight')