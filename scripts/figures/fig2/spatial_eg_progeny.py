import pandas as pd
import numpy as np
import scanpy as sc
import os
import matplotlib.pyplot as plt
import decoupler as dc


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


# Define path
fig_path = 'figures/manuscript/fig2/'
fig_name = 'spatial_eg_progeny.pdf'
plt.rcParams['font.sans-serif'] = 'Arial'

# Define samples
samples = {
    'MS371': ([750, 4000, 3000, 1000], 'OPC', 'Astros_c', 'PI3K'),
    'MS197U': ([5000, 8500, 7500, 3500], 'Stroma', 'Endothelia', 'TGFb'),
}

# Read
adata = sc.read_h5ad('data/prc/visium/props_niches.h5ad')
adata = adata[adata.obs['leiden'] != '7']
meta = pd.read_csv(os.path.join('data','metadata.csv')).set_index('sample_id')


# Define samples
samples = {
    'CO74': ([2000, 5000, 3000, 4500], 'Oligos', 'OPC', 'Androgen'),
    'MS371': ([2000, 5000, 1250, 2750], 'Microglia', 'Macrophages_f', 'TNFa'),
    'MS377T': ([2000, 5000, 1500, 3000], 'Endothelia', 'Astros', 'WNT'),
    'MS549T': ([2500, 8500, 1900, 4900], 'Stroma', 'Endothelia', 'TGFb'),
}

# Plot
fig, axes = plt.subplots(4, 5, figsize=(16, 6), facecolor='white', tight_layout=True, dpi=150)
for i, sample_id in enumerate(samples):
    slide = get_ann_slide(sample_id, adata)
    crop_coord, ctype_a, ctype_b, pathway = samples[sample_id]
    lesion = meta.loc[sample_id, 'lesion_type']
    i_axes = axes[i]
    sc.pl.spatial(slide, color=None, size=1.5, return_fig=False, show=False, ax=i_axes[0], frameon=False, crop_coord=crop_coord, legend_loc=False)
    sc.pl.spatial(slide, color='leiden', size=1.5, return_fig=False, show=False, ax=i_axes[1], frameon=False, crop_coord=crop_coord, legend_loc=False)
    i_axes[0].set_title('{0} {1}'.format(lesion, sample_id))
    acts = dc.get_acts(slide, 'props')
    sc.pl.spatial(acts, color=ctype_a, cmap='magma', size=1.5, return_fig=True, show=False, ax=i_axes[2], frameon=False, crop_coord=crop_coord)
    sc.pl.spatial(acts, color=ctype_b, cmap='magma', size=1.5, return_fig=True, show=False, ax=i_axes[3], frameon=False, crop_coord=crop_coord)
    acts = dc.get_acts(slide, 'progeny')
    sc.pl.spatial(acts, color=pathway, cmap='coolwarm', size=1.5, return_fig=True, show=False, ax=i_axes[4], frameon=False, crop_coord=crop_coord)

os.makedirs(fig_path, exist_ok=True)
fig.savefig(os.path.join(fig_path, fig_name), bbox_inches='tight')