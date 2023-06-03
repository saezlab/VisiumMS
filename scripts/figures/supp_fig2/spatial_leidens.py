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
fig_path = 'figures/manuscript/supp_fig2/'
plt.rcParams['font.sans-serif'] = 'Arial'

# Read
adata = sc.read_h5ad('data/prc/visium/props_niches.h5ad')
adata = adata[adata.obs['leiden'] != '7']
meta = pd.read_csv(os.path.join('data','metadata.csv')).set_index('sample_id')

# Plot slides
for i, sample_id in enumerate(meta.index):
    slide = get_ann_slide(sample_id, adata)
    lesion = meta.loc[sample_id, 'lesion_type']
    fig, ax = plt.subplots(1, 1, figsize=(4, 4), facecolor='white', tight_layout=True, dpi=150)
    sc.pl.spatial(slide, color='leiden', size=1.5, return_fig=False, show=False, ax=ax, frameon=False, legend_loc=False)
    ax.set_title('{0} {1}'.format(lesion, sample_id))
    os.makedirs(fig_path, exist_ok=True)
    fig.savefig(os.path.join(fig_path, 'niches_{0}.pdf'.format(sample_id)), bbox_inches='tight')
