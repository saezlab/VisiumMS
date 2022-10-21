import pandas as pd
import numpy as np
import scanpy as sc

from scipy import stats

import seaborn as sns
import matplotlib.pyplot as plt
import decoupler as dc
import os


# Defina path
fig_path = 'figures/manuscript/supp_fig3/'
fig_name = 'spatial_microglia_examples.pdf'
plt.rcParams['font.sans-serif'] = 'Arial'

# Read
deg_net = pd.read_csv('figures/manuscript/fig3/markers.csv').groupby('group').head(100)
slides = sc.read_h5ad('data/prc/visium/props_niches.h5ad')

# Remove spots without Microglia
msk = slides.obsm['props']['Microglia'].values > 0.1111
slides = slides[msk]

# Enrich by cellstate
dc.run_ulm(slides, deg_net, source='group', target='names', weight=None, use_raw=False)
slides.obsm['states'] = slides.obsm['ulm_estimate'].copy()


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
    msk = adata.obs['sample_id'] == sample_id
    obs = adata[msk].obs.copy()
    progeny = adata[msk].obsm['progeny'].copy()
    states = adata[msk].obsm['states'].copy()
    
    index = [''.join([idx.split('-')[0], '-1']) for idx in obs.index]
    obs.index = index
    progeny.index = index
    states.index = index
    
    slide = read_slide(sample_id)
    inter = obs.index.intersection(slide.obs.index)
    obs = obs.loc[inter]
    slide = slide[inter]
    progeny = progeny.loc[inter]
    states = states.loc[inter]
    slide.obs = obs
    slide.obsm['progeny'] = progeny
    slide.obsm['states'] = states
    return slide

def plot_slide(sample_id, axes, slides):

    slide = get_ann_slide(sample_id, slides)
    slide = slide[slide.obsm['props']['Microglia'] > 0.1111]

    ax = axes[0]
    props = dc.get_acts(slide, 'props')
    sc.pl.spatial(props, color='Macrophages_f', return_fig=False, show=False, ax=ax, size=1.6, frameon=False, cmap='magma')

    ax = axes[1]
    props = dc.get_acts(slide, 'states')
    sc.pl.spatial(props, color='2', return_fig=False, show=False, ax=ax, size=1.6, frameon=False, cmap='coolwarm')

    ax = axes[2]
    props = dc.get_acts(slide, 'states')
    sc.pl.spatial(props, color='3', return_fig=False, show=False, ax=ax, size=1.6, frameon=False, cmap='coolwarm')

# Write
fig, axes = plt.subplots(1, 9, figsize=(3.5*9, 3), facecolor='white', tight_layout=True, dpi=150)
ax = axes[0:3]
plot_slide('CO40', ax, slides)

ax = axes[3:6]
plot_slide('MS371', ax, slides)

ax = axes[6:9]
plot_slide('MS497I', ax, slides)

# Write
os.makedirs(fig_path, exist_ok=True)
fig.savefig(os.path.join(fig_path, fig_name), bbox_inches='tight')