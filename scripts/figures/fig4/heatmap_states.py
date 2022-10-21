import pandas as pd
import numpy as np
import scanpy as sc

from scipy import stats

import seaborn as sns
import matplotlib.pyplot as plt
import decoupler as dc
import os


# Defina path
fig_path = 'figures/manuscript/fig4/'
plt.rcParams['font.sans-serif'] = 'Arial'

# Read
deg_net = pd.read_csv('figures/manuscript/fig4/markers.csv').groupby('group').head(100)
slides = sc.read_h5ad('data/prc/visium/props_niches.h5ad')

# Remove spots without Astros
msk = (slides.obsm['props']['Astros'].values > 0.1111) | (slides.obsm['props']['Astros_c'].values > 0.1111)
slides = slides[msk]


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


def get_corrs_slide(slide, mode='progeny'):

    # Extract row and col names
    ctyps_props = slide.obsm['states'].columns.values
    ctyps_signs = slide.obsm[mode].columns.values

    # Init empty mats
    corrs = np.zeros((ctyps_signs.size, ctyps_props.size))
    pvals = np.zeros((ctyps_signs.size, ctyps_props.size))
    for i, c_a in enumerate(ctyps_signs):
        for j, c_b in enumerate(ctyps_props):
            
            if c_a == c_b:
                corrs[i, j], pvals[i, j] = np.nan, np.nan
            else:
                # Compute pearson
                corrs[i, j], pvals[i, j] = stats.pearsonr(slide.obsm[mode][c_a].values, slide.obsm['states'][c_b].values)

    # Transform to dfs
    corrs = pd.DataFrame(corrs, index=ctyps_signs, columns=ctyps_props)
    pvals = pd.DataFrame(pvals, index=ctyps_signs, columns=ctyps_props)

    # Flip to have same order as misty
    corrs = corrs.loc[np.flip(corrs.index)].T
    pvals = pvals.loc[np.flip(pvals.index)].T

    return corrs, pvals


def aggregate(lst):
    out = np.zeros(lst[0].shape)
    for i in range(lst[0].shape[0]):
        for j in range(lst[0].shape[1]):
            vals = np.array([lst[k].iloc[i, j] for k in range(len(lst))])
            if np.all(~np.isfinite(vals)):
                out[i, j] = np.nan
            else:
                out[i, j] = np.mean(vals[np.isfinite(vals)])
    out = pd.DataFrame(out, index=lst[0].index, columns=lst[0].columns)
    return out.T


def plot_heatmap_corr(corrs, pvals, cmap='coolwarm', p_thr=0.05, c_thr=0.15):
    
    pvals = pvals.copy()

    # Transform to asterisks
    pvals[np.isfinite(pvals)] = np.where((pvals < p_thr) & (np.abs(corrs) > c_thr), '*', '')

    # Define color map
    cmap = plt.get_cmap(cmap).copy()
    cmap.set_bad(color='gray')

    # Plot
    fig, ax = plt.subplots(1, 1, figsize=(4, 4), facecolor='white', dpi=125)
    htm = sns.heatmap(corrs, cmap=cmap, square=True, center=0, vmax=0.5, vmin=-0.5, ax=ax, cbar_kws={"shrink": .4, "aspect": 5},
                      annot=pvals.values.astype('U'), fmt='', annot_kws={'fontweight': 'black', 'color': 'black'})
    i = 0
    for _, spine in htm.spines.items():
        if i % 2 == 0:
            spine.set_visible(True)
        i += 1
    return fig


# Enrich by cellstate
dc.run_ulm(slides, deg_net, source='group', target='names', weight=None, use_raw=False)
slides.obsm['states'] = slides.obsm['ulm_estimate'].copy()

meta = pd.read_csv(os.path.join('data','metadata.csv')).set_index('sample_id')
ids = meta.index.values[meta['lesion_type'] == 'Chronic Active']

corrs_progeny = []
pvals_progeny = []
corrs_states = []
pvals_states = []
corrs_props = []
pvals_props = []
for sample_id in ids:
    slide = get_ann_slide(sample_id, slides)
    # Filter spots without astros
    msk = (slide.obsm['props']['Astros'].values > 0.1111) | (slide.obsm['props']['Astros_c'].values > 0.1111)
    slide = slide[msk]
    slide.obsm['props'] = slide.obsm['props'].drop(['Astros', 'Astros_c'], axis=1)
    corr, pval = get_corrs_slide(slide, 'progeny')
    corrs_progeny.append(corr)
    pvals_progeny.append(pval)
    
    corr, pval = get_corrs_slide(slide, 'states')
    corrs_states.append(corr)
    pvals_states.append(pval)

    corr, pval = get_corrs_slide(slide, 'props')
    corrs_props.append(corr)
    pvals_props.append(pval)

corrs_progeny = aggregate(corrs_progeny)
pvals_progeny = aggregate(pvals_progeny)
corrs_states = aggregate(corrs_states)
pvals_states = aggregate(pvals_states)
corrs_props = aggregate(corrs_props)
pvals_props = aggregate(pvals_props)

# Write
os.makedirs(fig_path, exist_ok=True)
fig = plot_heatmap_corr(corrs_props, pvals_props, cmap='PiYG')
fig.savefig(os.path.join(fig_path, 'heatmap_props.pdf'), bbox_inches='tight')

#plot_heatmap_corr(corrs_states, pvals_states)
fig = plot_heatmap_corr(corrs_progeny, pvals_progeny)
fig.savefig(os.path.join(fig_path, 'heatmap_progeny.pdf'), bbox_inches='tight')

sample_id = 'MS197U'

slide = get_ann_slide(sample_id, slides)
msk = (slide.obsm['props']['Astros'].values > 0.1111) | (slide.obsm['props']['Astros_c'].values > 0.1111)
slide = slide[msk]

# Write
fig, axes = plt.subplots(2, 3, figsize=(9, 6), facecolor='white', tight_layout=True, dpi=150)
axes = axes.ravel()

ax = axes[0]
props = dc.get_acts(slide, 'states')
sc.pl.spatial(props, color='2', return_fig=False, show=False, ax=ax, size=1.6, frameon=False, cmap='coolwarm')

ax = axes[1]
props = dc.get_acts(slide, 'props')
sc.pl.spatial(props, color='Macrophages_f', return_fig=False, show=False, ax=ax, size=1.6, frameon=False, cmap='magma')

ax = axes[2]
props = dc.get_acts(slide, 'progeny')
sc.pl.spatial(props, color='TNFa', return_fig=False, show=False, ax=ax, size=1.6, frameon=False, cmap='coolwarm')

ax = axes[3]
props = dc.get_acts(slide, 'states')
sc.pl.spatial(props, color='5', return_fig=False, show=False, ax=ax, size=1.6, frameon=False, cmap='coolwarm')

ax = axes[4]
props = dc.get_acts(slide, 'props')
sc.pl.spatial(props, color='Endothelia', return_fig=False, show=False, ax=ax, size=1.6, frameon=False, cmap='magma')

ax = axes[5]
props = dc.get_acts(slide, 'progeny')
sc.pl.spatial(props, color='PI3K', return_fig=False, show=False, ax=ax, size=1.6, frameon=False, cmap='coolwarm')

os.makedirs(fig_path, exist_ok=True)
fig.savefig(os.path.join(fig_path, 'spatial_eg_props.pdf'), bbox_inches='tight')
