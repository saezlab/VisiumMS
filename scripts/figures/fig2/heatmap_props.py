import pandas as pd
import numpy as np
import scanpy as sc
import os
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats


# Define path
fig_path = 'figures/manuscript/fig2/'
fig_name = 'heatmap_props.pdf'
plt.rcParams['font.sans-serif'] = 'Arial'

# Read
adata = sc.read_h5ad('data/prc/visium/props_niches.h5ad')
adata = adata[adata.obs['leiden'] != '7']


def get_corrs_slide(slide):

    # Extract row and col names
    ctyps_props = slide.obsm['props'].columns.values
    ctyps_signs = slide.obsm['props'].columns.values

    # Init empty mats
    corrs = np.zeros((ctyps_signs.size, ctyps_props.size))
    pvals = np.zeros((ctyps_signs.size, ctyps_props.size))
    for i, c_a in enumerate(ctyps_signs):
        for j, c_b in enumerate(ctyps_props):
            
            if c_a == c_b:
                corrs[i, j], pvals[i, j] = np.nan, np.nan
            else:
                # Compute pearson
                corrs[i, j], pvals[i, j] = stats.pearsonr(slide.obsm['props'][c_a].values, slide.obsm['props'][c_b].values)

    # Transform to dfs
    corrs = pd.DataFrame(corrs, index=ctyps_signs, columns=ctyps_props)
    pvals = pd.DataFrame(pvals, index=ctyps_signs, columns=ctyps_props)

    # Flip to have same order as misty
    corrs = corrs.loc[np.flip(corrs.index)]
    pvals = pvals.loc[np.flip(pvals.index)]

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
    return out


# Read meta
meta = pd.read_csv(os.path.join('data','metadata.csv')).set_index('sample_id')

# Compute corrs
ids = meta.index.values
slides = []
corrs = []
pvals = []
for sample_id in ids:
    slide = adata[adata.obs['sample_id'] == sample_id]
    corr, pval = get_corrs_slide(slide)
    corrs.append(corr)
    pvals.append(pval)

corrs = aggregate(corrs)
pvals = aggregate(pvals)
corrs.to_csv(os.path.join(fig_path, 'props_corrs_table.csv'))
pvals.to_csv(os.path.join(fig_path, 'props_pvals_table.csv'))

# Transform to asterisks
pvals[np.isfinite(pvals)] = np.where((pvals < 0.10) & (np.abs(corrs) > 0.15), '*', '')
corrs = corrs.T
pvals = pvals.T


# Define color map
cmap = plt.get_cmap('PiYG').copy()
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

# Write
os.makedirs(fig_path, exist_ok=True)
fig.savefig(os.path.join(fig_path, fig_name), bbox_inches='tight')