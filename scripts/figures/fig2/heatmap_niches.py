import pandas as pd
import numpy as np
import scanpy as sc
import os
import seaborn as sns
import matplotlib.pyplot as plt
from anndata import AnnData
from scipy.stats import ranksums
from scipy import stats
import decoupler as dc


# Defina path
fig_path = 'figures/manuscript/fig2/'
fig_name = 'heatmap_niches.pdf'
plt.rcParams['font.sans-serif'] = 'Arial'

# Read
adata = sc.read_h5ad('data/prc/visium/props_niches.h5ad')
adata = adata[adata.obs['leiden'] != '7']

# Remove cluster 7 due to low quality
subadata = adata[adata.obs['leiden'] != '7']

# Compute sample mean props
leidens = np.unique(subadata.obs['leiden'].values)
samples = np.unique(subadata.obs['sample_id'].values)
props = []
for leiden in leidens:
    for sample in samples:
        msk = (subadata.obs['leiden'] == leiden) & (subadata.obs['sample_id'] == sample)
        prop = subadata.obsm['props'].loc[msk].mean()
        prop = pd.DataFrame(prop, columns=[leiden + '_' + sample]).T
        props.append(prop)
props = pd.concat(props)
props = AnnData(props)
props.obs['leiden'] = [i.split('_')[0] for i in props.obs_names]
props.obs['sample_id'] = [i.split('_')[1] for i in props.obs_names] 

# Compute ranksums
ctypes = props.var_names
pvals = []
for leiden in leidens:
    row = []
    for ctype in ctypes:
        msk = props.obs['leiden'] == leiden
        w, p = ranksums(props[:, ctype][msk].X, props[:, ctype][~msk].X, alternative='greater')
        row.append(p[0])
    row = dc.p_adjust_fdr(row)
    pvals.append(row)
pvals = pd.DataFrame(pvals, columns=ctypes, index=leidens)

# Transform to asterisks
pvals.loc[:, :] = np.where(pvals.values < 0.05, '*', '')

# Compute mean props per niche
mean_props_leiden = []
for leiden in np.unique(subadata.obs['leiden'].values):
    msk = subadata.obs['leiden'] == leiden
    mean_props_leiden.append(pd.DataFrame(subadata.obsm['props'].loc[msk].mean(), columns=[leiden]).T)
mean_props_leiden = pd.concat(mean_props_leiden)

# Scale
scaled_mean_props_leiden = mean_props_leiden.copy()
scaled_mean_props_leiden.loc[:, :] = stats.zscore(scaled_mean_props_leiden.values, axis=0, ddof=1)

# Sort names
sorted_niches = ['0', '1', '2', '3', '4', '5', '6']
sorted_celltypes = np.sort(['Macrophages_f', 'Oligos', 'OPC',
                    'Oligos_d', 'Astros', 'Astros_c', 'B-cells', 'Microglia', 'Stroma', 'Endothelia', 'T-cells'])

scaled_mean_props_leiden = scaled_mean_props_leiden.loc[sorted_niches, sorted_celltypes]
pvals = pvals.loc[sorted_niches, sorted_celltypes]

# Define color map
cmap = plt.get_cmap('PiYG').copy()
cmap.set_bad(color='gray')

scaled_mean_props_leiden = scaled_mean_props_leiden.T
pvals = pvals.T

fig, ax = plt.subplots(1, 1, figsize=(4, 4), facecolor='white', dpi=125)
htm = sns.heatmap(scaled_mean_props_leiden, cmap=cmap, square=True, center=0, vmax=1, vmin=-1, ax=ax, cbar_kws={"shrink": .4, "aspect": 5},
                  annot=pvals.values.astype('U'), fmt='', annot_kws={'fontweight': 'black', 'color': 'black'})
i = 0
for _, spine in htm.spines.items():
    if i % 2 == 0:
        spine.set_visible(True)
    i += 1

os.makedirs(fig_path, exist_ok=True)
fig.savefig(os.path.join(fig_path, fig_name), bbox_inches='tight')
