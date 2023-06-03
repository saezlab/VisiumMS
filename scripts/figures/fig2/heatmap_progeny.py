import pandas as pd
import numpy as np
import scanpy as sc
import os
import matplotlib.pyplot as plt
import decoupler as dc
from scipy import stats
from scipy.stats import ranksums
import seaborn as sns


# Define path
fig_path = 'figures/manuscript/fig2/'
fig_name = 'heatmap_progeny.pdf'
plt.rcParams['font.sans-serif'] = 'Arial'

# Read
adata = sc.read_h5ad('data/prc/visium/props_niches.h5ad')
adata = adata[adata.obs['leiden'] != '7']

# Plot
acts = dc.get_acts(adata, 'progeny')
mean_acts = dc.summarize_acts(acts, groupby='leiden', min_std=0)
scaled_mean_acts = mean_acts.copy()
scaled_mean_acts.loc[:, :] = stats.zscore(mean_acts.values, axis=0, ddof=1)

# Compute ranksums
leidens = np.unique(acts.obs['leiden'].values)
pathways = acts.var_names
pvals = []
for leiden in leidens:
    row = []
    msk = acts.obs['leiden'] == leiden
    for pathway in pathways:
        w, p = ranksums(acts[:, pathway][msk].X, acts[:, pathway][~msk].X, alternative='greater')
        row.append(p[0])
    row = dc.p_adjust_fdr(row)
    pvals.append(row)
pvals = pd.DataFrame(pvals, columns=pathways, index=leidens)

# Transform to asterisks
pvals.loc[:, :] = np.where(pvals.values < 0.05, '*', '')

# Define color map
cmap = plt.get_cmap('coolwarm').copy()
cmap.set_bad(color='gray')

# Plot
fig, ax = plt.subplots(1, 1, figsize=(4, 4), facecolor='white', dpi=125)
htm = sns.heatmap(scaled_mean_acts.T, cmap=cmap, square=True, center=0, vmax=1, vmin=-1, ax=ax, cbar_kws={"shrink": .4, "aspect": 5},
                  annot=pvals.T.values.astype('U'), fmt='', annot_kws={'fontweight': 'black', 'color': 'black'})
i = 0
for _, spine in htm.spines.items():
    if i % 2 == 0:
        spine.set_visible(True)
    i += 1

os.makedirs(fig_path, exist_ok=True)
fig.savefig(os.path.join(fig_path, fig_name), bbox_inches='tight')
