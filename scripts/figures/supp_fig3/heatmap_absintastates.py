import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import scipy
import seaborn as sns

# Define path
fig_path = 'figures/manuscript/supp_fig3/'
fig_name = 'heatmap_absintastates.pdf'
plt.rcParams['font.sans-serif'] = 'Arial'

adata = sc.read_h5ad('data/prc/sc/microglia.h5ad')
adata = adata[~np.isin(adata.obs['leiden'], ['5'])]

iron = sc.read_h5ad('data/prc/sc/annotated.h5ad')
iron = iron[iron.obs['leiden'] == 'Macrophages_f']

adata = adata.concatenate(iron, join='outer', index_unique=None)

import decoupler as dc

deg = pd.read_csv('data/absinta2021/deg.csv', index_col=0)
names_dict = {
    1: 'MIMS-foamy',
    8: 'MIMS-iron',
    4: 'Cell-state 4'
}
deg = deg[np.isin(deg['cluster'], [1, 4, 8])]
deg['cluster'] = [names_dict[x] for x in deg['cluster']]
dc.run_ulm(adata, deg, source='cluster', target='gene', weight=None)
adata.obsm['absinta2021'] = adata.obsm['ulm_estimate'].copy()

acts = dc.get_acts(adata, 'absinta2021')
mean_acts = dc.summarize_acts(acts, 'leiden', obs=None, mode='mean', min_std=0)
mean_acts.loc[:, :] = scipy.stats.zscore(mean_acts.values, axis=0, ddof=1)

# Compute ranksums
a_groups = acts.var_names
l_groups = np.unique(acts.obs['leiden'].values)
pvals = []
for a_g in a_groups:
    row = []
    for l_g in l_groups:
        msk = acts.obs['leiden'] == l_g
        w, p = scipy.stats.ranksums(acts[:, a_g][msk].X, acts[:, a_g][~msk].X, alternative='greater')
        row.append(p[0])
    row = dc.p_adjust_fdr(row)
    pvals.append(row)
pvals = pd.DataFrame(pvals, columns=l_groups, index=a_groups).T
# Transform to asterisks
pvals.loc[:, :] = np.where(pvals.values < 0.05, '*', '')

# Define color map
cmap = plt.get_cmap('coolwarm').copy()
cmap.set_bad(color='gray')

mean_acts = mean_acts.T
pvals = pvals.T

fig, ax = plt.subplots(1, 1, figsize=(4, 4), facecolor='white', dpi=125)
htm = sns.heatmap(mean_acts, cmap=cmap, square=True, center=0, vmax=1, vmin=-1, ax=ax, cbar_kws={"shrink": .4, "aspect": 5},
                  annot=pvals.values.astype('U'), fmt='', annot_kws={'fontweight': 'black', 'color': 'black'})
i = 0
for _, spine in htm.spines.items():
    if i % 2 == 0:
        spine.set_visible(True)
    i += 1
    
os.makedirs(fig_path, exist_ok=True)
fig.savefig(os.path.join(fig_path, fig_name), bbox_inches='tight')
