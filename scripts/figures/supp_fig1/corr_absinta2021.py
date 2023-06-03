import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import decoupler as dc
import seaborn as sns
from scipy import stats


# Defina path
fig_path = 'figures/manuscript/supp_fig1/'
fig_name = 'corr_absinta2021.pdf'

# Read data
absinta2021 = sc.read_h5ad('data/absinta2021/absinta2021.h5ad')
lerma2022 = sc.read_h5ad('data/prc/sc/raw.h5ad')
plt.rcParams['font.sans-serif'] = 'Arial'

# Get pseudo-bulk profile
absinta2021 = dc.get_pseudobulk(absinta2021, sample_col='cell_type', groups_col=None, min_prop=0.2, min_smpls=0, use_raw=False)
lerma2022 = dc.get_pseudobulk(lerma2022, sample_col='leiden', groups_col=None, min_prop=0.2, min_smpls=0, use_raw=False)

# Filter genes by intersection
inter = absinta2021.var_names.intersection(lerma2022.var_names)
absinta2021 = absinta2021[:, inter]
lerma2022 = lerma2022[:, inter]

# Normalize
sc.pp.normalize_total(absinta2021, target_sum=1e4)
sc.pp.log1p(absinta2021)
sc.pp.normalize_total(lerma2022, target_sum=1e4)
sc.pp.log1p(lerma2022)

# Find unique labels
clusters_absinta2021 = np.unique(absinta2021.obs['cell_type'].values)
clusters_lerma2022 = np.unique(lerma2022.obs['leiden'].values)

# Compute corrs
corrs = []
pvals = []
for c_a in clusters_absinta2021:
    c_row = []
    p_row = []
    for c_b in clusters_lerma2022:
        x = absinta2021[c_a, :].X.ravel()
        y = lerma2022[c_b, :].X.ravel()
        c, p = stats.pearsonr(x, y)
        c_row.append(c)
        p_row.append(p)
    corrs.append(c_row)
    pvals.append(p_row)

corrs = pd.DataFrame(corrs, index=clusters_absinta2021, columns=clusters_lerma2022)
pvals = pd.DataFrame(pvals, index=clusters_absinta2021, columns=clusters_lerma2022)

# Sort labels to match
sorted_absinta2021 = np.array(['astrocytes', 'immune', 'lymphocytes', 'oligodendrocytes', 'opc', 'vascular_cells', 'neurons'])
sorted_lerma2021 = np.array(['Astros', 'Astros_c', 'Microglia', 'Macrophages_f', 
                             'B-cells', 'T-cells', 'Oligos', 'OPC', 'Oligos_d', 'Endothelia', 'Stroma'])
corrs = corrs.loc[sorted_absinta2021, sorted_lerma2021]
pvals = pvals.loc[sorted_absinta2021, sorted_lerma2021]
pvals.loc[:, :] = np.where((pvals < 0.05) & (np.abs(corrs) > 0.75), '*', '')

# Transpose
corrs = corrs.T
pvals = pvals.T

# Plot
cmap = plt.get_cmap('Blues').copy()
cmap.set_bad(color='gray')

fig, ax = plt.subplots(1, 1, figsize=(4, 4), facecolor='white', dpi=125)
htm = sns.heatmap(corrs, cmap=cmap, square=True, vmax=1, vmin=0, ax=ax, cbar_kws={"shrink": .4, "aspect": 5},
                  annot=pvals.values.astype('U'), fmt='', annot_kws={'fontweight': 'black', 'color': 'black'})
i = 0
for _, spine in htm.spines.items():
    if i % 2 == 0:
        spine.set_visible(True)
    i += 1
ax.set_xlabel('absinta2021')
ax.set_ylabel('lerma2022')
os.makedirs(fig_path, exist_ok=True)
fig.savefig(os.path.join(fig_path, fig_name), bbox_inches='tight')
