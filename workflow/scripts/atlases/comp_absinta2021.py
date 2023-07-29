import decoupler as dc
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
import os
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-c','--lerma_path', required=True)
parser.add_argument('-a','--absinta_path', required=True)
parser.add_argument('-p','--plot_path', required=True)
args = vars(parser.parse_args())

lerma_path = args['lerma_path']
absinta_path = args['absinta_path']
plot_path = args['plot_path']

# Read and process lerma2023
lerma = sc.read_h5ad(lerma_path)
del lerma.X
lerma.X = lerma.layers['counts']
del lerma.layers['counts']
lerma = dc.get_pseudobulk(
    adata=lerma,
    sample_col='leiden',
    groups_col=None,
    mode='sum',
    min_cells=0,
    min_counts=0
)
g_expr = dc.filter_by_expr(lerma, group='leiden', min_count=10, min_total_count=15)
g_prop = dc.filter_by_prop(lerma, min_prop=0.2, min_smpls=2)
genes = np.intersect1d(g_expr, g_prop)
lerma = lerma[:, genes]

# Read and process absinta2021
absinta = sc.read_h5ad(absinta_path)
del absinta.X
absinta.X = absinta.layers['counts'].copy()
del absinta.layers['counts']
absinta = dc.get_pseudobulk(
    adata=absinta,
    sample_col='leiden',
    groups_col=None,
    mode='sum',
    min_cells=0,
    min_counts=0
)
g_expr = dc.filter_by_expr(absinta, group='leiden', min_count=10, min_total_count=15)
g_prop = dc.filter_by_prop(absinta, min_prop=0.2, min_smpls=2)
genes = np.intersect1d(g_expr, g_prop)
absinta = absinta[:, genes]

# Find shared genes
inter = lerma.var_names.intersection(absinta.var_names)
lerma = lerma[:, inter].copy()
absinta = absinta[:, inter].copy()

# Normalize
sc.pp.normalize_total(lerma)
sc.pp.log1p(lerma)
sc.pp.normalize_total(absinta)
sc.pp.log1p(absinta)

# Compute corrs
lerma_order = ['AS', 'BC', 'TC', 'EC', 'SC', 'MG', 'NEU', 'OL', 'OPC']
absinta_order = ['AS', 'BC & TC', 'EC', 'MG', 'NEU', 'OL', 'OPC']
corrs = []
pvals = []
for c_a in absinta_order:
    c_row = []
    p_row = []
    for c_b in lerma_order:
        x = absinta[c_a, :].X.ravel()
        y = lerma[c_b, :].X.ravel()
        c, p = stats.pearsonr(x, y)
        c_row.append(c)
        p_row.append(p)
    corrs.append(c_row)
    pvals.append(p_row)

corrs = pd.DataFrame(corrs, index=absinta_order, columns=lerma_order)
pvals = pd.DataFrame(pvals, index=absinta_order, columns=lerma_order)

# Adjust by fdr
pvals.loc[:, :] = dc.p_adjust_fdr(pvals.values.ravel()).reshape(pvals.shape)

# Find significant
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

ax.set_ylabel('Lerma-Martin (2023)')
ax.set_xlabel('Absinta (2021)')
ax.tick_params(axis='x', labelrotation=90)

fig.savefig(plot_path, bbox_inches='tight')

