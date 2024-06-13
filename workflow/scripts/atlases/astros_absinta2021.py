import pandas as pd
import numpy as np
import decoupler as dc
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import seaborn as sns
from scipy import stats
import scipy
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--absinta_path', required=True)
parser.add_argument('-l','--lerma_path', required=True)
parser.add_argument('-n','--ann_lerma_path', required=True)
parser.add_argument('-b','--deg_path', required=True)
parser.add_argument('-c','--plot_path', required=True)
args = vars(parser.parse_args())

absinta_path = args['absinta_path']
lerma_path = args['lerma_path']
ann_lerma_path = args['ann_lerma_path']
deg_path = args['deg_path']
plot_path = args['plot_path']

# Read
absinta = sc.read_h5ad(absinta_path)

# Filter for AS and remove empry genes
absinta.obs['cluster'] = ['Cluster_' + str(i) for i in absinta.obs['Seurat cluster']]
absinta = absinta[absinta.obs['leiden'] == 'AS'].copy()
absinta.obs['cs'] = pd.Categorical([c if c == 'Cluster_12' else 'rest'  for c in absinta.obs['cluster']], categories=['rest', 'Cluster_12'], ordered=True)
sc.pp.filter_genes(absinta, min_cells=3)

# Compute deg across clusters
sc.tl.rank_genes_groups(absinta, groupby='cluster', method='t-test_overestim_var')
deg_abs = sc.get.rank_genes_groups_df(absinta, group=None)
deg_abs = deg_abs[(deg_abs['pvals_adj'] < 0.05) & (deg_abs['logfoldchanges'] > 0.5)]
deg_ler = pd.read_csv(deg_path)

# Compute jaccard between markers
groups_a = deg_abs['group'].unique()
groups_b = deg_ler['group'].unique()
df = []
for group_a in groups_a:
    for group_b in groups_b:
        genes_a = deg_abs[deg_abs['group'] == group_a].set_index('names').sort_values('scores', ascending=False).head(100)
        genes_b = deg_ler[deg_ler['group'] == group_b].set_index('names').sort_values('scores', ascending=False).head(100)
        inter = np.intersect1d(genes_a.index, genes_b.index)
        diff_ab = np.setdiff1d(genes_a.index, genes_b.index)
        diff_ba = np.setdiff1d(genes_b.index, genes_a.index)
        union = np.union1d(genes_a.index, genes_b.index)
        back = np.union1d(deg_abs['names'].unique().astype('U'), deg_ler['names'].unique().astype('U'))
        back = np.setdiff1d(back, np.union1d(genes_a.index, genes_b.index).astype('U'))
        jacc = inter.size / union.size
        odds, pval = scipy.stats.fisher_exact([[inter.size, diff_ab.size], [diff_ba.size, back.size]], alternative='greater')
        df.append([group_a, group_b, inter.size, union.size, jacc, odds, pval])
df = pd.DataFrame(df, columns=['group_a', 'group_b', 'inter', 'union', 'jacc', 'odds', 'pval'])
df['padj'] = scipy.stats.false_discovery_control(df['pval'], method='bh')
jacc = df.pivot(index='group_a', columns='group_b', values='jacc')
order = jacc.max(0).sort_values(ascending=False).index
jacc = jacc.loc[:, order]

sign = pd.DataFrame(index=jacc.index, columns=jacc.columns)
sign.loc[:, :] = ''

for row in df[df['padj'] < 0.01][['group_a', 'group_b']].iterrows():
    i, row = row
    group_a, group_b = row['group_a'], row['group_b']
    sign.loc[group_a, group_b] = '*'

sign = sign.loc[:, order]

# Plot
fig1, ax = plt.subplots(1, 1, figsize=(4, 4), dpi=150)
htm = sns.heatmap(jacc, cmap='Blues', square=True, cbar_kws={"shrink": 0.25, 'aspect': 5, 'label': 'Jaccard index'}, ax=ax, annot=sign, fmt='', vmax=1, annot_kws={'fontweight': 'black', 'color': 'black'})
ax.set_xlabel('Lerma et al.')
ax.set_ylabel('Absinta et al.')
i = 0
for _, spine in htm.spines.items():
    if i % 2 == 0:
        spine.set_visible(True)
    i += 1

def get_expr_df(adata, gene):
    from scipy.sparse import issparse
    if issparse(adata.X):
        gex = adata[:, gene].copy().X.A.ravel()
    else:
        gex = adata[:, gene].copy().X.ravel()
    df = pd.DataFrame()
    df['group'] = adata.obs['cs']
    df['sample'] = adata.obs_names
    df['value'] = gex
    return df.reset_index(drop=True)


def plot_violin(adata, gene, ax, palette):
    df = get_expr_df(adata, gene)
    sns.violinplot(data=df, x='group', y='value', hue='group', ax=ax, palette=palette)
    ax.set_ylabel(gene)
    ax.set_xlabel('')
    #ax.get_legend().remove()

# Plot top marker genes
fig2, axes = plt.subplots(2, 3, figsize=(5, 3), tight_layout=True, dpi=150)
axes = axes.ravel()
glist = ['DNAH11', 'CFAP299', 'ZBBX', 'SPAG17', 'CFAP54', 'DNAH6']
palette={'rest': 'tab:gray', 'Cluster_12': 'tab:green'}
for i in range(len(glist)):
    plot_violin(adata=absinta, gene=glist[i], ax=axes[i], palette=palette)

# Compute correlation at psbulk level
del absinta.X
absinta.X = absinta.layers['counts'].copy()
del absinta.layers['counts']
absinta = dc.get_pseudobulk(
    adata=absinta,
    sample_col='cluster',
    groups_col=None,
    mode='sum',
    min_cells=0,
    min_counts=0
)
g_expr = dc.filter_by_expr(absinta, group='cluster', min_count=10, min_total_count=15)
g_prop = dc.filter_by_prop(absinta, min_prop=0.2, min_smpls=2)
genes = np.intersect1d(g_expr, g_prop)
absinta = absinta[:, genes]

lerma = sc.read_h5ad(lerma_path)
lerma.obs['cell_states'] = pd.read_csv(ann_lerma_path, index_col=0)
del lerma.X
lerma.X = lerma.layers['counts']
del lerma.layers['counts']
lerma = dc.get_pseudobulk(
    adata=lerma,
    sample_col='cell_states',
    groups_col=None,
    mode='sum',
    min_cells=0,
    min_counts=0
)
g_expr = dc.filter_by_expr(lerma, group='cell_states', min_count=10, min_total_count=15)
g_prop = dc.filter_by_prop(lerma, min_prop=0.2, min_smpls=2)
genes = np.intersect1d(g_expr, g_prop)
lerma = lerma[:, genes]
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
lerma_order = ['AS_C', 'AS_GM', 'AS_R', 'AS_Dis1', 'AS_TC', 'AS_Homeo', 'AS_Dis2', 'AS_Phago', 'AS_Dis3', 'AS_stress']
absinta_order = ['Cluster_12', 'Cluster_14', 'Cluster_4', 'Cluster_9']
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
ann = pd.DataFrame(index=pvals.index, columns=pvals.columns)
ann.loc[:, :] = np.where((pvals.values < 0.05) & (np.abs(corrs.values) > 0.75), '*', '')

# Plot
cmap = plt.get_cmap('Blues').copy()
cmap.set_bad(color='gray')

fig_corrs, ax = plt.subplots(1, 1, figsize=(4, 4), facecolor='white', dpi=150)
htm = sns.heatmap(corrs, cmap=cmap, square=True, vmax=1, vmin=0, ax=ax, cbar_kws={"shrink": 0.25, 'aspect': 5, 'label': 'Pearson correlation'},
                  annot=ann, fmt='', annot_kws={'fontweight': 'black', 'color': 'black'})
i = 0
for _, spine in htm.spines.items():
    if i % 2 == 0:
        spine.set_visible(True)
    i += 1

ax.set_xlabel('Lerma et al.')
ax.set_ylabel('Absinta et al.')
ax.tick_params(axis='x', labelrotation=90)

# Save to pdf
pdf = matplotlib.backends.backend_pdf.PdfPages(plot_path)
for fg in [fig1, fig_corrs, fig2]:
    pdf.savefig(fg, bbox_inches='tight')
pdf.close()
