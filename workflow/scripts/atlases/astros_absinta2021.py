import pandas as pd
import numpy as np
import decoupler as dc
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import seaborn as sns
import scipy
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--absinta_path', required=True)
parser.add_argument('-b','--deg_path', required=True)
parser.add_argument('-c','--plot_path', required=True)
args = vars(parser.parse_args())

absinta_path = args['absinta_path']
deg_path = args['deg_path']
plot_path = args['plot_path']

# Read
adata = sc.read_h5ad(absinta_path)

# Filter for AS and remove empry genes
adata.obs['cluster'] = ['Cluster_' + str(i) for i in adata.obs['Seurat cluster']]
adata = adata[adata.obs['leiden'] == 'AS'].copy()
adata.obs['cs'] = pd.Categorical([c if c == 'Cluster_12' else 'rest'  for c in adata.obs['cluster']], categories=['rest', 'Cluster_12'], ordered=True)
sc.pp.filter_genes(adata, min_cells=3)

# Compute deg across clusters
sc.tl.rank_genes_groups(adata, groupby='cluster', method='t-test_overestim_var')
deg_abs = sc.get.rank_genes_groups_df(adata, group=None)
deg_abs = deg_abs[(deg_abs['pvals_adj'] < 0.05) & (deg_abs['logfoldchanges'] > 0.5)]
deg_ler = pd.read_csv(deg_path)

# Compute jaccard between markers
groups_a = deg_abs['group'].unique()
groups_b = deg_ler['group'].unique()
df = []
for group_a in groups_a:
    for group_b in groups_b:
        genes_a = deg_abs[deg_abs['group'] == group_a].set_index('names').head(100)
        genes_b = deg_ler[deg_ler['group'] == group_b].set_index('names').head(100)
        inter = np.intersect1d(genes_a.index, genes_b.index)
        diff_ab = np.setdiff1d(genes_a.index, genes_b.index)
        diff_ba = np.setdiff1d(genes_b.index, genes_a.index)
        union = np.union1d(genes_a.index, genes_b.index)
        jacc = inter.size / union.size
        odds, pval = scipy.stats.fisher_exact([[inter.size, diff_ab.size], [diff_ba.size, 20000]], alternative='greater')
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

fig1, ax = plt.subplots(1, 1, figsize=(6, 3), dpi=150)
sns.heatmap(jacc, cmap='Blues', square=True, cbar_kws={"shrink": 0.25, 'aspect': 5, 'label': 'Jaccard index'}, ax=ax, annot=sign, fmt='')
ax.set_xlabel('Lerma et al.')
ax.set_ylabel('Absinta et al.')

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
    plot_violin(adata=adata, gene=glist[i], ax=axes[i], palette=palette)

# Save to pdf
pdf = matplotlib.backends.backend_pdf.PdfPages(plot_path)
for fg in [fig1, fig2]:
    pdf.savefig(fg, bbox_inches='tight')
pdf.close()
