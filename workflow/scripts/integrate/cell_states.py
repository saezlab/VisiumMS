import argparse
import os
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import decoupler as dc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import scanpy.external as sce

n_hvg = 4096
resolution = 1.

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--inp_path', required=True)
parser.add_argument('-p','--plot_path', required=True)
args = vars(parser.parse_args())

inp_path = args['inp_path']
plot_path = args['plot_path']

ctype = os.path.basename(plot_path).replace('states_', '').replace('.pdf', '')

adata = sc.read_h5ad(inp_path)
del adata.uns['log1p']

adata = adata[adata.obs['leiden'] == ctype]
sc.pp.filter_genes(adata, min_cells=3)

adata.layers['lognorm'] = adata.X.copy()

# Compute HVG
msk = adata.obs.groupby('Sample id').count()['total_counts'] >= 5
keep = msk[msk].index.values.astype('U')
adata = adata[np.isin(adata.obs['Sample id'], keep), :].copy()

sc.pp.highly_variable_genes(adata, batch_key='Sample id')
batch_msk = np.array(adata.var['highly_variable_nbatches'] > 1)
hvg = adata.var[batch_msk].sort_values(['highly_variable_nbatches', 'dispersions_norm'], ascending=[False, False]).head(n_hvg).index
adata.var['highly_variable'] = np.isin(adata.var.index, hvg)
del adata.uns['hvg']

# PCA
sc.pp.scale(adata)
sc.tl.pca(adata, svd_solver='arpack', use_highly_variable=True)
pca_var = adata.uns['pca']['variance_ratio']
del adata.uns['pca']
del adata.varm
adata.X = adata.layers['lognorm'].copy()
del adata.layers['lognorm']
adata.var = adata.var[['highly_variable','highly_variable_nbatches']]

# Plots
fig1, ax = plt.subplots(2, 2, figsize=(8, 6), dpi=150, tight_layout=True)
ax = ax.ravel()

# Plot var explained
x = np.arange(1, pca_var.size + 1)
ax[0].scatter(x, pca_var)
ax[0].grid(True)
ax[0].set_axisbelow(True)
ax[0].set_ylabel('Variance explained')
ax[0].set_xlabel('PCs')

# Plot pca
sc.pl.pca(adata, color='Sample id', ax=ax[1], return_fig=False, show=False, frameon=True)
ax[1].set_title('')

# Plot HVG
_ = ax[2].hist(adata.var['highly_variable_nbatches'].values, cumulative=-1)
_ = ax[2].hist(adata.var[adata.var['highly_variable']]['highly_variable_nbatches'].values, cumulative=-1)
ax[2].axhline(y=n_hvg, linestyle='--', c='black')
ax[2].grid(True)
ax[2].set_axisbelow(True)
ax[2].set_xlabel('Number of samples')
ax[2].set_ylabel('Number of HVG')

# Plot cummulative cells
sample_df = adata.obs[['Sample id', 'total_counts']].groupby(['Sample id']).count().reset_index()
bottom = 0
i = 0
for label, height in sample_df.values:
    ax[3].bar(x='Sample id', height=height, bottom=bottom, label=label, color=adata.uns['Sample id_colors'][i])
    bottom += height
    i += 1
ax[3].grid(False)
ax[3].set_axisbelow(True)
ax[3].margins(0)
ax[3].set_ylabel('Number of cells')

# Run harmony
sce.pp.harmony_integrate(adata, key='Sample id', adjusted_basis='X_pca', max_iter_harmony=30)

# Run umap with updated connectivity
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# Plot
fig2, axes = plt.subplots(2, 3, figsize=(12, 6), tight_layout=True, dpi=150)
axes = axes.ravel()

ax = axes[0]
sc.pl.umap(adata, color='doublet_score', ax=ax, return_fig=False, show=False, vmin=0, vmax=0.5)
ax.set_title('Doublet score')

ax = axes[1]
sc.pl.umap(adata, color='pct_counts_mt', ax=ax, return_fig=False, show=False, vmin=0, vmax=10)
ax.set_title('Fraction MT counts')

ax = axes[2]
sc.pl.umap(adata, color='Sample id', ax=ax, return_fig=False, show=False)

ax = axes[3]
adata.obs['log_n_genes_by_counts'] = np.log10(adata.obs['n_genes_by_counts'])
sc.pl.umap(adata, color='log_n_genes_by_counts', ax=ax, return_fig=False, show=False, vmin=0)
ax.set_title('Number of genes expr (log10)')
del adata.obs['log_n_genes_by_counts']

ax = axes[4]
adata.obs['log_total_counts'] = np.log10(adata.obs['total_counts'])
sc.pl.umap(adata, color='log_total_counts', ax=ax, return_fig=False, show=False, vmin=0)
ax.set_title('Total counts (log10)')
del adata.obs['log_total_counts']

ax = axes[5]
sc.pl.umap(adata, color='Lesion type', ax=ax, return_fig=False, show=False)

# Cluster cells
sc.tl.leiden(adata, resolution=resolution, key_added='leiden')

# Find markers
markers = sc.tl.rank_genes_groups(adata, 'leiden')
markers = sc.get.rank_genes_groups_df(adata, group=None)
msk = ~markers['names'].str.contains('MT-|LINC|MIR', case=False)
markers = markers.loc[msk].groupby('group').head(5)
markers = markers.groupby('group')['names'].apply(lambda x: np.array(x, dtype='U')).to_dict()


# Cluster props
fig3, axes = plt.subplots(1, 2, figsize=(12, 3), tight_layout=True, dpi=150)
axes = axes.ravel()

def stackbar(y, type_names, title, level_names, cmap, ax):

    n_bars, n_types = y.shape

    r = np.array(range(n_bars))
    sample_sums = np.sum(y, axis=1)

    barwidth = 0.85
    cum_bars = np.zeros(n_bars)

    for n in range(n_types):
        bars = [i / j * 100 for i, j in zip([y[k][n] for k in range(n_bars)], sample_sums)]
        ax.bar(r, bars, bottom=cum_bars, color=cmap[n], width=barwidth, label=type_names[n], linewidth=0)
        cum_bars += bars

    ax.set_title(title)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, ncol=2)
    ax.set_xticks(r)
    ax.set_xticklabels(level_names, rotation=90)
    ax.set_ylabel("Proportion")
    ax.grid(False)
    ax.margins(0)

    return ax

# Barplot
ax = axes[0]
df = adata.obs.groupby(['leiden', 'Sample id'])[['total_counts']].count().reset_index()
df = df.pivot(index='leiden', columns='Sample id', values='total_counts')
stackbar(y=df.values, type_names=df.columns, title='Cluster representation', level_names=df.index, cmap=adata.uns['Sample id_colors'], ax=ax)

ax = axes[1]
df = adata.obs.groupby(['leiden', 'Lesion type'])[['total_counts']].count().reset_index()
df = df.pivot(index='leiden', columns='Lesion type', values='total_counts')
stackbar(y=df.values, type_names=df.columns, title='Cluster representation', level_names=df.index, cmap=adata.uns['Lesion type_colors'], ax=ax)

# Umap and dotplot
fig4, axes = plt.subplots(1, 2, figsize=(12, 5), tight_layout=True, dpi=150, gridspec_kw={'width_ratios': [1, 4]})
axes = axes.ravel()

# UMAP
ax = axes[0]
sc.pl.umap(adata, color='leiden', return_fig=False, show=False, ax=ax, frameon=False, size=5)
ax.set_title('')

# Dotplot
ax = axes[1]
if 'dendrogram_leiden' in adata.uns:
    del adata.uns['dendrogram_leiden']
sc.pl.dotplot(adata, markers, 'leiden', dendrogram=True, ax=ax, return_fig=False, show=False)

# Save to pdf
pdf = matplotlib.backends.backend_pdf.PdfPages(plot_path)
for fig in [fig1, fig2, fig3, fig4]:
    pdf.savefig(fig)
pdf.close()

