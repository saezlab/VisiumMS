import argparse
import os
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-d','--adata_path', required=True)
parser.add_argument('-c','--cstates_path', required=True)
parser.add_argument('-p','--plot_path', required=True)
parser.add_argument('-a','--ann_path', required=True)
parser.add_argument('-g','--deg_path', required=True)
args = vars(parser.parse_args())

adata_path = args['adata_path']
cstates_path = args['cstates_path']
plot_path = args['plot_path']
ann_path = args['ann_path']
deg_path = args['deg_path']

# Read df
df_ct = pd.read_csv(cstates_path)

# Find cell type from path and subset
cell_type = os.path.basename(adata_path).split('.h5ad')[0]
msk = df_ct['name'].str.startswith(cell_type)
df_ct = df_ct[msk]

# Get names_dict
names_dict = dict()
for name, leidens in zip(df_ct['name'], df_ct['leidens']):
    for leiden in leidens.split(';'):
        names_dict[leiden] = name

# Get markers list
markers = []
for r in df_ct['markers'].dropna().str.split(';'):
    markers.extend(r)

# Get custom order
custom_order = df_ct['name'].values

# Read adata
adata = sc.read_h5ad(adata_path)

# Update metadata
adata.obs['cell_states'] = [names_dict[num] for num in adata.obs['leiden']]
adata.obs['cell_states'] = pd.Categorical(
    adata.obs['cell_states'],
    categories=custom_order,
    ordered=True
)

# Compute DEG
sc.tl.rank_genes_groups(adata, groupby='cell_states', method='t-test_overestim_var')
deg = sc.get.rank_genes_groups_df(adata, group=None)

# Filter deg to keep marker genes
deg = deg[(deg['pvals_adj'] < 0.05) & (deg['logfoldchanges'] > 0.5)]

# Umap and dotplot
fig1, axes = plt.subplots(1, 2, figsize=(12, 5), tight_layout=True, dpi=150, gridspec_kw={'width_ratios': [1, 2]})
axes = axes.ravel()

# UMAP
ax = axes[0]
sc.pl.umap(adata, color='cell_states', return_fig=False, show=False, ax=ax, frameon=False, size=20)
ax.set_title('')

# Dotplot
ax = axes[1]
if 'dendrogram_cell_states' in adata.uns:
    del adata.uns['dendrogram_cell_states']
sc.pl.dotplot(adata, markers, 'cell_states', dendrogram=False, ax=ax, return_fig=False, show=False, dot_max=0.5, standard_scale='var')

# Cluster props
fig2, axes = plt.subplots(1, 2, figsize=(12, 3), tight_layout=True, dpi=150)
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
df = adata.obs.groupby(['cell_states', 'Sample id'])[['total_counts']].count().reset_index()
df = df.pivot(index='cell_states', columns='Sample id', values='total_counts')
stackbar(y=df.values, type_names=df.columns, title='Cluster representation', level_names=df.index, cmap=adata.uns['Sample id_colors'], ax=ax)

ax = axes[1]
df = adata.obs.groupby(['cell_states', 'Lesion type'])[['total_counts']].count().reset_index()
df = df.pivot(index='cell_states', columns='Lesion type', values='total_counts')
stackbar(y=df.values, type_names=df.columns, title='Cluster representation', level_names=df.index, cmap=adata.uns['Lesion type_colors'], ax=ax)

# Save to pdf
pdf = matplotlib.backends.backend_pdf.PdfPages(plot_path)
for fig in [fig1, fig2]:
    pdf.savefig(fig)
pdf.close()

# Save results
adata.obs[['cell_states']].to_csv(ann_path)
deg.to_csv(deg_path, index=False)
