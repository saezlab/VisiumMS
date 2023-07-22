import scanpy as sc
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--input_path', required=True)
parser.add_argument('-m','--markers_path', required=True)
parser.add_argument('-r','--resolution', required=True)
parser.add_argument('-a', '--annotation', required=True)
parser.add_argument('-p','--plot_path', required=True)
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

input_path = args['input_path']
markers_path = args['markers_path']
resolution = float(args['resolution'])
annotation = args['annotation']
plot_path = args['plot_path']
out_path = args['out_path']

# Read merged object
adata = sc.read_h5ad(input_path)

# Read markers
markers= pd.read_csv(markers_path)
markers = dict(markers.groupby('cell_type')['gene'].apply(list))

# Cluster cells
sc.tl.leiden(adata, resolution=resolution, key_added='leiden')

# Process annotation
if annotation != 'None':
    annotation = {key: value for key, value in (pair.split(":") for pair in annotation.split(";"))}
    # Remove cells that are not in annot
    msk = np.isin(adata.obs['leiden'].values.astype('U'), np.array(list(annotation.keys()), dtype='U'))
    adata = adata[msk].copy()
    # Rename obs
    adata.obs['leiden'] = [annotation[c] for c in adata.obs['leiden']]

# Cluster props
fig1, axes = plt.subplots(1, 2, figsize=(12, 3), tight_layout=True, dpi=150)
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
fig2, axes = plt.subplots(1, 2, figsize=(12, 5), tight_layout=True, dpi=150, gridspec_kw={'width_ratios': [1, 2]})
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
for fig in [fig1, fig2]:
    pdf.savefig(fig)
pdf.close()

# Write
adata.write(out_path)
