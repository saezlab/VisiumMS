import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import os
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import seaborn as sns
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-s','--sn_inp_path', required=True)
parser.add_argument('-v','--vs_inp_path', required=True)
parser.add_argument('-m','--meta_path', required=True)
parser.add_argument('-p','--plot_path', required=True)
args = vars(parser.parse_args())

sn_inp_path = args['sn_inp_path']
vs_inp_path = args['vs_inp_path']
meta_path = args['meta_path']
plot_path = args['plot_path']

# Read meta
print('Read meta')
meta = pd.read_csv(meta_path)

# SN
# Define data and cmap
print('Read adata')
adata = sc.read_h5ad(sn_inp_path)
df = adata.obs.groupby('Sample id').count()[['n_genes']].rename({'n_genes': 'ncells'}, axis=1).reset_index()
df['ncells'] = np.log10(df['ncells'])
df = pd.merge(df, meta, on='Sample id')
df = df.sort_values(['Lesion type'], ascending=False)
_, idxs = np.unique(df['Lesion type'].values.astype('U'), return_inverse=True)
cmap = plt.get_cmap("tab10")
palette = {k:tuple(v) for k,v in zip(df['Sample id'], cmap(idxs))}

print('Plot sn')
# Plot
fig1, axes = plt.subplots(3, 1, figsize=(8,8), facecolor='white', sharex=True, tight_layout=True, dpi=150)

ax = axes[0]
sns.barplot(data=df, x="Sample id", y="ncells", ax=ax, palette=palette, order=palette.keys())
ax.set_ylabel('Number of cells (log10)')
ax.tick_params(axis='x', rotation=90)
ax.set_xlabel('')

ax = axes[1]
df = adata.obs[['Sample id', 'total_counts']].copy().reset_index()
df['log_total_counts'] = np.log10(df['total_counts'].values)
sns.boxplot(data=df, x="Sample id", y="log_total_counts", ax=ax, fliersize=0, palette=palette, order=palette.keys())
ax.set_ylabel('Number of UMIs (log10)')
ax.tick_params(axis='x', rotation=90)
ax.set_ylim(0, None)
ax.set_xlabel('')

ax = axes[2]
df = adata.obs[['Sample id', 'n_genes_by_counts']].copy().reset_index()
df['log_n_genes_by_counts'] = np.log10(df['n_genes_by_counts'].values)
sns.boxplot(data=df, x="Sample id", y="log_n_genes_by_counts", ax=ax, fliersize=0, palette=palette, order=palette.keys())
ax.set_ylabel('Number of genes (log10)')
ax.tick_params(axis='x', rotation=90)
ax.set_ylim(0, None)
ax.set_xlabel('')
fig1.suptitle('QC of sn samples')


# VS
# Read slides
del adata
vs_samples = meta[~meta['Batch vs'].isnull()]['Sample id'].values.astype('U')
adata = []
for sample in vs_samples:
    tmp = sc.read_h5ad(os.path.join(vs_inp_path, sample, 'adata.h5ad'))
    adata.append(tmp)
    del tmp
adata = ad.concat(adata, join='outer')

# Define data
df = adata.obs.groupby('Sample id').count()[['n_genes']].rename({'n_genes': 'ncells'}, axis=1).reset_index()
df['ncells'] = np.log10(df['ncells'])
df = pd.merge(df, meta, on='Sample id')
df = df.sort_values(['Lesion type'], ascending=False)
_, idxs = np.unique(df['Lesion type'].values.astype('U'), return_inverse=True)
cmap = plt.get_cmap("tab10")
palette = {k:tuple(v) for k,v in zip(df['Sample id'], cmap(idxs))}

print('Plot vs')
# Plot
fig2, axes = plt.subplots(3, 1, figsize=(8,8), facecolor='white', sharex=True, tight_layout=True, dpi=150)

ax = axes[0]
sns.barplot(data=df, x="Sample id", y="ncells", ax=ax, palette=palette, order=palette.keys())
ax.set_ylabel('Number of cells (log10)')
ax.tick_params(axis='x', rotation=90)
ax.set_xlabel('')

ax = axes[1]
df = adata.obs[['Sample id', 'total_counts']].copy().reset_index()
df['log_total_counts'] = np.log10(df['total_counts'].values)
sns.boxplot(data=df, x="Sample id", y="log_total_counts", ax=ax, fliersize=0, palette=palette, order=palette.keys())
ax.set_ylabel('Number of UMIs (log10)')
ax.tick_params(axis='x', rotation=90)
ax.set_ylim(0, None)
ax.set_xlabel('')

ax = axes[2]
df = adata.obs[['Sample id', 'n_genes_by_counts']].copy().reset_index()
df['log_n_genes_by_counts'] = np.log10(df['n_genes_by_counts'].values)
sns.boxplot(data=df, x="Sample id", y="log_n_genes_by_counts", ax=ax, fliersize=0, palette=palette, order=palette.keys())
ax.set_ylabel('Number of genes (log10)')
ax.tick_params(axis='x', rotation=90)
ax.set_ylim(0, None)
ax.set_xlabel('')
fig2.suptitle('QC of vs samples')

# Save to pdf
pdf = matplotlib.backends.backend_pdf.PdfPages(plot_path)
for fig in [fig1, fig2]:
    pdf.savefig(fig)
pdf.close()
