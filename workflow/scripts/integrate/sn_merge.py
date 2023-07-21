import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--input_path', required=True)
parser.add_argument('-m','--meta_path', required=True)
parser.add_argument('-n','--n_hvg', required=True)
parser.add_argument('-p','--plot_path', required=True)
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

input_path = args['input_path']
meta_path = args['meta_path']
n_hvg = int(args['n_hvg'])
plot_path = args['plot_path']
out_path = args['out_path']

# Load meta data
meta = pd.read_csv(meta_path)
sn_samples = meta[~meta['Batch sn'].isnull()]['Sample id'].values.astype('U')

# Iter each sample
adata = []
for sample in sn_samples:

    # Read adata
    tmp = sc.read_h5ad(os.path.join(input_path, sample, 'adata.h5ad'))
    tmp.obs_names = [sample + '-' + b.split('-')[0] for b in tmp.obs_names]
    print(sample, tmp.shape)

    # Fetch sample metadata
    m = meta[meta['Sample id'] == sample]

    # Add metadata to adata
    for col in m.columns:
        tmp.obs[col] = m[col].values[0]

    # Append
    adata.append(tmp)
    del tmp

# Merge objects and delete list
adata = ad.concat(adata, join='outer')

# Store raw counts
adata.layers['counts'] = adata.X.copy()

# Log-normalize expression
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.layers['lognorm'] = adata.X.copy()

# Compute HVG
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
fig, ax = plt.subplots(2, 2, figsize=(8, 6), dpi=150, tight_layout=True)
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

fig.savefig(plot_path, bbox_inches='tight')

# Save
adata.write(out_path)

