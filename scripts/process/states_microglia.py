import scanpy as sc
import pandas as pd
import numpy as np
import os
import scanpy.external as sce


# Read adata
adata = sc.read('data/prc/sc/raw.h5ad')

# Filtr Microglia
msk = np.array(['Micro' in c for c in adata.obs['cell_type'].values])
adata = adata[msk]

# Process
adata.var_names_make_unique()
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_genes=200)

# Normalize
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Filter by HVG
sc.pp.highly_variable_genes(adata, batch_key='batch')
adata.var = adata.var[['highly_variable','highly_variable_nbatches']]

num_hvg_genes = 3000
batch_msk = np.array(adata.var.highly_variable_nbatches > 1)
hvg = adata.var[batch_msk].sort_values('highly_variable_nbatches').tail(num_hvg_genes).index
adata.var['highly_variable'] = [g in hvg for g in adata.var.index]
adata.var = adata.var[['highly_variable','highly_variable_nbatches']]
adata.raw = adata
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack', use_highly_variable=True)

# Integrate
sce.pp.harmony_integrate(adata, 'batch', basis='X_pca', max_iter_harmony=30)

# Generate UMAP features
sc.pp.neighbors(adata, use_rep='X_pca_harmony')
sc.tl.umap(adata)

# Run leiden clustering algorithm
sc.tl.leiden(adata, resolution=0.5)

# Write
adata.write('data/prc/sc/microglia.h5ad')
