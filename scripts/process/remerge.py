import scanpy as sc

import numpy as np
import pandas as pd

import os
import argparse

'''
Delete list of clusters and remerge.
'''

# Read command line and set args
parser = argparse.ArgumentParser(prog='qc', description='Run QC per sample')
parser.add_argument('-i', '--input_dir', help='Input path to annotated object', required=True)
parser.add_argument('-c', '--clusters', help='Clusters to delete', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory', required=True)
args = vars(parser.parse_args())

input_path = args['input_dir']
clusters = args['clusters'].split(',')
output_path = args['output_dir']
###############################

# Read merged object
adata = sc.read_h5ad(input_path).raw.to_adata()
adata.uns['log1p'] = {'base': None} # Fix anndata bug https://github.com/scverse/scanpy/issues/2181

# Delete clusters
for clust in clusters:
    adata = adata[adata.obs['leiden'] != clust]

# Compute HVG
sc.pp.highly_variable_genes(adata, batch_key='batch')
adata.var = adata.var[['highly_variable', 'highly_variable_nbatches']]
adata.raw = adata

# Filter by HVG
num_hvg_genes = 3000
batch_msk = np.array(adata.var.highly_variable_nbatches > 1)
hvg = adata.var[batch_msk].sort_values('highly_variable_nbatches').tail(num_hvg_genes).index
adata.var['highly_variable'] = [g in hvg for g in adata.var.index]
adata.var = adata.var[['highly_variable','highly_variable_nbatches']]
adata = adata[:,hvg]

# Run PCA
sc.pp.scale(adata)
sc.tl.pca(adata, svd_solver='arpack')

# Run UMAP
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# Write to file
adata.write(os.path.join(output_path, 'annotated.h5ad'))
