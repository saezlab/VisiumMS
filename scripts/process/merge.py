import scanpy as sc
import scanpy.external as sce

import numpy as np
import pandas as pd

import os
import argparse

'''
Open all samples QC processed files, concatenate them and run integration
'''

# Read command line and set args
parser = argparse.ArgumentParser(prog='qc', description='Run QC per sample')
parser.add_argument('-i', '--input_dir', help='Input directory containing all sample directories', required=True)
parser.add_argument('-m', '--metadata', help='Path to metadata', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the processed object', required=True)
args = vars(parser.parse_args())

input_path = args['input_dir']
meta_path = args['metadata']
output_path = args['output_dir']
###############################

# Load meta data
meta = pd.read_csv(meta_path)
samples = np.unique(meta['sample_id'])

adata = []
for sample in os.listdir(input_path):
    if not os.path.isdir(os.path.join(input_path, sample)) or sample.startswith('.'):
        continue
    print(sample)
    
    # Read adata
    tmp = sc.read_h5ad(os.path.join(input_path, sample, sample+'.h5ad'))
    
    # Fetch sample metadata
    m = meta[meta['sample_id'] == sample]
    
    # Add metadata to adata
    for col in m.columns:
        tmp.obs[col] = m[col].values[0]
    
    # Append
    adata.append(tmp)
    del tmp
    
# Merge objects and delete list
adata = adata[0].concatenate(adata[1:], join='outer')

# Log-normalize expression
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# Compute HVG
sc.pp.highly_variable_genes(adata, batch_key='batch')
adata.var = adata.var[['highly_variable','highly_variable_nbatches']]
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
adata.write(os.path.join(output_path, 'merged.h5ad'))
