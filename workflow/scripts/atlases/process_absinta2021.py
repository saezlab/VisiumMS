import pandas as pd
import numpy as np
from anndata import AnnData
import scanpy as sc
import os
import argparse
from scipy.sparse import csr_matrix

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--ann_path', required=True)
parser.add_argument('-c','--cnt_path', required=True)
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

ann_path = args['ann_path']
cnt_path = args['cnt_path']
out_path = args['out_path']

# Create obs
obs = pd.read_csv(ann_path, sep='\t', index_col=0)
obs.index.name = 'obs_names'
obs = obs.rename(columns={'NBB_case': 'Sample id', 'pathology': 'Lesion type', 'seurat_cluster': 'Seurat cluster', 'cell_type': 'leiden'})

# Rename lesion type
lt_dict = {
    'chronic_active_MS_lesion_edge': 'CA',
    'chronic_inactive_MS_lesion_edge': 'CI',
    'MS_lesion_core': 'CA & CI',
    'MS_periplaque_white_matter': 'CA & CI',
    'control_white_matter': 'Ctrl'
}
obs['Lesion type'] = [lt_dict[l] for l in obs['Lesion type']]

# Rename leiden
ct_dict = {
    'oligodendrocytes': 'OL',
    'opc': 'OPC',
    'astrocytes': 'AS',
    'immune': 'MG',
    'neurons': 'NEU',
    'vascular_cells': 'EC',
    'lymphocytes': 'BC & TC'
}
obs['leiden'] = [ct_dict[l] for l in obs['leiden']]
print(obs.shape)

# Create adata
adata = AnnData(pd.read_csv(cnt_path).T, dtype=float, obs=obs)
adata.X = csr_matrix(adata.X)

# QC
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Normalize
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# Write
print(adata.shape)
adata.write(out_path)
