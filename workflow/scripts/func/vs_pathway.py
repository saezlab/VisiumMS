import scanpy as sc
import pandas as pd
import numpy as np
import decoupler as dc
import liana as li
import os
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--slide_path', required=True)
parser.add_argument('-g', '--gmt_path', required=True)
parser.add_argument('-n', '--db_name', required=True)
parser.add_argument('-o', '--out_path', required=True)
args = vars(parser.parse_args())

slide_path = args['slide_path']
gmt_path = args['gmt_path']
db_name = args['db_name']
out_path = args['out_path']

# Read data
slide = sc.read_h5ad(slide_path)

def spatial_para(adata):
    from scipy.sparse import csr_matrix
    # Smooth by connectivity
    X = adata.X
    if isinstance(X, csr_matrix):
        X = X.A
    conn = adata.obsp['spatial_connectivities']
    X = (conn @ X) / conn.sum(axis=1).A  # Normalize by spatial weights
    
    return X

# Spatially weight
li.ut.spatial_neighbors(slide, bandwidth=150, cutoff=0.1, kernel='gaussian', set_diag=True)
slide.X = spatial_para(slide)

# Read gene set database
if db_name != 'progeny':
    gmt = dc.read_gmt(gmt_path)
    gmt['source'] = [s.split(db_name)[1].replace('_', ' ').lstrip() for s in gmt['source']]
    weight = None
else:
    gmt = pd.read_csv(gmt_path).groupby('source', observed=True).head(1000)
    weight = 'weight'
msk = ~gmt['source'].str.contains('FETAL|INFECTION|SARS', case=False)
gmt = gmt[msk]


# Run enrichment analysis
dc.run_ulm(slide, gmt, weight=weight, use_raw=False)

# Save
slide.obsm['ulm_estimate'].to_csv(out_path)
