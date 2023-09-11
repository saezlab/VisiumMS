import numpy as np
import pandas as pd
import scanpy as sc
import liana as li
import anndata as ad
from anndata import AnnData
from scipy.sparse import csr_matrix
from tqdm import tqdm
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-s','--sample_path', required=True)
parser.add_argument('-p','--props_path', required=True)
parser.add_argument('-n','--sn_lr_path', required=True)
parser.add_argument('-b','--bandwidth', required=True)
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

sample_path = args['sample_path']
props_path = args['props_path']
sn_lr_path = args['sn_lr_path']
bandwidth = int(args['bandwidth'])
out_path = args['out_path']

# Read func results for sn
sn_lr = pd.read_csv(sn_lr_path)

# Generate feature universes
lr_universe = (
    sn_lr
    .rename(columns={'receptor_complex': 'receptor', 'ligand_complex': 'ligand'})
    [['ligand', 'receptor']]
    .assign(ligand=lambda x: x['ligand'].str.split('_').str[0], receptor=lambda x: x['receptor'].str.split('_').str[0])
    .drop_duplicates()
)
sn_ctlr = (
    sn_lr
    .rename(columns={'receptor_complex': 'receptor', 'ligand_complex': 'ligand'})
    [['source', 'target', 'ligand', 'receptor']]
    .assign(ligand=lambda x: x['ligand'].str.split('_').str[0], receptor=lambda x: x['receptor'].str.split('_').str[0])
    .drop_duplicates()
)
lr_universe = [x[1:] for x in lr_universe.itertuples()]
ct_universe = [x[1:] for x in sn_lr[['source', 'target']].drop_duplicates().itertuples()]

# Read slide
slide = sc.read_h5ad(sample_path)
props = AnnData(pd.read_csv(props_path, index_col=0), dtype=float)
props.obsm['spatial'] = slide.obsm['spatial'].copy()
props.uns['spatial'] = slide.uns['spatial'].copy()

# Compute local spatial product score
def spatial_prod(adata, interactions, sep='^'):
    # Smooth by connectivity
    X = adata.X
    if isinstance(X, csr_matrix):
        X = X.A
    conn = adata.obsp['spatial_connectivities']
    X = (conn @ X) / conn.sum(axis=1).A  # Normalize by spatial weights
    new_X = np.zeros((X.shape[0], len(interactions)))
    vars = list(adata.var_names)
    cols = []
    for i, interaction in enumerate(tqdm(interactions)):
        a, b = interaction
        if (a in vars) and (b in vars):
            idx_a = vars.index(a)
            idx_b = vars.index(b)
            new_X[:, i] = X[:, idx_a] * X[:, idx_b]  # Product
        cols.append('{0}{1}{2}'.format(a, sep, b))

    new_adata = AnnData(pd.DataFrame(
        new_X,
        columns=cols,
        index=adata.obs_names
    ), dtype=float)
    new_adata.obsm['spatial'] = slide.obsm['spatial'].copy()
    new_adata.uns['spatial'] = slide.uns['spatial'].copy()
    return new_adata

# Compute cell type scores
li.ut.spatial_neighbors(props, bandwidth=bandwidth, cutoff=0.1, kernel='gaussian', set_diag=True)
ct_res = spatial_prod(props, ct_universe)

# Compute ligand receptor scores
li.ut.spatial_neighbors(slide, bandwidth=bandwidth, cutoff=0.1, kernel='gaussian', set_diag=True)
lr_res = spatial_prod(slide, lr_universe)

# Merge results
ctlr_slide = ad.concat([ct_res, lr_res], axis=1)
ctlr_slide.obsm['spatial'] = slide.obsm['spatial'].copy()
ctlr_slide.uns['spatial'] = slide.uns['spatial'].copy()
it_universe = [('^'.join(x[1:3]), '^'.join(x[3:])) for x in sn_ctlr.itertuples()]

# Compute cell type + ligand-receptor scores
li.ut.spatial_neighbors(ctlr_slide, bandwidth=bandwidth, cutoff=0.1, kernel='gaussian', set_diag=True)
ctlr_res = spatial_prod(ctlr_slide, it_universe, sep='|')

# Save to csv
ctlr_res.to_df().to_csv(out_path)
