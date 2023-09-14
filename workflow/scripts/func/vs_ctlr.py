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

# Generate feature universe
sn_ctlr = (
    sn_lr
    .rename(columns={'receptor_complex': 'receptor', 'ligand_complex': 'ligand'})
    [['source', 'target', 'ligand', 'receptor']]
    .assign(ligand=lambda x: x['ligand'].str.split('_').str[0], receptor=lambda x: x['receptor'].str.split('_').str[0])
    .drop_duplicates()
)
it_universe = [x[1:] for x in sn_ctlr.itertuples()]

# Read slide
slide = sc.read_h5ad(sample_path)
slide.X = (slide.X > 0) * 1
props = AnnData(pd.read_csv(props_path, index_col=0), dtype=float)
props.X = (props.X > 0.11) * 1
props.obsm['spatial'] = slide.obsm['spatial'].copy()
props.uns['spatial'] = slide.uns['spatial'].copy()

# Merge features
merged = ad.concat([slide, props], axis=1)
merged.obsm['spatial'] = slide.obsm['spatial'].copy()
merged.uns['spatial'] = slide.uns['spatial'].copy()

# Compute local spatial product score
def spatial_prod(adata, interactions, sep='^'):
    # Smooth by connectivity
    X = adata.X
    if isinstance(X, csr_matrix):
        X = X.A
    conn = adata.obsp['spatial_connectivities']
    X = (conn @ X) / conn.sum(axis=1).A  # Normalize by spatial weights
    X = X / np.max(X, axis=0)  # Make features comparable
    new_X = np.zeros((X.shape[0], len(interactions)))
    vars = list(adata.var_names)
    cols = []
    for i, interaction in enumerate(tqdm(interactions)):
        if np.all(np.isin(interaction, vars)):
            idxs = np.unique([vars.index(a) for a in interaction])  # Make unique
            new_X[:, i] = np.prod(X[:, idxs], axis=1)
        cols.append('{0}'.format(sep).join(interaction))

    new_adata = AnnData(pd.DataFrame(
        new_X,
        columns=cols,
        index=adata.obs_names
    ), dtype=float)
    new_adata.obsm['spatial'] = slide.obsm['spatial'].copy()
    new_adata.uns['spatial'] = slide.uns['spatial'].copy()
    return new_adata

# Compute cell type + ligand-receptor scores
li.ut.spatial_neighbors(merged, bandwidth=bandwidth, cutoff=0.1, kernel='gaussian', set_diag=True)
ctlr_res = spatial_prod(merged, it_universe, sep='^')

# Save to csv
ctlr_res.to_df().to_csv(out_path)
