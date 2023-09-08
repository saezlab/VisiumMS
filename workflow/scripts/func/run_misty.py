import numpy as np
import pandas as pd
import scanpy as sc
import liana as li
from liana.method import MistyData
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-d','--sn_lr_path', required=True)
parser.add_argument('-w','--sn_pw_path', required=True)
parser.add_argument('-s','--sample_path', required=True)
parser.add_argument('-l','--lr_path', required=True)
parser.add_argument('-p','--prps_path', required=True)
parser.add_argument('-r','--rctm_path', required=True)
parser.add_argument('-a','--padj', required=True)
parser.add_argument('-b','--bandwidth', required=True)
parser.add_argument('-i','--inters_path', required=True)
parser.add_argument('-t','--metrics_path', required=True)
args = vars(parser.parse_args())

sn_lr_path = args['sn_lr_path']
sn_pw_path = args['sn_pw_path']
sample_path = args['sample_path']
lr_path = args['lr_path']
prps_path = args['prps_path']
rctm_path = args['rctm_path']
padj = float(args['padj'])
bandwidth = int(args['bandwidth'])
inters_path = args['inters_path']
metrics_path = args['metrics_path']

# Read func results for sn
sn_lr = pd.read_csv(sn_lr_path)
sn_pw = pd.read_csv(sn_pw_path)

# Read side
slide = sc.read_h5ad(sample_path)
slide.obsm['lr_scores'] = pd.read_csv(lr_path, index_col=0)
slide.obsm['props'] = pd.read_csv(prps_path, index_col=0)
slide.obsm['pathway'] = pd.read_csv(rctm_path, index_col=0)

# Extract adatas
pathway = li.ut.obsm_to_adata(slide, 'pathway')
props = li.ut.obsm_to_adata(slide, 'props')
lr_scores = li.ut.obsm_to_adata(slide, 'lr_scores')

# Add spatial context
li.ut.spatial_neighbors(pathway, bandwidth=bandwidth, cutoff=0.1, kernel='gaussian', set_diag=True)
li.ut.spatial_neighbors(props, bandwidth=bandwidth, cutoff=0.1, kernel='gaussian', set_diag=True)

# Subset by significant lr and pathways
sn_lr = sn_lr[sn_lr['interaction_padj'] < padj]['interaction'].unique().astype(str)
sn_pw = sn_pw[sn_pw['adj_pvals'] < padj]['source'].unique().astype(str)
sn_lr = sn_lr[np.isin(sn_lr, lr_scores.var_names)]
sn_pw = sn_pw[np.isin(sn_pw, pathway.var_names)]
lr_scores = lr_scores[:, sn_lr].copy()
pathway = pathway[:, sn_pw].copy()

# Build model
misty = MistyData(
    data={"intra": lr_scores, "props": props, "pathway": pathway},
)

# Run
misty(
    model='linear',
    verbose=True,
    bypass_intra=True,
)

# Extract and save
misty.uns['interactions'].to_csv(inters_path, index=False)
misty.uns['target_metrics'].to_csv(metrics_path, index=False)
