import numpy as np
import pandas as pd
import scanpy as sc
import liana as li
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-s','--sample_path', required=True)
parser.add_argument('-t','--thr_gexprop', required=True)
parser.add_argument('-b','--bandwidth', required=True)
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

sample_path = args['sample_path']
thr_gexprop = float(args['thr_gexprop'])
bandwidth = int(args['bandwidth'])
out_path = args['out_path']

# Read slide
slide = sc.read_h5ad(sample_path)

# Run local bivariate scores
li.ut.spatial_neighbors(slide, bandwidth=bandwidth, cutoff=0.1, kernel='gaussian', set_diag=True)
li.mt.lr_bivar(
    slide,
    function_name='cosine', # Name of the function
    n_perms=None, # Number of permutations to calculate a p-value
    mask_negatives=False, # Whether to mask LowLow/NegativeNegative interactions
    add_categories=False, # Whether to add local categories to the results
    expr_prop=thr_gexprop, # Minimum expr. proportion for ligands/receptors and their subunits
    use_raw=False,
    verbose=False
)

# Save
slide_lr = slide.obsm['local_scores'].to_df()
slide_lr.to_csv(out_path)

