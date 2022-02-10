"""
Script to generate raw data from a processed, annotated object. 
"""

import numpy as np
import pandas as pd

import scanpy as sc

import os

from anndata import AnnData
import argparse


# Read command line and set args
parser = argparse.ArgumentParser(prog='qc', description='Generates raw data object')
parser.add_argument('-s', '--path_samples', help='Input path to raw samples', required=True)
parser.add_argument('-a', '--path_ann_obj', help='Annotated AnnData object', required=True)
parser.add_argument('-o', '--path_output', help='Path were to save output', required=True)
args = vars(parser.parse_args())

path_samples = args['path_samples']
path_ann_obj = args['path_ann_obj']
path_output = args['path_output']
###############################

meta = sc.read(path_ann_obj).obs

adata = []
for dname in os.listdir(path_samples):
    # Read
    tmp = sc.read_10x_mtx(os.path.join(path_samples, dname, 'filtered_feature_bc_matrix'),
                            var_names='gene_symbols', cache=True)
    tmp.var_names_make_unique()
    
    # Store
    adata.append(tmp)
    
adata = adata[0].concatenate(adata[1:], join='outer')
adata = adata[meta.index]
adata.obs = meta

# Write
adata.write(path_output)
