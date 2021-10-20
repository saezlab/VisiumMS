import scanpy as sc
import scanpy.external as sce
import numpy as np
import pandas as pd

import argparse
import os

"""
Script to plot different QC metrics after filtering the data.
"""

# Read command line and set args
parser = argparse.ArgumentParser(prog='qc', description='Run QC per sample')
parser.add_argument('-i', '--input_path', help='Input path to merged object', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
args = vars(parser.parse_args())

input_path = args['input_path']
output_path = args['output_dir']
###############################

# Read merged object
adata = sc.read_h5ad(input_path)

# Run harmony
sce.pp.harmony_integrate(adata, 'batch', adjusted_basis='X_pca', max_iter_harmony=30)
sc.pp.neighbors(adata)

# Run umap with updated connectivity
sc.tl.umap(adata)

# Write to file
adata.write(os.path.join(output_path, 'integrated.h5ad'))
