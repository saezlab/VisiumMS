import scanpy as sc
import numpy as np

import os
from scvi.model import CondSCVI, DestVI

import argparse


# Read command line and set args
parser = argparse.ArgumentParser(prog='qc', description='Generates regression model from raw single nuc data')
parser.add_argument('-r', '--path_raw', help='Path to raw single nuclei adata', required=True)
parser.add_argument('-s', '--path_slides', help='Path to raw visium slides', required=True)
parser.add_argument('-n', '--label_name', help='Label of cell type from raw data', required=True)
parser.add_argument('-d', '--sample_id', help='Label of sample id from raw data', required=True)
parser.add_argument('-o', '--path_output', help='Path were to save model', required=True)
args = vars(parser.parse_args())

path_raw = args['path_raw']
path_slides = args['path_slides']
label_name = args['label_name']
sample_id = args['sample_id']
path_output = args['path_output']

# Read raw data
adata = sc.read_h5ad(input_path)

# Filter empty genes
sc.pp.filter_genes(adata, min_counts=10)

# Store raw expression in layer
adata.layers["counts"] = adata.X.copy()

# Filter HVG genes
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, batch_key='batch')
adata.var = adata.var[['highly_variable','highly_variable_nbatches']]
adata.raw = adata
batch_msk = np.array(adata.var.highly_variable_nbatches > 1)
hvg = adata.var[batch_msk].sort_values('highly_variable_nbatches').tail(num_hvg_genes).index
adata.var['highly_variable'] = [g in hvg for g in adata.var.index]
adata.var = adata.var[['highly_variable','highly_variable_nbatches']]
adata = adata[:,hvg]

# Read Visium slides
slides = []
for dname in os.listdir(data_path):
    if not dname.startswith('.'):
        tmp = sc.read_visium(os.path.join(path_slides, dname,'outs'))
        tmp.var_names_make_unique()
        slides.append(tmp)
    
slides = slides[0].concatenate(slides[1:], join='inner')

# Filter intersection of genes
intersect = np.intersect1d(slides.var_names, adata.var_names)
adata = adata[:, intersect].copy()

# Set up scLVM model
CondSCVI.setup_anndata(adata, layer="counts", labels_key=label_name, batch_key=sample_id)
model = CondSCVI(adata, weight_obs=False)
model.view_anndata_setup()

# Train
model.train()

# Save model
model.save(path_output, overwrite=True)
history = model.history["elbo_train"]






