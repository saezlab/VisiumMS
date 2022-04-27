import scanpy as sc
import numpy as np

from scvi.model import CondSCVI, DestVI

import os
import argparse


# Read command line and set args
parser = argparse.ArgumentParser(prog='stlvm', description='Deconvolution with stLVM')
parser.add_argument('-m', '--model_path', help='Directory containing sclvm model', required=True)
parser.add_argument('-s', '--sample_path', help='Path to visium slide', required=True)
parser.add_argument('-o', '--output_path', help='Output path where to store the model', required=True)
args = vars(parser.parse_args())

model_path = args['model_path']
sample_path = args['sample_path']
sample_id = os.path.basename(sample_path)
output_path = args['output_path']

# Load sclvm model
sc_model = CondSCVI.load(model_path)

# Load individual slide
st_adata = sc.read_visium(os.path.join(sample_path, 'outs'))
st_adata.var_names_make_unique()
st_adata.layers["counts"] = st_adata.X.copy()
sc.pp.normalize_total(st_adata, target_sum=10e4)
sc.pp.log1p(st_adata)
st_adata.raw = st_adata

# Intersect genes
intersect = np.intersect1d(sc_model.adata.var_names, st_adata.var_names)
st_adata = st_adata[:, intersect].copy()

# Set up model
DestVI.setup_anndata(st_adata, layer="counts")
st_model = DestVI.from_rna_model(st_adata, sc_model)
st_model.view_anndata_setup()

# Train
st_model.train(max_epochs=2500)

# Save model
st_model.save(os.path.join(output_path, sample_id), overwrite=True, save_anndata=True)
st_model.history["elbo_train"].to_csv(os.path.join(output_path, sample_id, 'history.csv'))

