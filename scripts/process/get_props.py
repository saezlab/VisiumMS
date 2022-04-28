import scanpy as sc
import numpy as np

from scvi.model import DestVI
import destvi_utils

import os
import argparse

parser = argparse.ArgumentParser(prog='get_prop', description='Gets proportions from deconv slide')
parser.add_argument('-s', '--slide_path', help='Directory containing deconvoluted slide', required=True)
args = vars(parser.parse_args())

slide_path = args['slide_path']

# Load slide
st_model = DestVI.load(slide_path)

# Get props
st_model.adata.obsm["proportions"] = st_model.get_proportions()

# Compute optimal props thresholds
thrs = destvi_utils.automatic_proportion_threshold(st_model.adata, kind_threshold='secondary', output_file=os.path.join(slide_path, 'thr.html'))

# Filter and recalculate props
st_model.adata.obsm["f_proportions"] = st_model.adata.obsm["proportions"].copy()
for cell_type in thrs.keys():
    thr = thrs[cell_type]
    vals = st_model.adata.obsm["f_proportions"][cell_type]
    st_model.adata.obsm["f_proportions"][cell_type] = vals * (vals > thr)
st_model.adata.obsm["f_proportions"] = st_model.adata.obsm["f_proportions"] / st_model.adata.obsm["f_proportions"].sum(1).values.reshape(-1,1)

# Save
st_model.adata.write(os.path.join(slide_path, 'adata.h5ad'))

