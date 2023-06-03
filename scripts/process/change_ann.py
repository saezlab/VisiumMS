import pandas as pd
import numpy as np
import scanpy as sc
import os

import argparse

# Read command line and set args
parser = argparse.ArgumentParser()
parser.add_argument('-a', '--path_ann_obj', help='Annotated AnnData object', required=True)
args = vars(parser.parse_args())

path_ann_obj = args['path_ann_obj']
###############################

# Read annotated adata
adata = sc.read_h5ad(path_ann_obj)

# Change cell_type names
ann_dict = {
   'Astros_c': 'Astros',
   'Astros_r': 'Astros',
   'Oligos_d': 'OPC'
}
adata.obs['cell_type'] = adata.obs['leiden'].copy()
adata.obs['leiden'] = [ann_dict[cell_type] if cell_type in ann_dict else cell_type for cell_type in adata.obs['leiden']]

# Save
adata.write(path_ann_obj)

