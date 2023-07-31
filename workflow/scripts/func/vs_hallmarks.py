import scanpy as sc
import pandas as pd
import numpy as np
import decoupler as dc
import os
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-s','--slide_path', required=True)
parser.add_argument('-g','--hallm_path', required=True)
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

slide_path = args['slide_path']
hallm_path = args['hallm_path']
out_path = args['out_path']

# Read data
slide = sc.read_h5ad(slide_path)
hallmarks = pd.read_csv(hallm_path)

# Run enrichment analysis
dc.run_ulm(slide, hallmarks, weight=None, use_raw=False)

# Save
slide.obsm['ulm_estimate'].to_csv(out_path)

