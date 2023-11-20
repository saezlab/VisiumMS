import scanpy as sc
import pandas as pd
import numpy as np
import decoupler as dc
import os
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--slide_path', required=True)
parser.add_argument('-g', '--gmt_path', required=True)
parser.add_argument('-n', '--db_name', required=True)
parser.add_argument('-o', '--out_path', required=True)
args = vars(parser.parse_args())

slide_path = args['slide_path']
gmt_path = args['gmt_path']
db_name = args['db_name']
out_path = args['out_path']

# Read data
slide = sc.read_h5ad(slide_path)

# Read reactome and filter by sign
if db_name != 'progeny':
    gmt = dc.read_gmt(gmt_path)
    gmt['source'] = [s.split(db_name)[1].replace('_', ' ').lstrip() for s in gmt['source']]
    weight = None
else:
    gmt = pd.read_csv(gmt_path).groupby('source', observed=True).head(1000)
    weight = 'weight'

# Run enrichment analysis
dc.run_ulm(slide, gmt, weight=weight, use_raw=False)

# Save
slide.obsm['ulm_estimate'].to_csv(out_path)
