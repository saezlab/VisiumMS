import scanpy as sc
import pandas as pd
import numpy as np
import decoupler as dc
import os
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-s','--slide_path', required=True)
parser.add_argument('-g','--gmt_path', required=True)
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

slide_path = args['slide_path']
gmt_path = args['gmt_path']
out_path = args['out_path']

# Read data
slide = sc.read_h5ad(slide_path)
reactome = dc.read_gmt(gmt_path)
reactome['source'] = [s.split('REACTOME')[1].replace('_', ' ').lstrip() for s in reactome['source']]

# Keep relevant gene sets
msk = ~reactome['source'].str.contains('FETAL|INFECTION|SARS', case=False)
reactome = reactome[msk]

# Run enrichment analysis
dc.run_ulm(slide, reactome, weight=None, use_raw=False)

# Save
slide.obsm['ulm_estimate'].to_csv(out_path)
