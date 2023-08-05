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
parser.add_argument('-p', '--snp_path', required=True)
parser.add_argument('-o', '--out_path', required=True)
args = vars(parser.parse_args())

slide_path = args['slide_path']
gmt_path = args['gmt_path']
snp_path = args['snp_path']
out_path = args['out_path']

# Read data
slide = sc.read_h5ad(slide_path)

# Read sign pathways
#snp = pd.read_csv(snp_path)
#snp = (
#    snp[(snp['adj_pvals'] < 0.05) & (snp['score'] > 0.)]
#    .groupby(['contrast', 'cell_type'])
#    .head(25)
#    ['source']
#    .unique()
#    .astype('U')
#)

# Read reactome and filter by sign
gmt = dc.read_gmt(gmt_path)
gmt['source'] = [s.split('HALLMARK')[1].replace('_', ' ').lstrip() for s in gmt['source']]
#reactome = reactome[np.isin(reactome['source'].values.astype('U'), snp)]

# Run enrichment analysis
dc.run_ulm(slide, gmt, weight=None, use_raw=False)

# Save
slide.obsm['ulm_estimate'].to_csv(out_path)
