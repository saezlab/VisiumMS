import numpy as np
import pandas as pd
import decoupler as dc
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--inp_path', required=True)
parser.add_argument('-g','--gmt_path', required=True)
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

inp_path = args['inp_path']
gmt_path = args['gmt_path']
out_path = args['out_path']

# Read data
df = pd.read_csv(inp_path, index_col=0)
reactome = dc.read_gmt(gmt_path)
reactome['source'] = [s.split('REACTOME')[1].replace('_', ' ').lstrip() for s in reactome['source']]

# Keep relevant gene sets
msk = ~reactome['source'].str.contains('FETAL|INFECTION|SARS', case=False)
reactome = reactome[msk]

# Drop duplicates
reactome = reactome.drop_duplicates(['source', 'target'])

# Get unique values
contrasts = df['contrast'].unique()
cell_types = df['cell_type'].unique()

# Iterate and compute enrichment score with ulm
res = []
for contrast in contrasts:
    for cell_type in cell_types:
        cdf = df[(df['cell_type'] == cell_type) & (df['contrast'] == contrast)].copy()
        if cdf.shape[0] > 0:
            print(contrast, cell_type)
            cdf = cdf.reset_index(names='genes').pivot(index='cell_type', columns='genes', values='stat')
            acts, pvals = dc.run_ulm(cdf, reactome, weight=None)
            cdf = (
                dc.melt([acts, pvals])
                .rename(columns={'sample': 'cell_type'})
                .drop(columns=['method'])
            )
            cdf['adj_pvals'] = dc.p_adjust_fdr(cdf['pvals'].values)
            cdf['contrast'] = contrast
            cdf = cdf.sort_values('adj_pvals')
            cdf = cdf[['contrast', 'cell_type', 'source', 'score', 'pvals', 'adj_pvals']]
            res.append(cdf)
res = pd.concat(res).reset_index(drop=True)

# Write
res.to_csv(out_path, index=False)

