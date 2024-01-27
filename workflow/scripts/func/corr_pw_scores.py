import pandas as pd
import numpy as np
import decoupler as dc
import scanpy as sc
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--ctlr_path', required=True)
parser.add_argument('-b','--meta_path', required=True)
parser.add_argument('-c','--out_path', required=True)
args = vars(parser.parse_args())

ctlr_path = args['ctlr_path']
meta_path = args['meta_path']
out_path = args['out_path']

# Read
meta = pd.read_csv(meta_path)
res = pd.read_csv(ctlr_path)

def read_slide(sample_id):
    slide = sc.read_h5ad('data/prc/vs/{0}/adata.h5ad'.format(sample_id))
    slide.obsm['score'] = pd.read_csv('data/prc/vs/{0}/ctlr_scores.csv'.format(sample_id), index_col=0)
    slide.obsm['pathway'] = pd.read_csv('data/prc/vs/{0}/reactome.csv'.format(sample_id), index_col=0)
    return dc.get_acts(slide, 'score'), dc.get_acts(slide, 'pathway')


def compute_corr(sample_id, inters):
    # Open
    score, pathw = read_slide(sample_id)
    score = score[:, inters].copy()
    # Compute correlation
    corr = np.corrcoef(score.X, pathw.X, rowvar=False)[: score.shape[1], score.shape[1] :]
    corr = pd.DataFrame(corr, index=score.var_names, columns=pathw.var_names)
    corr = corr.reset_index().melt(id_vars='index').dropna().sort_values('value', ascending=False).reset_index(drop=True).rename(columns={'index': 'inter', 'variable': 'pathway', 'value': 'corr'})
    corr['sample_id'] = sample_id
    return corr

# Extract
vs_samples = meta[~meta['Batch vs'].isnull()]['Sample id'].values.astype('U')
inters = res['names'].unique().astype('U')

# Compute
corr = []
for sample_id in vs_samples:
    print(sample_id)
    corr.append(compute_corr(sample_id, inters))
corr = pd.concat(corr)

# Write
corr.to_csv(out_path, index=False)
