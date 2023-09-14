import numpy as np
import pandas as pd
import anndata as ad
from anndata import AnnData
import decoupler as dc
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-n','--sn_lr_path', required=True)
parser.add_argument('-m','--meta_path', required=True)
parser.add_argument('-p','--thr_adjpval', required=True)
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())


sn_lr_path = args['sn_lr_path']
meta_path = args['meta_path']
thr_adjpval = float(args['thr_adjpval'])
out_path = args['out_path']

# Read func results for sn
sn_lr = pd.read_csv(sn_lr_path)

# Subset by significant lr and pathways
sn_lr = sn_lr[sn_lr['interaction_padj'] < thr_adjpval]

# Finds signs
sn_lr['names'] = ['{0}^{1}|{2}^{3}'.format(s, t, l.split('_')[0], r.split('_')[0]) for i, s,t,l,r in \
                  sn_lr[['source', 'target', 'ligand_complex', 'receptor_complex']].itertuples()]
signs = sn_lr[['contrast', 'names', 'sign']].assign(names=lambda x: x['names'].str.replace('|', '^', regex=False)).drop_duplicates(['contrast', 'names'])

# Read meta
meta = pd.read_csv(meta_path)
vs_samples = meta[~meta['Batch vs'].isnull()]['Sample id'].values.astype('U')
meta = meta.set_index('Sample id')

# Gather results
scores = []
for sample_id in vs_samples:
    print(sample_id)
    slide = AnnData(pd.read_csv('data/prc/vs/{0}/ctlr_scores.csv'.format(sample_id), index_col=0), dtype=float)
    slide.obs['Sample id'] = sample_id
    slide.obs['Lesion type'] = meta.loc[sample_id, 'Lesion type']
    slide.obs_names = [sample_id + '|' + i for i in slide.obs_names]
    scores.append(slide)
scores = ad.concat(scores, join='outer')
scores.X[np.isnan(scores.X)] = 0.

# Compute mean scores
mean_scores = dc.get_pseudobulk(
    adata=scores,
    sample_col='Sample id',
    groups_col=None,
    mode='mean',
    min_cells=0,
    min_counts=0,
)

def run_test(adata, signs, contrast, thr_adjpval=0.15, method='wilcoxon'):
    cond_a, cond_b = contrast.split('vs')
    sub_mean_ctlr = adata[adata.obs['Lesion type'].str.contains('{0}|{1}'.format(cond_a, cond_b))].copy()
    sub_signs = signs[signs['contrast'] == contrast]
    sub_mean_ctlr = sub_mean_ctlr[:, sub_signs['names']].copy()
    res = dc.rank_sources_groups(
        sub_mean_ctlr,
        groupby='Lesion type',
        reference=cond_b,
        method=method,
    ).sort_values('pvals_adj')
    res = (
        pd.merge(res, sub_signs)
        .assign(both=lambda x: np.sign(x['meanchange']) == np.sign(x['sign']))
    )
    res = res[(res['both']) & (res['pvals_adj'] < thr_adjpval)]
    res = res.reset_index(drop=True).drop(columns=['group', 'reference', 'both'])
    res = res[['contrast', 'names', 'statistic', 'meanchange', 'pvals', 'pvals_adj', 'sign']]
    return res

# Run contrasts
res = []
for contrast in ['CAvsCtrl', 'CIvsCtrl', 'CAvsCI']:
    tmp = run_test(mean_scores, signs, contrast, thr_adjpval=thr_adjpval)
    if tmp.shape[0] > 0:
        res.append(tmp)
res = pd.concat(res)
res.to_csv(out_path, index=False)
