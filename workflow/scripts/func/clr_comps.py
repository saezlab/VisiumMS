import anndata as ad
from composition_stats import closure, clr
from scipy import stats
import decoupler as dc
import scanpy as sc
import numpy as np
import pandas as pd
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--ann_path', required=True)
parser.add_argument('-m','--meta_path', required=True)
parser.add_argument('-k','--ks_path', required=True)
parser.add_argument('-w','--wt_path', required=True)
args = vars(parser.parse_args())

ann_path = args['ann_path']
meta_path = args['meta_path']
thr_padj = 0.05
ks_path = args['ks_path']
wt_path = args['wt_path']

adata = ad.read_h5ad(ann_path)
meta = pd.read_csv(meta_path)

def get_cnts(df, group_col):
    msk = [(not s.endswith('_NA')) and (not s.startswith('NEU')) for s in df[group_col]]
    cnts = (
        df
        .loc[msk]
        .reset_index()
        .groupby(['Sample id', group_col])
        .count()
        [['index']]
        .reset_index()
        .pivot(index='Sample id', columns=group_col, values='index')
        .fillna(0)
    )
    cnts.loc[:, :] = closure(cnts + 0.5)
    cnts = ad.AnnData(
        cnts.copy(),
        dtype=float,
    )
    cnts.layers['clr'] = clr(cnts.X)
    return cnts

def test_ks(df):
    ks = []
    ltypes = df.obs['Lesion type'].unique()
    for ctype in df.var_names:
        vars = [df[df.obs['Lesion type'] == ltype, ctype].layers['clr'].ravel().copy() for ltype in ltypes]
        s, p = stats.kruskal(*vars)
        ks.append([ctype, s, p])
    ks = pd.DataFrame(ks, columns=['ctype', 'stat', 'pval'])
    ks['padj'] = dc.p_adjust_fdr(ks['pval'].values)
    ks = ks.sort_values('padj')
    return ks

def test_wt(df, groups, contrasts=[['CA', 'Ctrl'], ['CI', 'Ctrl'], ['CA', 'CI']]):
    wt = []
    for ctype in groups:
        for contrast in contrasts:
            x, y = [df[df.obs['Lesion type'] == ltype, ctype].layers['clr'].ravel().copy() for ltype in contrast]
            s, p = stats.ranksums(x, y)
            wt.append([ctype, contrast[0], contrast[1], s, p])
    wt = pd.DataFrame(wt, columns=['ctype', 'group', 'ref', 'stat', 'pval'])
    wt['padj'] = dc.p_adjust_fdr(wt['pval'].values)
    wt = wt.sort_values('padj')
    return wt

def read_ctype(ctype, group_col):
    df = pd.read_csv('data/prc/ctypes/{ctype}_ann.csv'.format(ctype=ctype), index_col=0).reset_index()
    df[['Sample id', 'barcode']] = df['index'].str.split('-', n=1, expand=True)
    df = get_cnts(df, group_col=group_col)
    return df

def test_ctype(ctype, meta, thr):
    df = read_ctype(ctype, group_col='cell_states')
    df.obs = pd.merge(df.obs.reset_index(names='Sample id'), meta).set_index('Sample id')
    ks = test_ks(df)
    wt = test_wt(df, groups=ks[ks['padj'] < thr]['ctype'].values)
    return ks, wt

def read_niches(meta, group_col):
    vs_samples = meta[~meta['Batch vs'].isnull()]['Sample id'].values.astype('U')
    df = []
    for sample_id in vs_samples:
        tmp = pd.read_csv('data/prc/vs/{sample_id}/niches.csv'.format(sample_id=sample_id), index_col=0)
        # Remove ependym and GM
        tmp = tmp.loc[~np.isin(tmp[group_col].values.astype('U'), ['Ependym', 'GM']), :]
        tmp['Sample id'] = sample_id
        df.append(tmp)
    df = pd.concat(df)
    df = get_cnts(df, group_col=group_col)
    return df

def read_niches_prop(meta, groups_col='leiden'):
    if groups_col == 'leiden':
        msk = meta['Condition'] == 'MS'
    else:
        msk = ~meta['Batch vs'].isnull()
    vs_samples = meta[(~meta['Batch vs'].isnull()) & msk]['Sample id'].values.astype('U')
    df = []
    for sample_id in vs_samples:
        # Remove ependym, GM and neurons
        obs = pd.read_csv('data/prc/vs/{sample_id}/niches.csv'.format(sample_id=sample_id), index_col=0)
        obs['Sample id'] = sample_id
        msk = ~ np.isin(obs['leiden'].values.astype('U'), ['Ependym', 'GM'])
        tmp = pd.read_csv('data/prc/vs/{sample_id}/abunds.csv'.format(sample_id=sample_id), index_col=0)
        tmp = tmp.drop(columns=['NEU'])
        tmp, obs = tmp.loc[msk, :].copy(), obs.loc[msk].copy()
        # Transform to anndata
        tmp = ad.AnnData(tmp, obs=obs, dtype=float)
        # Sum per niche
        tmp = dc.get_pseudobulk(
            adata=tmp,
            sample_col='Sample id',
            groups_col=groups_col,
            mode='sum',
            min_cells=0,
            min_counts=0,
            skip_checks=True,
        )
        df.append(tmp)
    # Merge
    df = ad.concat(df)
    df.obs = pd.merge(df.obs.reset_index(), meta[['Sample id', 'Lesion type']]).set_index('index')
    # Compute prop and clr
    df.X = closure(df.X + 0.5)
    df.layers['clr'] = clr(df.X)
    return df

def test_wt_niches(df, groups, contrasts=[['CA', 'CI']]):
    wt = []
    niches = ['PPWM', 'LR', 'LC', 'VI']
    for ctype in groups:
        for niche in niches:
            for contrast in contrasts:
                x, y = [df[(df.obs['Lesion type'] == ltype) & (df.obs['leiden'] == niche), ctype].layers['clr'].ravel().copy() for ltype in contrast]
                s, p = stats.ranksums(x, y)
                wt.append([ctype, niche, contrast[0], contrast[1], s, p])
    wt = pd.DataFrame(wt, columns=['ctype', 'niche', 'group', 'ref', 'stat', 'pval'])
    wt['padj'] = dc.p_adjust_fdr(wt['pval'].values)
    wt = wt.sort_values('padj')
    return wt

# For cell types
df = adata.obs[~np.isin(adata.obs['leiden'].values, ['BC', 'SC'])].copy()
df['leiden'] = df['leiden'].astype(str)
df = get_cnts(df, group_col='leiden')
df.obs = pd.merge(df.obs.reset_index(names='Sample id'), meta).set_index('Sample id')
ks_df = test_ks(df)
wt_df = test_wt(df, groups=ks_df[ks_df['padj'] < thr_padj]['ctype'].values)
ks_df['type'], wt_df['type'] = 'sn', 'sn'

# For cell states
for ctype in df.var_names:
    ks, wt = test_ctype(ctype, meta, thr=thr_padj)
    ks['type'], wt['type'] = 'cs', 'cs'
    ks_df = pd.concat([ks_df, ks])
    wt_df = pd.concat([wt_df, wt])

# For niches in ST
df = read_niches(meta, group_col='leiden')
df.obs = pd.merge(df.obs.reset_index(names='Sample id'), meta).set_index('Sample id')
df = df[df.obs['Condition'] == 'MS'].copy()
ks = test_ks(df)
wt = test_wt(df, groups=ks[ks['padj'] < thr_padj]['ctype'].values, contrasts=[['CA', 'CI']])
ks['type'], wt['type'] = 'ns_st', 'ns_st'
ks_df = pd.concat([ks_df, ks])
wt_df = pd.concat([wt_df, wt])

# For cell types in ST at niche level
df = read_niches_prop(meta)
ks = test_ks(df)
wt = test_wt_niches(df, groups=ks[ks['padj'] < thr_padj]['ctype'].values, contrasts=[['CA', 'CI']])
ks['type'], wt['type'] = 'ct_ns_st', 'ct_ns_st'
ks_df = pd.concat([ks_df, ks])
wt_df = pd.concat([wt_df, wt])

# For cell types in ST at whole slide
df = read_niches_prop(meta, groups_col=None)
ks = test_ks(df)
wt = test_wt(df, groups=ks[ks['padj'] < thr_padj]['ctype'].values)
ks['type'], wt['type'] = 'ct_st', 'ct_st'
ks_df = pd.concat([ks_df, ks])
wt_df = pd.concat([wt_df, wt])

# Filter
ks_df = ks_df[ks_df['padj'] < thr_padj]
wt_df = wt_df[wt_df['padj'] < thr_padj]

# Save
ks_df.to_csv(ks_path, index=False)
wt_df.to_csv(wt_path, index=False)