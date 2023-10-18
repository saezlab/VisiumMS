import numpy as np
import pandas as pd
import anndata as ad
from anndata import AnnData
import scanpy as sc
import argparse
import liana as li
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-c','--vs_ctrl_path', required=True)
parser.add_argument('-p','--sn_pw_path', required=True)
parser.add_argument('-m','--meta_path', required=True)
parser.add_argument('-t','--thr_score', required=True)
parser.add_argument('-a','--thr_padj', required=True)
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

vs_ctrl_path = args['vs_ctrl_path']
sn_pw_path = args['sn_pw_path']
meta_path = args['meta_path']
thr_score = float(args['thr_score'])
thr_padj = float(args['thr_padj'])
out_path = args['out_path']

def get_names(pair, contrasts, ctrl, sn_pw):
    a, b = pair.split('^')
    sub_ctrl = ctrl[ctrl['pairs'] == pair]
    sub_phtw = sn_pw[sn_pw['cell_type'].str.contains('|'.join([a, b]))]
    lr_names = sub_ctrl['names'].unique().astype(str)
    # Extract consisten pathways between cell type pairs
    pw_names = []
    for contrast in contrasts.keys():
        tmp = sub_phtw[sub_phtw['contrast'] == contrast].pivot(index='cell_type', columns='source', values='score').dropna(axis=1).T
        if a in tmp.columns and b in tmp.columns:
            tmp = tmp.assign(sign=lambda x: (np.sign(x[a]) == np.sign(x[b])) * np.sign(x[a]))
            pw_names.extend(tmp[tmp['sign'] != 0].index)
    pw_names = np.unique(pw_names)
    pw_names = pw_names[np.isin(pw_names, pw_scores.var_names)]
    return lr_names, pw_names


def spatial_para(adata):
    from scipy.sparse import csr_matrix
    # Smooth by connectivity
    X = adata.X
    if isinstance(X, csr_matrix):
        X = X.A
    conn = adata.obsp['spatial_connectivities']
    X = (conn @ X) / conn.sum(axis=1).A  # Normalize by spatial weights
    return X


def get_corr(sub_lr_scores, sub_pw_scores):
    from scipy import stats
    corr_df = []
    for lr in sub_lr_scores.var_names:
        x = sub_lr_scores[:, lr].X.ravel()
        msk = ~np.isnan(x)
        if np.sum(msk) >= 5:
            for pw in sub_pw_scores.var_names:
                y = sub_pw_scores[:, pw].X.ravel()
                r, p = stats.pearsonr(x[msk], y[msk])
                corr_df.append([lr, pw, r, p])
    corr_df = pd.DataFrame(corr_df, columns=['ctlr', 'pw', 'r2', 'pval'])
    return corr_df


def read_adatas(sample_id):
    slide_path = 'data/prc/vs/{0}/adata.h5ad'.format(sample_id)
    ctlr_path = 'data/prc/vs/{0}/ctlr_scores.csv'.format(sample_id)
    pthw_path = 'data/prc/vs/{0}/reactome.csv'.format(sample_id)
    slide = sc.read_h5ad(slide_path)
    lr_scores = AnnData(pd.read_csv(ctlr_path, index_col=0), dtype=float)
    lr_scores.X[np.isnan(lr_scores.X)] = 0.
    lr_scores.obsm['spatial'] = slide.obsm['spatial'].copy()
    lr_scores.uns['spatial'] = slide.uns['spatial'].copy()
    pw_scores = AnnData(pd.read_csv(pthw_path, index_col=0), dtype=float)
    pw_scores.obsm['spatial'] = slide.obsm['spatial'].copy()
    pw_scores.uns['spatial'] = slide.uns['spatial'].copy()
    return slide, lr_scores, pw_scores

def get_pairs_corrs(pairs, contrasts, ctrl, sn_pw, thr_score=0.10):
    corrs = []
    for pair in tqdm(pairs):
        # Subset by cell type specific lr and pw
        lr_names, pw_names = get_names(pair, contrasts, ctrl, sn_pw)
        if lr_names.size > 0 and pw_names.size > 0:
            sub_lr_scores = lr_scores[:, lr_names].copy()
            sub_pw_scores = pw_scores[:, pw_names].copy()
            
            # Compute para
            li.ut.spatial_neighbors(sub_pw_scores, cutoff=0.1, bandwidth=150, set_diag=True)
            sub_pw_scores.X = spatial_para(sub_pw_scores)
            
            # Clip low values
            sub_lr_scores.X[sub_lr_scores.X < thr_score] = np.nan
            
            # Compute corr
            corr = get_corr(sub_lr_scores, sub_pw_scores)
            corrs.append(corr)
    if corrs:
        corrs = pd.concat(corrs)
        return corrs


lt_dict = {
    'Ctrl': {'CAvsCtrl': -1, 'CIvsCtrl': -1},
    'CA': {'CAvsCtrl': +1, 'CAvsCI': +1},
    'CI': {'CIvsCtrl': +1, 'CAvsCI': -1}
}

# Read func results for sn
all_sn_pw = pd.read_csv(sn_pw_path)
all_ctrl = pd.read_csv(vs_ctrl_path)

meta = pd.read_csv(meta_path)
vs_samples = meta[~meta['Batch vs'].isnull()]['Sample id'].values.astype('U')
corrs = []
for sample_id in vs_samples:
    print(sample_id)
    # Read slide
    slide, lr_scores, pw_scores = read_adatas(sample_id)
    
    # Find contrasts
    lesion_type = slide.obs['Lesion type'].unique()[0]
    contrasts = lt_dict[lesion_type]
    
    # Filter by conditions
    ctrl = all_ctrl[all_ctrl['contrast'].str.contains('|'.join(contrasts.keys()))].copy()
    sn_pw = all_sn_pw[all_sn_pw['contrast'].str.contains('|'.join(contrasts.keys()))].copy()
    ctrl = ctrl[ctrl['pvals_adj'] < thr_padj]
    sn_pw = sn_pw[sn_pw['adj_pvals'] < thr_padj]
    
    # Keep coherent signs
    ctrl[['source', 'target']] = ctrl['names'].str.split('^', n=2, expand=True, regex=False)[[0, 1]]
    ctrl['keep'] = [contrasts[x] == s for x,s in zip(ctrl['contrast'], ctrl['sign'])]
    ctrl = ctrl[ctrl['keep']].drop(columns='keep')
    
    # Identify unique pairs
    ctrl['pairs'] = ['^'.join(np.sort(np.unique([a, b]))) for a, b in zip(ctrl['source'], ctrl['target'])]
    ctrl['pairs'] = [x if '^' in x else "^".join([x, x]) for x in ctrl['pairs']]
    pairs = ctrl['pairs'].unique()
    
    corr  = get_pairs_corrs(pairs, contrasts, ctrl, sn_pw, thr_score=thr_score)
    if corr is not None:
        corr['sample_id'] = sample_id
        corr['lesion_type'] = lesion_type
        corrs.append(corr)
corrs = pd.concat(corrs)

# Summarize pws across slides
corrs = corrs.groupby(['ctlr', 'pw', 'lesion_type']).median(numeric_only=True).reset_index()
corrs = corrs[(corrs['pval'] < 0.05) & (corrs['r2'] > 0.1)].reset_index(drop=True)
corrs[['source', 'target']] = corrs['ctlr'].str.split('^', n=2, expand=True, regex=False)[[0, 1]]

# Make reference table
cell_types = all_sn_pw['cell_type'].unique().astype(str)
table_dict = dict()
for contrast in ['CAvsCtrl', 'CIvsCtrl', 'CAvsCI']:
    tble = (
        all_sn_pw[(all_sn_pw['adj_pvals'] < 0.15) & (all_sn_pw['contrast'] == contrast)]
        .pivot(index='source', columns='cell_type', values='score')
    )
    msk = ~np.isin(cell_types, tble.columns)
    tble.loc[:, cell_types[msk]] = np.nan
    table_dict[contrast] = tble

def subset_ref_ctype_table(df, source, target, cond):
    df = (
        df
        .loc[:, [source, target]]
        .dropna(how='all')
    )
    df.columns = ['source_cond_{0}'.format(cond), 'target_cond_{0}'.format(cond)]
    return df


def get_sn_dirs(inter, cond, source, target, pw, r2):
    i_corrs = corrs[(corrs['ctlr'] == inter) & (corrs['lesion_type'] == cond)]
    c_a, c_b = lt_dict[cond]
    tmp = pd.concat([
        subset_ref_ctype_table(table_dict[c_a], source, target, 'a'),
        subset_ref_ctype_table(table_dict[c_b], source, target, 'b'),
        ], axis=1)
    tmp = tmp.loc[[pw]]

    dir_c_a = np.all(np.sign(tmp[['source_cond_a', 'target_cond_a']].dropna(axis=1).values) == lt_dict[cond][c_a])
    dir_c_b = np.all(np.sign(tmp[['source_cond_b', 'target_cond_b']].dropna(axis=1).values) == lt_dict[cond][c_b])
    dir_bool = np.all([dir_c_a, dir_c_b])
    if np.sign(r2) < 0:
        dir_bool = ~dir_bool
        
    tmp['dir_bool'] = dir_bool
    
    tmp = tmp[['source_cond_a', 'source_cond_b', 'target_cond_a', 'target_cond_b', 'dir_bool']]
    return tmp.values[0]

# Chekc that direction of r2 and pws is consistent
cols = []
for row in tqdm(corrs.values):
    inter, pw, cond, r2, p, source, target = row
    cols.append(get_sn_dirs(inter, cond, source, target, pw, r2))
cols = pd.DataFrame(cols, columns=['source_cond_a', 'source_cond_b', 'target_cond_a', 'target_cond_b', 'dir_bool'])
cols = pd.concat([corrs, cols], axis=1)

# Filter by sign consistency
cols = cols[cols['dir_bool']].drop(columns=['dir_bool'])

# Write
cols.to_csv(out_path, index=False)
