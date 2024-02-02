import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import liana as li
import decoupler as dc
from scipy.sparse import csr_matrix


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--meta_path', required=True)
parser.add_argument('-c','--plot_path', required=True)
args = vars(parser.parse_args())

meta_path = args['meta_path']
plot_path = args['plot_path']


def read_annots(sample_id):
    slide = sc.read_h5ad('data/prc/vs/{0}/adata.h5ad'.format(sample_id))
    slide.obs['niches'] = pd.read_csv('data/prc/vs/{0}/niches.csv'.format(sample_id), index_col=0)
    return slide

def extract_niches(slide):
    tmp = slide.obs.pivot(columns='niches', values='Age')
    tmp.values[~np.isnan(tmp)] = 1
    slide.obsm['niches'] = tmp.fillna(0).copy()
    return dc.get_acts(slide, 'niches')

def update_w(adata):
    w = adata.obsp['spatial_connectivities'].A
    w[w != 0] = 1
    w = csr_matrix(w / w.sum(1))
    return w

def get_score(adata, col):
    w = update_w(adata)
    score = w @ adata.X
    score = pd.DataFrame(score, index=adata.obs_names, columns=adata.var_names)
    adata.obsm['score'] = score
    score = dc.get_acts(adata, 'score')
    score = dc.get_pseudobulk(
        score,
        sample_col=col,
        groups_col=None,
        mode='mean',
        min_cells=0,
        min_counts=0
    )
    return score


meta = pd.read_csv(meta_path)
vs_samples = meta[~meta['Batch vs'].isnull()]['Sample id'].values.astype('U')
lst = ['PPWM', 'LR', 'LC', 'VI']

dfs = []
for sample_id in vs_samples:
    if sample_id.startswith('CO'):
        continue
    print(sample_id)
    slide = read_annots(sample_id)
    slide = slide[slide.obs['niches'].isin(lst)].copy()
    niches = extract_niches(slide)
    li.ut.spatial_neighbors(niches, bandwidth=150, cutoff=0.1, kernel='gaussian', set_diag=True)
    score = get_score(niches, 'niches')
    df = (
        score
        .to_df()
        .reset_index()
        .melt(id_vars='index')
        .rename(columns={'index': 'source', 'niches': 'target'})
    )
    df['sample_id'] = sample_id
    dfs.append(df)
dfs = pd.concat(dfs)

msk = dfs['source'].isin(lst) & dfs['target'].isin(lst)
df = dfs.loc[msk].groupby(['source', 'target'])['value'].mean().reset_index()
df = dfs.groupby(['source', 'target'])['value'].mean().reset_index()
tmp = df.pivot(index='source', columns='target', values='value').fillna(0).loc[lst, lst]
fig, ax = plt.subplots(1, 1, figsize=(4, 3), dpi=150)
sns.heatmap(tmp, cmap='Blues', annot=True, ax=ax, cbar_kws={"shrink": 0.5, 'label': 'Mean proportion'})

# Save
fig.savefig(plot_path, bbox_inches='tight')
