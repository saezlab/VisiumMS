import pandas as pd
import numpy as np
import scanpy as sc

from scipy import stats

import seaborn as sns
import matplotlib.pyplot as plt

# Read results
res_ch = pd.read_csv('data/prc/liana/ChronicActive.csv')
res_ct = pd.read_csv('data/prc/liana/Control.csv')
res = pd.concat([res_ch, res_ct])
res = res.drop_duplicates(['source', 'target', 'ligand.complex', 'receptor.complex'])
pos = pd.read_csv('data/prc/sign/deg/pos.csv')
neg = pd.read_csv('data/prc/sign/deg/neg.csv')
deg = pd.concat([pos, neg])
cell_types = deg.groupby('contrast')['name'].apply(lambda x: set(x))


def filter_by_deg(df, cell_types):
    df = df.copy()
    df['receptor.complex'] = [r.split('_') for r in df['receptor.complex']]
    df = df.explode('receptor.complex')
    msk = []
    for s, t, l, r in zip(df['source'], df['target'], df['ligand.complex'], df['receptor.complex']):
        m = True
        if s in cell_types.index:
            if l not in cell_types.loc[s]:
                m = False
        if t in cell_types.index:
            if r not in cell_types.loc[t]:
                m = False
        msk.append(m)
    msk = np.array(msk)
    return df.loc[msk]

def read_slide(sample_id):

    # Read rna-seq
    slide = sc.read_visium('data/raw/visium/{0}/outs/'.format(sample_id))
    slide.var_names_make_unique()
    sc.pp.filter_genes(slide, min_cells=3)
    sc.pp.filter_cells(slide, min_genes=200)

    # QC
    slide.var['mt'] = slide.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(slide, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # Store raw counts
    slide.layers['counts'] = slide.X

    # Normalize
    sc.pp.normalize_total(slide, target_sum=1e4)
    sc.pp.log1p(slide)

    # Read props
    props = pd.read_csv('data/prc/visium/{0}/cell_props.csv'.format(sample_id), index_col=0)
    inter = slide.obs.index.intersection(props.index)
    slide.obsm['props'] = props.loc[inter]

    return slide


def compute_lr_corrs(df, sample_id):
    
    # Read slide
    slide = read_slide(sample_id)

    # Make a copy
    df = df.copy()
    min_prop = 1 / slide.obsm['props'].shape[1]
    corrs = []
    pvals = []
    for s, t, l, r in zip(df['source'], df['target'], df['ligand.complex'], df['receptor.complex']):

        # Filter by min_prop for both cells
        msk_s = slide.obsm['props'][s].values > min_prop
<<<<<<< HEAD
        msk_t = slide.obsm['props'][s].values > min_prop
=======
        msk_t = slide.obsm['props'][t].values > min_prop
>>>>>>> Added plotting scripts
        tmp = slide[msk_s & msk_t]

        # Compute corrs
        if l in tmp.var_names and r in tmp.var_names:
            x_l = tmp[:, l].X.A.ravel()
            x_r = tmp[:, r].X.A.ravel()
            c, p = stats.pearsonr(x_l, x_r)
        else:
            c, p = np.nan, np.nan
        corrs.append(c)
        pvals.append(p)

    # Format df
    df['corrs'] = corrs
    df['pvals'] = pvals
    df['sample_id'] = sample_id

    return df


# Filter by deg
res = filter_by_deg(res, cell_types)

# Run per sample
meta = pd.read_csv('data/metadata.csv').set_index('sample_id')
dfs = []
for sample_id in meta.index:
    dfs.append(compute_lr_corrs(res, sample_id))
dfs = pd.concat(dfs)
dfs['lesion_type'] = [meta.loc[s, 'lesion_type'] for s in dfs['sample_id']]

# Compute mean corrs
mean_corrs = (
    dfs
    .dropna()
    .groupby(['source', 'target', 'ligand.complex', 'receptor.complex', 'lesion_type'])[['corrs']]
    .mean()
    .reset_index()
    .set_index(['source', 'target', 'ligand.complex', 'receptor.complex'])
    .pivot(columns='lesion_type', values='corrs')
    .reset_index()
    .assign(diff=lambda x: x['Chronic Active'] - x['Control'])
    .dropna()
    .sort_values('diff')
)

f_mean_corrs = mean_corrs[np.abs(mean_corrs['diff']) > 0.05]

# Write
<<<<<<< HEAD
f_mean_corrs.to_csv('data/prc/liana/mapped_inters.csv', index=False)

=======
f_mean_corrs.to_csv('data/prc/liana/mapped_inters.csv', index=False)
>>>>>>> Added plotting scripts
