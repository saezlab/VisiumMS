import pandas as pd
import numpy as np
import scanpy as sc
import os
from scipy import stats

import seaborn as sns
import matplotlib.pyplot as plt
import decoupler as dc

# Define path
fig_path = 'figures/manuscript/fig3/'
plt.rcParams['font.sans-serif'] = 'Arial'


# Read results
res = pd.read_csv('figures/manuscript/fig3/microglia_inters.csv')
res['receptor.complex'] = [r.split('_') for r in res['receptor.complex']]
res = res.explode('receptor.complex')
deg_net = pd.read_csv('figures/manuscript/fig3/markers.csv').groupby('group').head(100)

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


def compute_lr_corrs(df, sample_id, deg):
    
    # Read slide
    slide = read_slide(sample_id)
    dc.run_ulm(slide, deg, source='group', target='names', weight=None, use_raw=False)
    slide.obsm['props'].loc[:, '2'] = slide.obsm['ulm_estimate'].loc[:, '2']
    slide.obsm['props'].loc[:, '3'] = slide.obsm['ulm_estimate'].loc[:, '3']

    # Make a copy
    df = df.copy()
    min_prop = 1 / slide.obsm['props'].shape[1]
    corrs = []
    pvals = []
    for s, t, l, r in zip(df['source'], df['target'], df['ligand.complex'], df['receptor.complex']):

        # Filter by min_prop for both cells
        msk_s = slide.obsm['props'][s].values > min_prop
        msk_t = slide.obsm['props'][t].values > min_prop
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

# Run per sample
meta = pd.read_csv('data/metadata.csv').set_index('sample_id')
msk = np.isin(meta['lesion_type'], ['Control', 'Chronic Active'])
meta = meta.loc[msk]
dfs = []
for sample_id in meta.index:
    dfs.append(compute_lr_corrs(res, sample_id, deg_net))
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
    .sort_values('diff', ascending=False)
)

os.makedirs(fig_path, exist_ok=True)
mean_corrs.to_csv(os.path.join(fig_path, 'corrs_inters.csv'), index=False)

# Example
sample_id = 'MS377T'
slide = read_slide(sample_id)
fig, axes = plt.subplots(1, 3, figsize=(9, 3), facecolor='white', tight_layout=True)
sc.pl.spatial(slide, color='FTL', size=1.6, frameon=False, ax=axes[0], return_fig=False, show=False)
sc.pl.spatial(slide, color='APOE', size=1.6, frameon=False, ax=axes[1], return_fig=False, show=False)
sc.pl.spatial(slide, color='TREM2', size=1.6, frameon=False, ax=axes[2], return_fig=False, show=False)
fig.savefig(os.path.join(fig_path, 'spatial_eg_inter.pdf'), bbox_inches='tight')
