import scanpy as sc
import pertpy as pt
import numpy as np
import pandas as pd
from anndata import AnnData

from composition_stats import closure
from composition_stats import clr

import decoupler as dc

import os
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-m','--meta_path', required=True)
parser.add_argument('-a', '--ann_path', required=True)
parser.add_argument('-s','--cstates_path', required=True)
parser.add_argument('-p','--plot_path', required=True)
parser.add_argument('-t','--table_path', required=True)
args = vars(parser.parse_args())

meta_path = args['meta_path']
ann_path = args['ann_path']
cstates_path = args['cstates_path']
plot_path = args['plot_path']
table_path = args['table_path']

def get_credibility(adata, fdr=0.25):
    df = []
    for ctye in adata.var_names:
        print(ctye)
        sccoda_model = pt.tl.Sccoda()
        sccoda_data = sccoda_model.prepare(
                adata,
                formula="condition",
                reference_cell_type=ctye
            )
        sccoda_model.run_nuts(sccoda_data, rng_key=42)
        sccoda_model.set_fdr(sccoda_data, est_fdr=fdr)
        df.append(sccoda_model.credible_effects(sccoda_data).reset_index().assign(ref_ctype=ctye))
    df = pd.concat(df)
    df = (
        df
        .groupby('Cell Type')
        [['Final Parameter']]
        .sum()
        .reset_index()
        .assign(pct_credible=lambda x: x['Final Parameter'] / sccoda_data.shape[1])
        .drop(columns='Final Parameter')
    )
    return df

def contrast_compositions(adata, groupby, contrasts, fdr=0.25):
    df = []
    for contrast in contrasts:
        print(contrast)
        # Subset
        cntr = adata[adata.obs[groupby].str.contains('{0}|{1}'.format(*contrast))].copy()
    
        # Define reference
        cntr.obs['condition'] = pd.Categorical(cntr.obs[groupby].astype('str'), categories=contrast)

        # Check min samples
        if (cntr.obs['condition'].value_counts() >= 3).all():

            # Filter noisy features
            keep_f = (np.sum(cntr.X != 0, axis=0) / cntr.shape[0]) > 0.5
            cntr = cntr[:, keep_f].copy()

            # Check credibility
            cred = get_credibility(cntr, fdr=fdr)
            cred = cred[cred['pct_credible'] >= 0.5]

            # If enough to test
            if cred.shape[0] > 0:

                # Compute props and clr transformation
                cntr.X = clr(closure(cntr.X + 0.5))
                cntr = cntr[:, cred['Cell Type']].copy()

                # Compute t-test
                res = dc.rank_sources_groups(cntr, groupby='condition', reference=contrast[0], method='t-test')
                res = res.loc[res['group'] == contrast[1]]
                df.append(res)
    if len(df) > 0:
        df = pd.concat(df)
    else:
        df = None
    return df


# Read sn
meta = pd.read_csv(meta_path)
ctypes = sc.read_h5ad(ann_path).obs.copy()
ctypes = (
    ctypes
    .reset_index(drop=True)
    .groupby(['Sample id', 'leiden'])
    .count()
    [['total_counts']]
    .reset_index()
    .pivot(index='Sample id', columns='leiden', values='total_counts')
)
ctypes = AnnData(X=ctypes, obs=meta.set_index('Sample id').loc[ctypes.index], dtype=np.float64)

# Read vs
vs_samples = meta[~meta['Batch vs'].isnull()]['Sample id'].values.astype('U')

niches = []
for vs_sample in vs_samples:
    sample_niches = pd.read_csv('data/prc/vs/{0}/niches.csv'.format(vs_sample), index_col=0)
    sample_niches['Sample id'] = vs_sample
    sample_niches['Lesion type'] = meta.set_index('Sample id').loc[vs_sample, 'Lesion type']
    niches.append(sample_niches)
niches = pd.concat(niches)
niches = (
    niches
    .groupby(['Sample id', 'leiden'])
    .count()
    [['Lesion type']]
    .reset_index()
    .pivot(index='Sample id', columns='leiden', values='Lesion type')
    .fillna(0)
)
niches = AnnData(X=niches, obs=meta.set_index('Sample id').loc[niches.index], dtype=np.float64)
niches.var['n_cells'] = niches.X.sum(0)
niches = niches[niches.obs['Lesion type'].str.contains('CA|CI')].copy()
niches = niches[:, niches.var_names.str.contains('LC|LR|PPWM|VI')].copy()
niches.obs = (
    niches.obs
    .reset_index()
    [['Lesion type', 'Sample id']]
    .assign(scCODA_sample_id=lambda x: x['Sample id'])
    .set_index('scCODA_sample_id')
)
niches = niches[niches.obs.sort_values('Lesion type').index].copy()

# Plot
g1 = pt.pl.coda.boxplots(
    ctypes,
    feature_name="Lesion type",
    add_dots=True,
    plot_facets=True,
)

g2 = pt.pl.coda.boxplots(
    niches,
    feature_name="Lesion type",
    add_dots=True,
    plot_facets=True
)

# Save to pdf
pdf = matplotlib.backends.backend_pdf.PdfPages(plot_path)
for fig in [g1.fig, g2.fig]:
    pdf.savefig(fig, bbox_inches='tight')
pdf.close()

# Coda
res_sn = contrast_compositions(
    ctypes,
    groupby='Lesion type',
    contrasts=[['Ctrl', 'CA'], ['Ctrl', 'CI'], ['CA', 'CI']],
    fdr=0.25
)
res_sn['category'] = 'sn'

res_vs = contrast_compositions(
    niches,
    groupby='Lesion type',
    contrasts=[['CA', 'CI']],
    fdr=0.25
)
res_vs['category'] = 'vs'

res = [res_sn, res_vs]
for ctype in ctypes.var_names:
    obs = pd.read_csv('data/prc/ctypes/{0}_ann.csv'.format(ctype), index_col=0)
    obs['sample_id'] = [i.split('-')[0] for i in obs.index]
    obs['n'] = 1
    obs = (
        obs
        .groupby(['sample_id', 'cell_states'])
        .count()
        .reset_index()
        .pivot(index='sample_id', columns='cell_states', values='n')
    ).fillna(0)

    # Remove NA
    vars = obs.columns[~obs.columns.str.contains('_NA')]
    obs = AnnData(X=obs.loc[:, vars].copy(), obs=meta.set_index('Sample id').loc[obs.index])

    # Run comospition analysis
    res_obs = contrast_compositions(
        obs,
        groupby='Lesion type',
        contrasts=[['Ctrl', 'CA'], ['Ctrl', 'CI'], ['CA', 'CI']],
        fdr=0.25
    )
    if res_obs is not None:
        res_obs['category'] = 'cs'
        res.append(res_obs)
res = pd.concat(res)

# Write
res.to_csv(table_path, index=False)

