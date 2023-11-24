import pandas as pd
import numpy as np
import scanpy as sc
import decoupler as dc
import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-m','--meta_path', required=True)
parser.add_argument('-p','--plot_path', required=True)
args = vars(parser.parse_args())

meta_path = args['meta_path']
plot_path = args['plot_path']


def read_slide(sample_id):
    slide = sc.read_h5ad('data/prc/vs/{0}/adata.h5ad'.format(sample_id))
    slide.obsm['props'] = pd.read_csv('data/prc/vs/{0}/props.csv'.format(sample_id), index_col=0)
    slide.obs['niches'] = pd.read_csv('data/prc/vs/{0}/niches.csv'.format(sample_id), index_col=0)
    slide.obs_names = [o + '|' + sample_id for o in slide.obs_names]
    return slide

# Read meta
meta = pd.read_csv(meta_path)

# Read slides
vs_samples = meta[~meta['Batch vs'].isnull()]['Sample id'].values.astype('U')
adata = []
for vs_sample in vs_samples:
    print(vs_sample)
    slide = read_slide(vs_sample)
    adata.append(slide)
adata = ad.concat(adata, join='outer')

niches = ['GM', 'WM', 'PPWM', 'LR', 'LC', 'VI']
adata = adata[np.isin(adata.obs['niches'].astype('U'), niches)].copy()

adata.obs['niches'] = pd.Categorical(adata.obs['niches'], categories=niches)

markers = dict(
    GM=['SYT1', 'CEND1', 'BRSK2'],
    WM=['MBP', 'PLP1', 'CNP'],
    PPWM=['SUN2', 'BOK', 'SOX10'],
    LR=['CSF1R', 'APOC1', 'CHI3L1'],
    LC=['GJA1', 'SORBS1', 'DTNA'],
    VI=['VIM', 'CLDN5', 'VWF']
)
fig1 = sc.pl.dotplot(
    adata,
    markers,
    groupby='niches',
    standard_scale='var',
    dendrogram=False,
    categories_order=niches,
    figsize=(9, 3),
    cmap='Reds',
    size_title='Fraction of spots\n in group (%)',
    return_fig=True
)

fig1.show()
fig1.fig.set_dpi(150)
fig1 = fig1.fig

vars = ['NEU', 'OPC', 'OL', 'MG', 'AS', 'TC', 'SC', 'EC', 'BC']

pdata = dc.get_pseudobulk(
    adata=dc.get_acts(adata, 'props'),
    sample_col='Sample id',
    groups_col='niches',
    mode='mean',
    min_cells=0,
    min_counts=0
)

df = dc.rank_sources_groups(
    adata=pdata,
    groupby='niches',
    method='wilcoxon'
)
df = df[(df['pvals_adj'] < 0.05) & (df['statistic'] > 0)].reset_index(drop=True)

fig2 = sc.pl.matrixplot(
    pdata,
    vars,
    groupby='niches',
    cmap='Purples',
    standard_scale='var',
    categories_order=niches,
    figsize=(6, 3),
    colorbar_title='Scaled mean proportion',
    return_fig=True
)

fig2.show()
fig2.fig.set_dpi(150)
fig2 = fig2.fig

# Save to pdf
pdf = matplotlib.backends.backend_pdf.PdfPages(plot_path)
for fig in [fig1, fig2]:
    pdf.savefig(fig, bbox_inches='tight')
pdf.close()
