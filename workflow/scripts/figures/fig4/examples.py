import pandas as pd
import numpy as np
from tqdm import tqdm
from matplotlib_venn import venn3, venn3_circles
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import networkx as nx
import seaborn as sns
import pandas as pd
import anndata as ad
from anndata import AnnData
import decoupler as dc
import scanpy as sc
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--ctlr_path', required=True)
parser.add_argument('-b','--sn_annot_path', required=True)
parser.add_argument('-c','--corr_path', required=True)
parser.add_argument('-d','--meta_path', required=True)
parser.add_argument('-e','--colors_dict', required=True)
parser.add_argument('-f','--plot_path', required=True)
args = vars(parser.parse_args())

ctlr_path = args['ctlr_path']
sn_annot_path = args['sn_annot_path']
corr_path = args['corr_path']
meta_path = args['meta_path']
colors_dict = args['colors_dict']
plot_path = args['plot_path']

# Get palette
palette = dict(item.split(':') for item in colors_dict.strip("'").split(';'))

# Read LR results
res = pd.read_csv(ctlr_path)

# Read pathway corr results
corr = pd.read_csv(corr_path)

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
    slide.obs['niches'] = pd.read_csv('data/prc/vs/{0}/niches.csv'.format(sample_id), index_col=0)
    slide.obs_names = [sample_id + '|' + i for i in slide.obs_names]
    scores.append(slide)
scores = ad.concat(scores, join='outer')
scores.X[np.isnan(scores.X)] = 0.
scores = scores[(scores.obs['niches'] != 'Ependym') & (scores.obs['niches'] != 'GM')].copy()

# Compute mean scores
mean_scores = dc.get_pseudobulk(
    adata=scores,
    sample_col='niches',
    groups_col=None,
    mode='mean',
    min_cells=0,
    min_counts=0,
)

# Compute mean scores per sample
mean_scores_sample = dc.get_pseudobulk(
    adata=scores,
    sample_col='Sample id',
    groups_col=None,
    mode='mean',
    min_cells=0,
    min_counts=0,
)

# Compute mean scores per sample
mean_scores_niches = dc.get_pseudobulk(
    adata=scores[scores.obs['Lesion type'] != 'Ctrl'].copy(),
    sample_col='Sample id',
    groups_col='niches',
    mode='mean',
    min_cells=0,
    min_counts=0,
)

# Make psbulk of snRNA-seq
pdata = sc.read_h5ad(sn_annot_path)
pdata = dc.get_pseudobulk(
    pdata,
    sample_col='Sample id',
    groups_col='leiden',
    layer='counts',
    mode='sum',
    min_cells=10,
    min_counts=1000
)

cell_types = pdata.obs['leiden'].unique()
adata = []
for cell_type in cell_types:
    cdata = pdata[pdata.obs['leiden'] == cell_type].copy()
    genes = dc.filter_by_expr(cdata, group='Lesion type', min_count=10, min_total_count=15)
    print(cell_type, genes.size)
    
    cdata = cdata[:, genes].copy()
    sc.pp.normalize_total(cdata, target_sum=1e4)
    sc.pp.log1p(cdata)
    adata.append(cdata)
adata = ad.concat(adata, join='outer')

res['n_cs'] = [str(a).count('_') + str(b).count('_') for a, b in zip(res['source_cs'], res['target_cs'])]

intrs = (
    res
    .sort_values(['pvals_adj', 'n_cs'], ascending=[True, False])
    .drop_duplicates(['names', 'lesion_type'], keep='first')
    .groupby('lesion_type')
    .head(10)
)

def plot_inter_heatmap(adata, var_names, ax):
    sc.pl.matrixplot(
        adata=adata,
        var_names=var_names,
        groupby='niches',
        standard_scale='var',
        swap_axes=True,
        ax=ax,
        show=False,
        categories_order=['WM', 'PPWM', 'LR', 'LC', 'VI'],
        colorbar_title='Z-scaled interaction score',
    )

fig1, axes = plt.subplots(3, 1, figsize=(2, 9), dpi=150, sharex=True)
axes = axes.ravel()
plot_inter_heatmap(mean_scores, var_names=intrs[intrs['lesion_type'] == 'Ctrl']['names'].values, ax=axes[0])
plot_inter_heatmap(mean_scores, var_names=intrs[intrs['lesion_type'] == 'CA']['names'].values, ax=axes[1])
plot_inter_heatmap(mean_scores, var_names=intrs[intrs['lesion_type'] == 'CI']['names'].values, ax=axes[2])

def get_expr_df(adata, feature, group):
    vals = adata[:, feature].X.ravel().copy()
    df = pd.DataFrame()
    df['group'] = adata.obs.loc[:, group].values
    df['sample'] = adata.obs.loc[:, 'Sample id'].values
    df['value'] = vals
    return df.reset_index(drop=True)

def read_slide(sample_id):
    slide = sc.read_h5ad('data/prc/vs/{0}/adata.h5ad'.format(sample_id))
    slide.obsm['score'] = pd.read_csv('data/prc/vs/{0}/ctlr_scores.csv'.format(sample_id), index_col=0)
    return dc.get_acts(slide, 'score')

def plot_interaction(adata, scores, inter, sample_id, group):
    source, target, ligand, receptor = inter.split('^')
    lg_df = get_expr_df(adata[adata.obs['leiden'] == source], ligand, group)
    rp_df = get_expr_df(adata[adata.obs['leiden'] == target], receptor, group)
    it_df = get_expr_df(scores, inter, group)
    slide = read_slide(sample_id)
    
    fg1, ax = plt.subplots(1, 1, figsize=(2, 2), tight_layout=True, dpi=150)
    sns.boxplot(data=lg_df, x='group', y='value', hue='group', legend=False, ax=ax, fliersize=0, palette=palette)
    sns.stripplot(data=lg_df, x='group', y='value', hue='group', legend=False, ax=ax, linewidth=1, palette=palette)
    ax.set_title(source)
    ax.set_ylabel(ligand)
    ax.set_xlabel('')

    fg2, ax = plt.subplots(1, 1, figsize=(2, 2), tight_layout=True, dpi=150)
    sns.boxplot(data=rp_df, x='group', y='value', hue='group', legend=False, ax=ax, fliersize=0, palette=palette)
    sns.stripplot(data=rp_df, x='group', y='value', hue='group', legend=False, ax=ax, linewidth=1, palette=palette)
    ax.set_title(target)
    ax.set_ylabel(receptor)
    ax.set_xlabel('')

    fg3, ax = plt.subplots(1, 1, figsize=(2, 2), tight_layout=True, dpi=150)
    sns.boxplot(data=it_df, x='group', y='value', hue='group', legend=False, ax=ax, fliersize=0, palette=palette)
    sns.stripplot(data=it_df, x='group', y='value', hue='group', legend=False, ax=ax, linewidth=1, palette=palette)
    ax.set_title(inter)
    ax.set_ylabel('score')
    ax.set_xlabel('')

    fg4, ax = plt.subplots(1, 1, figsize=(3, 3), tight_layout=True, dpi=150)
    sc.pl.spatial(slide, color=inter, size=1.5, cmap='coolwarm', ax=ax, show=False, frameon=False)

    return [fg1, fg2, fg3, fg4]

def plot_corr(corr, inter):
    ltype = 'CA'
    fg, ax = plt.subplots(1, 1, figsize=(2, 2), dpi=150)
    data = pd.merge(corr[corr['inter'] == inter], meta[['Sample id', 'Lesion type']].rename(columns={'Sample id': 'sample_id'}))
    data = data[data['Lesion type'] == ltype]
    names = data.groupby(['pathway'])['corr'].median().sort_values(ascending=False).reset_index().head(10)['pathway'].values
    data = data[data['pathway'].isin(names)]
    sns.boxplot(data=data, y='pathway', x='corr', color=palette[ltype], fliersize=5, order=names, ax=ax)
    ax.set_ylabel('')
    ax.set_xlabel('Pearson correlation')
    ax.set_xlim(0, 1)
    ax.set_title(inter)
    return fg

inters = ['AS^MG^HMGB1^CD163', 'AS^MG^HMGB1^TLR2', 'MG^EC^CD14^ITGB1', 'MG^AS^CD14^ITGB1']
meta = meta.reset_index()
inters_figs = []
for inter in inters:
    inters_figs.extend(plot_interaction(
        adata=adata,
        scores=mean_scores_sample,
        inter=inter,
        sample_id='MS377T',
        group='Lesion type'
    ))
    inters_figs.append(plot_corr(corr, inter))

# Save to pdf
pdf = matplotlib.backends.backend_pdf.PdfPages(plot_path)
for fig in [fig1] + inters_figs:
    pdf.savefig(fig, bbox_inches='tight')
pdf.close()
