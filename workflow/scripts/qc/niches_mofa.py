import argparse
from composition_stats import closure, ilr, clr
import os
import pandas as pd
import numpy as np
import argparse
import scanpy as sc
import anndata as ad
import decoupler as dc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import scanpy.external as sce


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-m','--meta_path', required=True)
parser.add_argument('-d','--colors_dict', required=True)
parser.add_argument('-p','--plot_path', required=True)
args = vars(parser.parse_args())

meta_path = args['meta_path']
colors_dict = args['colors_dict']
plot_path = args['plot_path']

# Read meta
meta = pd.read_csv(meta_path)
vs_samples = meta[~meta['Batch vs'].isnull()]['Sample id'].values.astype('U')
colors_dict = dict(item.split('_') for item in colors_dict.strip("'").split(';'))

# Process hallmarks and abunds
adata = []
for sample_id in vs_samples:
    slide = sc.read_h5ad('data/prc/vs/{0}/adata.h5ad'.format(sample_id))
    slide.obsm['pathway'] = pd.read_csv('data/prc/vs/{0}/hallmarks.csv'.format(sample_id), index_col=0)
    slide.obsm['abunds'] = pd.read_csv('data/prc/vs/{0}/abunds.csv'.format(sample_id), index_col=0)
    slide.obsm['props'] = pd.read_csv('data/prc/vs/{0}/props.csv'.format(sample_id), index_col=0)
    slide.obs['niches'] = pd.read_csv('data/prc/vs/{0}/niches.csv'.format(sample_id), index_col=0)
    slide.obs_names = [sample_id + '_' + i for i in slide.obs_names]
    slide = slide[~slide.obs['areas'].isnull()].copy()
    adata.append(slide)
adata = ad.concat(adata, join='outer')
adata = adata[np.sort(adata.obs_names), :].copy()

adata.obsm['X_pca'] = adata.obsm['pathway'].values.copy()
sc.pl.pca(adata, color=["Sample id", "Lesion type"])

def stackbar(y, type_names, title, level_names, cmap, ax):

    n_bars, n_types = y.shape

    r = np.array(range(n_bars))
    sample_sums = np.sum(y, axis=1)

    barwidth = 0.85
    cum_bars = np.zeros(n_bars)

    for n in range(n_types):
        bars = [i / j * 100 for i, j in zip([y[k][n] for k in range(n_bars)], sample_sums)]
        ax.bar(r, bars, bottom=cum_bars, color=cmap[n], width=barwidth, label=type_names[n], linewidth=0)
        cum_bars += bars

    ax.set_title(title)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, ncol=2)
    ax.set_xticks(r)
    ax.set_xticklabels(level_names, rotation=90)
    ax.set_ylabel("Proportion")
    ax.grid(False)
    ax.margins(0)

    return ax


# Stacked bars
fig1, axes = plt.subplots(2, 2, figsize=(12, 6), tight_layout=True, dpi=150)
axes = axes.ravel()

# Barplot
ax = axes[0]
df = adata.obs.assign(total=lambda x: 0).groupby(['Sample id', 'niches']).count().reset_index()
df = df.pivot(index='niches', columns='Sample id', values='total')
stackbar(y=df.values, type_names=df.columns, title='Cluster representation', level_names=df.index, cmap=adata.uns['Sample id_colors'], ax=ax)

ax = axes[1]
df = adata.obs.groupby(['niches', 'Lesion type'])[['total_counts']].count().reset_index()
df = df.pivot(index='niches', columns='Lesion type', values='total_counts')
stackbar(y=df.values, type_names=df.columns, title='Cluster representation', level_names=df.index, cmap=adata.uns['Lesion type_colors'], ax=ax)

ax = axes[2]
df = adata.obs.groupby(['niches', 'areas'])[['total_counts']].count().reset_index()
df = df.pivot(index='niches', columns='areas', values='total_counts')
areas_cmap = [colors_dict[i] for i in df.columns]
stackbar(y=df.values, type_names=df.columns, title='Cluster representation', level_names=df.index, cmap=areas_cmap, ax=ax)

ax = axes[3]
df = adata.obs.groupby(['niches', 'Sample id'])[['total_counts']].count().reset_index()
df = df.pivot(index='Sample id', columns='niches', values='total_counts')
niches_cmap = [colors_dict[i] for i in df.columns]
stackbar(y=df.values, type_names=df.columns, title='Cluster representation', level_names=df.index, cmap=niches_cmap, ax=ax)

df = adata.obs[['Sample id', 'Lesion type', 'niches']].copy()
df['niches'] = df['niches'].astype('U').copy()
df = df[(df['niches'] != 'GM') & (df['niches'] != 'Ependym')]

total_counts = df.groupby('Sample id')['niches'].count()

# Calculate the count of each niche for each sample id
niche_counts = df.groupby(['Sample id', 'niches']).size()

# Calculate the proportion of each niche inside each sample
props = niche_counts.div(total_counts, level='Sample id').reset_index(name='proportion')

props = pd.merge(props, meta[['Sample id', 'Lesion type']])
#props.loc[props['Sample id'] == 'MS229', 'Lesion type'] = 'CI'

fig2, ax = plt.subplots(1, 1, figsize=(4, 3))
sns.boxplot(data=props, x='niches', y='proportion', hue='Lesion type', ax=ax, boxprops={'alpha': 0.6})
sns.stripplot(data=props, x='niches', y='proportion', hue='Lesion type', dodge=True, ax=ax, legend=False)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
ax.set_xlabel('Niches')
ax.set_ylabel('Proportions')

pbulk = dc.get_pseudobulk(
    adata=adata,
    sample_col='Sample id',
    groups_col='niches',
    layer='counts'
)

by_expr = dc.filter_by_expr(
    adata=pbulk,
    group=None,
)
by_prop = dc.filter_by_prop(
    pbulk,
    min_prop=0.2,
    min_smpls=2,
)

pbulk = pbulk[:, np.intersect1d(by_expr, by_prop)].copy()

sc.pp.normalize_total(pbulk)
sc.pp.log1p(pbulk)

gex_df = (
    pd.DataFrame(np.corrcoef(pbulk.X), index=pbulk.obs_names, columns=pbulk.obs_names)
    .reset_index()
    .melt(id_vars='index', value_vars=pbulk.obs_names)
    .rename(columns={'index': 'sample_a', 'variable': 'sample_b', 'value': 'corr'})
    .assign(niche_a=lambda x: [i.split('_')[1] for i in x['sample_a']])
    .assign(niche_b=lambda x: [i.split('_')[1] for i in x['sample_b']])
    .assign(id_a=lambda x: [i.split('_')[0] for i in x['sample_a']])
    .assign(id_b=lambda x: [i.split('_')[0] for i in x['sample_b']])
    .assign(is_same=lambda x: x['niche_a'] == x['niche_b'])
)
gex_df = gex_df[gex_df['id_a'] != gex_df['id_b']]
gex_df = pd.merge(gex_df, pbulk.obs.reset_index(names='sample_a')[['sample_a', 'Lesion type']])

props = dc.get_pseudobulk(
    adata=dc.get_acts(adata, 'props'),
    sample_col='Sample id',
    groups_col='niches',
    mode='mean',
    min_cells=0,
    min_counts=0
)
props.X = clr(closure(props.X))

prp_df = (
    pd.DataFrame(np.corrcoef(props.X), index=props.obs_names, columns=props.obs_names)
    .reset_index()
    .melt(id_vars='index', value_vars=props.obs_names)
    .rename(columns={'index': 'sample_a', 'variable': 'sample_b', 'value': 'corr'})
    .assign(niche_a=lambda x: [i.split('_')[1] for i in x['sample_a']])
    .assign(niche_b=lambda x: [i.split('_')[1] for i in x['sample_b']])
    .assign(id_a=lambda x: [i.split('_')[0] for i in x['sample_a']])
    .assign(id_b=lambda x: [i.split('_')[0] for i in x['sample_b']])
    .assign(is_same=lambda x: x['niche_a'] == x['niche_b'])
)
prp_df = prp_df[prp_df['id_a'] != prp_df['id_b']]
prp_df = pd.merge(prp_df, pbulk.obs.reset_index(names='sample_a')[['sample_a', 'Lesion type']])

fig3, axes = plt.subplots(1, 2, figsize=(8, 4), tight_layout=True, sharey=True)
axes = axes.ravel()

ax = axes[0]
sns.boxplot(data=gex_df, x='niche_a', y='corr', hue='is_same', ax=ax)
ax.legend([],[], frameon=False)
ax.set_title('Gene expression')
ax.set_xlabel('')
ax.set_ylabel('Pearson R2')

ax = axes[1]
ax.set_title('Cell type proportions')
sns.boxplot(data=prp_df, x='niche_a', y='corr', hue='is_same', ax=ax)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, title='is_same')
ax.set_xlabel('')
ax.set_ylabel('Pearson R2')

def plot_heatmap(obsm_key, title, ax, scale=False, flip=False, remove=False, cmap='Blues', figsize=(3, 9), cbar=True, square=False):
    # Compute sign
    has_groups = obsm_key.obs['niches'].unique().size > 1
    if has_groups:
        df = dc.rank_sources_groups(obsm_key, groupby='niches', reference='rest', method='wilcoxon')
        df = df[(df['statistic'] > 0) & (df['meanchange'] > 0) & (df['pvals_adj'] < 0.05)]
    else:
        df = pd.DataFrame(columns=['group', 'names'])

    # Compute mean
    obsm_key = dc.get_pseudobulk(
        adata=obsm_key,
        sample_col='niches',
        groups_col=None,
        min_cells=0,
        min_counts=0,
        mode='mean',
        skip_checks=True
    )
    if scale and has_groups:
        sc.pp.scale(obsm_key)

    sign = pd.DataFrame(np.empty(obsm_key.shape, dtype = str), index=obsm_key.obs_names, columns=obsm_key.var_names)
    for r, c in zip(df['group'], df['names']): 
        sign.loc[r, c] = '*'

    if remove:
        msk = np.any(sign.values.astype('U') == '*', axis=0)
        sign = sign.loc[:, msk]
        obsm_key = obsm_key[:, msk].copy()


    if flip:
        obsm_key = obsm_key.T
        sign = sign.T

    sns.heatmap(
        obsm_key.X,
        cmap=cmap,
        robust=True,
        xticklabels=obsm_key.var_names,
        yticklabels=obsm_key.obs_names,
        cbar_kws={"shrink": .4, "aspect": 5, 'label': title},
        annot=sign,
        fmt='',
        annot_kws={'fontweight': 'black', 'color': 'red'},
        square=square,
        ax=ax,
        cbar=cbar,
        vmin=-1,
        vmax=1,
    )

pathway = dc.get_acts(adata, 'pathway')
props = dc.get_acts(adata, 'props')

fig4, axes = plt.subplots(1, 2, figsize=(10, 10), tight_layout=True)
axes = axes.ravel()

ax = axes[0]
_ = plot_heatmap(pathway, title='Mean scaled scores', scale=True, cmap='RdBu_r', ax=ax, flip=True, square=True)

ax = axes[1]
_ = plot_heatmap(props, title='Mean scaled props', scale=True, cmap='RdBu_r', cbar=True, ax=ax, flip=True, square=True)

pdf = matplotlib.backends.backend_pdf.PdfPages(plot_path)
for fig in [fig1, fig2, fig3, fig4]:
    pdf.savefig(fig, bbox_inches='tight')
pdf.close()
