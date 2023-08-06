from mofapy2.run.entry_point import entry_point
import mofax as mfx
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

def plot_heatmap(obsm_key, title, ax, scale=False, flip=False, remove=False, cmap='Blues', figsize=(3, 9), cbar=True, square=False):
    # Compute sign
    has_groups = obsm_key.obs['leiden'].unique().size > 1
    if has_groups:
        df = dc.rank_sources_groups(obsm_key, groupby='leiden', reference='rest', method='wilcoxon')
        df = df[(df['statistic'] > 0) & (df['meanchange'] > 0) & (df['pvals_adj'] < 0.05)]
    else:
        df = pd.DataFrame(columns=['group', 'names'])

    # Compute mean
    obsm_key = dc.get_pseudobulk(
        adata=obsm_key,
        sample_col='leiden',
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
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, ncol=1)
    ax.set_xticks(r)
    ax.set_xticklabels(level_names, rotation=90)
    ax.set_ylabel("Proportion")
    ax.grid(False)
    ax.margins(0)

    return ax

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-s','--slide_path', required=True)
parser.add_argument('-n','--n_hvg', required=True)
parser.add_argument('-m','--mofa_path', required=True)
parser.add_argument('-r','--res', required=True)
parser.add_argument('-c','--colors_dict', required=True)
parser.add_argument('-a','--annotation', required=True)
parser.add_argument('-p','--plot_path', required=True)
parser.add_argument('-o','--out_path', required=True)

args = vars(parser.parse_args())

slide_path = args['slide_path']
n_hvg = int(args['n_hvg'])
mofa_path = args['mofa_path']
res = float(args['res'])
colors_dict = args['colors_dict']
annotation = args['annotation']
plot_path = args['plot_path']
out_path = args['out_path']

# Read adata
colors_dict = dict(item.split('_') for item in colors_dict.strip("'").split(';'))
sample_id = os.path.normpath(slide_path).split(os.path.sep)[-2]
print(sample_id)
adata = sc.read_h5ad(slide_path)
adata.uns['areas_colors'] = [colors_dict[c] for c in np.unique(adata.obs['areas'].dropna().values.astype('U'))]
del adata.uns['log1p']
adata.obsm['pathway'] = pd.read_csv('data/prc/vs/{0}/pathway.csv'.format(sample_id), index_col=0)
adata.obsm['props'] = pd.read_csv('data/prc/vs/{0}/props.csv'.format(sample_id), index_col=0)
adata.obsm['clr'] = pd.DataFrame(clr(closure(adata.obsm['props'].values)), index=adata.obs_names, columns=adata.obsm['props'].columns)

# Find HVG
sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg)
adata = adata[:, adata.var['highly_variable']].copy()

# Compute pca
sc.pp.scale(adata)
sc.tl.pca(adata, svd_solver='arpack', use_highly_variable=True)

# Generate input
clr_df = (
        adata.obsm['clr']
        .reset_index(names='sample')
        .melt(id_vars='sample', var_name="feature", value_name="value")
        .assign(view=lambda x: 'clr')
        [["sample", "feature", "value", "view"]]
    )

pca_df = (
        pd.DataFrame(adata.obsm['X_pca'], index=adata.obs_names, columns=['PC{0}'.format(i+1) for i in range(adata.obsm['X_pca'].shape[1])])
        .reset_index(names='sample')
        .melt(id_vars='sample', var_name="feature", value_name="value")
        .assign(view=lambda x: 'pca')
        [["sample", "feature", "value", "view"]]
    )

df = pd.concat([clr_df, pca_df])

# MOFA+
ent = entry_point()
ent.set_data_options(
    scale_views = False,
    center_groups = False
)
ent.set_data_df(df, likelihoods = ['gaussian', 'gaussian'])
ent.set_model_options(
    factors = 15,
    spikeslab_weights = False,
    ard_weights = False,
    ard_factors = False
)
ent.set_train_options(
    convergence_mode = "fast",
    gpu_mode = False,
    seed = 1
)
ent.build()
ent.run()

# Write and read model
ent.save(outfile=mofa_path)
model = mfx.mofa_model(mofa_path)

# Run leiden
adata.obsm['X_mofa'] = pd.DataFrame(model.get_factors(), index=model.get_cells()['cell'].values)
adata.obsm['X_umap'] = adata.obsm['X_mofa'].values
sc.pp.neighbors(adata, use_rep="X_mofa")
sc.tl.leiden(adata, resolution=res)

# Process annotation and color
sc.pl.umap(adata=adata, color='leiden')
if annotation != 'None':
    annotation = {key: value for key, value in (pair.split(":") for pair in annotation.split(";"))}
    # Remove cells that are not in annot
    msk = np.isin(adata.obs['leiden'].values.astype('U'), np.array(list(annotation.keys()), dtype='U'))
    adata = adata[msk].copy()
    # Rename obs
    adata.obs['leiden'] = [annotation[c] for c in adata.obs['leiden']]
    leiden_dict = colors_dict.copy()
else:
    leiden_dict = {k:v for k,v in zip(adata.obs['leiden'].unique(), adata.uns['leiden_colors'])}

# Extract obsms
props = dc.get_acts(adata, 'props')
pathway = dc.get_acts(adata, 'pathway')

figsize=(12, 12)
fig1, axes = plt.subplots(3, 3, figsize=figsize, tight_layout=True)
fig1.suptitle(sample_id)
axes = axes.ravel()

# R2
ax = axes[0]
r2 = model.get_r2()
r2['Factor'] = [r.replace('actor', '') for r in r2['Factor']]
sns.barplot(data=r2, x='Factor', y='R2', hue='View', ax=ax)
ax.legend(frameon=False)
ax.tick_params(axis='x', labelrotation=90)
ax.set_xlabel('')

# Plot mofa components
ax = axes[1]
sc.pl.umap(adata=adata, color='areas', return_fig=False, show=False, ax=ax, palette=colors_dict)
ax.set_xlabel('F1')
ax.set_ylabel('F2')

ax = axes[2]
sc.pl.umap(adata=adata, color='leiden', return_fig=False, show=False, ax=ax, palette=leiden_dict)
ax.set_xlabel('F1')
ax.set_ylabel('F2')

ax = axes[3]
df = adata.obs.groupby(['leiden', 'areas'])[['total_counts']].count().reset_index()
df = df.pivot(index='leiden', columns='areas', values='total_counts')
stackbar(y=df.values, type_names=df.columns, title='', level_names=df.index, cmap=adata.uns['areas_colors'], ax=ax)

ax = axes[4]
sc.pl.spatial(adata=adata, color='areas', return_fig=False, show=False, ax=ax, size=1.5, frameon=False, palette=colors_dict)

ax = axes[5]
sc.pl.spatial(adata=adata, color='leiden', return_fig=False, show=False, ax=ax, size=1.5, frameon=False, palette=leiden_dict)

ax = axes[6]
sc.pl.spatial(adata=props, color='OL', return_fig=False, show=False, ax=ax, size=1.5, frameon=False)

ax = axes[7]
sc.pl.spatial(adata=props, color='MG', return_fig=False, show=False, ax=ax, size=1.5, frameon=False)

ax = axes[8]
sc.pl.spatial(adata=props, color='AS', return_fig=False, show=False, ax=ax, size=1.5, frameon=False)


fig2, axes = plt.subplots(1, 2, figsize=figsize, tight_layout=True)
fig2.suptitle(sample_id)
axes = axes.ravel()

ax = axes[0]
_ = plot_heatmap(pathway, title='Mean scaled scores', scale=True, cmap='RdBu_r', ax=ax, flip=True, square=True)

ax = axes[1]
_ = plot_heatmap(props, title='Mean scaled props', scale=True, cmap='RdBu_r', cbar=True, ax=ax, flip=True, square=True)

pdf = matplotlib.backends.backend_pdf.PdfPages(plot_path)
for fig in [fig1, fig2]:
    pdf.savefig(fig)
pdf.close()

# Save obs
adata.obs[['leiden']].to_csv(out_path)
