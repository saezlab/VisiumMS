import pandas as pd
import scanpy as sc
import decoupler as dc
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import argparse
import seaborn as sns
import numpy as np
from numpy.random import default_rng
from tqdm import tqdm


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-o','--h5ad_path', required=True)
parser.add_argument('-a','--ann_path', required=True)
parser.add_argument('-m','--meta_path', required=True)
parser.add_argument('-r','--reactome_path', required=True)
parser.add_argument('-c','--collectri_path', required=True)
parser.add_argument('-p','--plot_path', required=True)
args = vars(parser.parse_args())

h5ad_path = args['h5ad_path']
ann_path = args['ann_path']
meta_path = args['meta_path']
reactome_path = args['reactome_path']
collectri_path = args['collectri_path']
plot_path = args['plot_path']

# Read AS atlas and define AS_C vs the rest
adata = sc.read_h5ad(h5ad_path)
adata.obs['cell_states'] = pd.read_csv(ann_path, index_col=0)
adata.obs['cs'] = ['AS_C' if cs == '11' else 'rest' for cs in adata.obs['leiden']]
adata.obs['cs'] = pd.Categorical(adata.obs.loc[:, 'cs'], categories=['rest', 'AS_C'])

# Remove samples with Ependym
adata = adata[~adata.obs['Sample id'].str.startswith('MS549')].copy()

# Compute marker genes
sc.tl.rank_genes_groups(adata, groupby='leiden', method='t-test_overestim_var')
df = sc.get.rank_genes_groups_df(adata, group='11')
df = df[(df['pvals_adj'] < 0.0001) & (df['logfoldchanges'] > 1)]
df['group'] = 'AS_C'

def get_expr_df(adata, gene):
    from scipy.sparse import issparse
    if issparse(adata.X):
        gex = adata[:, gene].copy().X.A.ravel()
    else:
        gex = adata[:, gene].copy().X.ravel()
    df = pd.DataFrame()
    df['group'] = adata.obs['cs']
    df['sample'] = adata.obs_names
    df['value'] = gex
    return df.reset_index(drop=True)


def plot_violin(adata, gene, ax, palette):
    df = get_expr_df(adata, gene)
    sns.violinplot(data=df, x='group', y='value', hue='group', ax=ax, palette=palette)
    ax.set_ylabel(gene)
    ax.set_xlabel('')
    ax.get_legend().remove()

# Plot top marker genes
fig1, axes = plt.subplots(2, 3, figsize=(5, 3), tight_layout=True, dpi=150)
axes = axes.ravel()
glist = list(df['names'].values[:6].astype('U'))
palette={'rest': 'tab:gray', 'AS_C': 'tab:green'}
for i in range(len(glist)):
    plot_violin(adata=adata, gene=glist[i], ax=axes[i], palette=palette)

def read_slide(sample_id):
    slide = sc.read_h5ad('data/prc/vs/{0}/adata.h5ad'.format(sample_id))
    slide.obs['niches'] = pd.read_csv('data/prc/vs/{0}/niches.csv'.format(sample_id), index_col=0)
    slide.obsm['reactome'] = pd.read_csv('data/prc/vs/{0}/reactome.csv'.format(sample_id), index_col=0)
    slide.obsm['props'] = pd.read_csv('data/prc/vs/{0}/props.csv'.format(sample_id), index_col=0)
    return slide


def get_asc_df(sample_id, markers):
    slide = read_slide(sample_id)
    dc.run_ulm(slide, net=markers, source='group', target='names', weight=None, use_raw=False)
    a_msk = ((slide.obsm['ulm_pvals'] < 0.05) & (slide.obsm['ulm_estimate'] > 0)).values.ravel()
    c_msk = (slide.obsm['props']['AS'] > 0.25).values.ravel()
    msk = a_msk * c_msk
    obs = slide.obs.loc[:, ['Sample id', 'Lesion type', 'niches']].copy()
    obs['has_ASC'] = msk
    return obs


# Read meta
meta = pd.read_csv(meta_path)

# Read slides
vs_samples = meta[~meta['Batch vs'].isnull()]['Sample id'].values.astype('U')

# Remove epedym slides
vs_samples = vs_samples[~np.isin(vs_samples, ['MS549H', 'MS549T'])]
obs = []
for vs_sample in vs_samples:
    print(vs_sample)
    obs.append(get_asc_df(vs_sample, markers=df))
obs = pd.concat(obs)

def test_cs(data, times=1000, seed=42):
    data = data.copy()
    rng = default_rng(seed=seed)
    estimate = data.groupby(['niches'])['has_ASC'].sum()
    counts = pd.Series(0, index=estimate.index)
    idx = np.arange(obs.shape[0])
    for i in tqdm(range(times)):
        np.random.shuffle(idx)
        data['has_ASC'] = data['has_ASC'].iloc[idx].values
        counts += data.groupby(['niches'])['has_ASC'].sum() > estimate
    res = pd.DataFrame(estimate)
    res['pval'] = counts / times
    return res

# Plot abundance of ASC in niches
res = test_cs(obs, times=1000, seed=42)
fig2, ax = plt.subplots(1, 1, figsize=(2, 2), dpi=150)
res.sort_values('has_ASC')['has_ASC'].plot.barh(ax=ax)
ax.set_xlabel('# of AS_C')
ax.set_ylabel('')
ax.set_xlim(0, 100)
plt.show()

def plot_where_asc(sample_id, markers):
    colors_dict = 'WM_tab:blue;GM_tab:orange;PPWM_tab:green;LR_tab:red;LC_tab:purple;Ependym_tab:brown;VI_tab:pink;Claustrum_tab:gray'
    colors_dict = dict(item.split('_') for item in colors_dict.strip("'").split(';'))
    slide = read_slide(sample_id)
    dc.run_ulm(slide, net=markers, source='group', target='names', weight=None, use_raw=False)
    act = dc.get_acts(slide, 'ulm_estimate')
    c_msk = (slide.obsm['props']['AS'] > 0.25).values.ravel()
    a_msk = ((slide.obsm['ulm_pvals'] < 0.05) & (slide.obsm['ulm_estimate'] > 0)).values.ravel()
    act.obs['asc'] = np.where(a_msk, 'asc', 'as')
    act.obs.loc[~c_msk, 'asc'] = 'nas'
    fig, axes = plt.subplots(1, 2, figsize=(3, 6), dpi=300)
    na_color = 'slategray'
    palette={
        'asc': '#e64040',
        'as': '#73b360',
        'nas': 'white' #slategray
    }
    sc.pl.spatial(act, color='asc', size=1.75, cmap='summer', ax=axes[0], colorbar_loc=None,
                  show=False, frameon=False, title='', na_color=na_color, palette=palette)
    sc.pl.spatial(act, color='niches', size=1.75, ax=axes[1], show=False, frameon=False, palette=colors_dict,
                  legend_loc=None, title='', na_color=na_color)
    return fig

fig3 = plot_where_asc('MS377T', df)
fig4 = plot_where_asc('MS411', df)

# Read REACTOME
gmt = dc.read_gmt(reactome_path)
gmt['source'] = [s.split('REACTOME')[1].replace('_', ' ').lstrip() for s in gmt['source']]
msk = ~gmt['source'].str.contains('FETAL|INFECTION|SARS', case=False)
gmt = gmt[msk]
gmt = gmt.drop_duplicates(['source', 'target'])

er = dc.get_ora_df(
    df=df.set_index('names'),
    net=gmt,
).sort_values(['FDR p-value', 'p-value'])
er['-log10(FDR)'] = -np.log10(er['FDR p-value'])

fig5 = dc.plot_dotplot(
    er.head(10),
    x='Combined score',
    y='Term',
    s='Odds ratio',
    c='FDR p-value',
    scale=0.15,
    figsize=(2, 4),
    return_fig=True,
    dpi=150
)

clt = pd.read_csv(collectri_path)
dc.run_ulm(adata, net=clt, use_raw=False)
act = dc.get_acts(adata, 'ulm_estimate')
tfs = dc.rank_sources_groups(
    adata=act,
    groupby='cs'
)
fig6, axes = plt.subplots(1, 3, figsize=(6, 2), tight_layout=True, dpi=150)
glist = ['FOXJ1', 'REST', 'TBX1']
for i in range(len(glist)):
    plot_violin(adata=act, gene=glist[i], ax=axes[i], palette=None)

def plot_where_func(sample_id, markers, net, name, is_tf=True):
    slide = read_slide(sample_id)
    dc.run_ulm(slide, net=markers, source='group', target='names', weight=None, use_raw=False)
    act = dc.get_acts(slide, 'ulm_estimate')
    c_msk = (slide.obsm['props']['AS'] > 0.25).values.ravel()
    a_msk = ((slide.obsm['ulm_pvals'] < 0.05) & (slide.obsm['ulm_estimate'] > 0)).values.ravel()
    if is_tf:
        # tfs
        if 'weight' in net.columns:
            weight = 'weight'
        else:
            weight = None
        dc.run_ulm(slide, net=net, weight=weight, use_raw=False)
        tf = dc.get_acts(slide, 'ulm_estimate')
        tf[~(c_msk * a_msk), name].X = np.nan
        fig, ax = plt.subplots(1, 1, figsize=(2, 2), dpi=300)
        sc.pl.spatial(tf, color=name, size=1.75, cmap='RdBu_r', ax=ax, colorbar_loc='right',
                      show=False, frameon=False, title=name, na_color='slategray', vcenter=0)
    else:
        slide[~(c_msk * a_msk), name].X = np.nan
        fig, ax = plt.subplots(1, 1, figsize=(2, 2), dpi=300)
        sc.pl.spatial(slide, color=name, size=1.75, cmap='viridis', ax=ax, colorbar_loc='right',
                      show=False, frameon=False, title=name, na_color='slategray')
    return fig

fig7 = plot_where_func('MS411', df, clt, name='FOXJ1')
fig8 = plot_where_func('MS411', df, gmt, name='CILIUM ASSEMBLY')
fig9 = plot_where_func('MS411', df, clt, name='FOXJ1', is_tf=False)
fig10 = plot_where_func('MS411', df, clt, name='CETN2', is_tf=False)
fig11 = plot_where_func('MS411', df, clt, name='SERPINA3', is_tf=False)

# Plot net
obs = dc.get_pseudobulk(
    adata=adata,
    sample_col='cs',
    groups_col=None,
    min_cells=0,
    min_counts=0,
    mode='mean'
)['AS_C'].to_df()
act, _ = dc.run_ulm(obs, net=clt)
fig12 = dc.plot_network(
    obs=obs,
    act=act,
    net=clt,
    n_sources=['FOXJ1'],
    n_targets=50,
    node_size=1.,
    figsize=(3, 4),
    c_pos_w='#9b2531',
    c_neg_w='#233c55',
    s_cmap='white',
    t_cmap='white',
    return_fig=True,
    vcenter=True,
    dpi=300
)

# Save to pdf
pdf = matplotlib.backends.backend_pdf.PdfPages(plot_path)
for fig in [fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8, fig9, fig10, fig11, fig12]:
    pdf.savefig(fig, bbox_inches='tight')
pdf.close()
