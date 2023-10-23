import anndata as ad
import decoupler as dc
import pandas as pd
import numpy as np
import scanpy as sc
import decoupler as dc
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import seaborn as sns
from tqdm import tqdm
import scipy
import argparse
sc.settings.set_figure_params(dpi=150, frameon=False, figsize=(3, 3))
plt.rcParams['axes.axisbelow'] = True

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-m','--meta_path', required=True)
parser.add_argument('-r','--pres_path', required=True)
parser.add_argument('-o','--out_path', required=True)
parser.add_argument('-p','--plot_path', required=True)
args = vars(parser.parse_args())

meta_path = args['meta_path']
pres_path = args['pres_path']
out_path = args['out_path']
plot_path = args['plot_path']

# Read meta
meta = pd.read_csv(meta_path)
vs_samples = meta[~meta['Batch vs'].isnull()]['Sample id'].values.astype('U')

# Read
adata = []
for sample_id in vs_samples:
    print(sample_id)
    slide = sc.read_h5ad('data/prc/vs/{0}/adata.h5ad'.format(sample_id))
    slide.obs['niches'] = pd.read_csv('data/prc/vs/{0}/niches.csv'.format(sample_id), index_col=0)
    slide.obs['Sample id'] = sample_id
    slide.obs_names = [sample_id + '|' + i for i in slide.obs_names]
    adata.append(slide)
adata = ad.concat(adata, join='outer')

# Get filtered pseudo-bulk profile
pdata = dc.get_pseudobulk(
    adata,
    sample_col='Sample id',
    groups_col='niches',
    layer='counts',
    mode='sum',
    min_cells=10,
    min_counts=1000
)
genes = dc.filter_by_expr(pdata, group=None, min_count=10, min_total_count=15)
pdata = pdata[:, genes].copy()
sc.pp.normalize_total(pdata, target_sum=1e4)
sc.pp.log1p(pdata)

def compute_spearman(adata, groupby, order):
    # Define ordered cats
    order_dict = {k: i for i, k in enumerate(order)}

    # SUbset by cats in order
    msk = np.isin(adata.obs[groupby].astype(str), order)
    
    # Compute x order
    x = np.array([order_dict[k] for k in adata.obs.loc[msk, groupby]])

    # Compute spearman with genes (y)
    df = []
    for name in tqdm(adata.var_names):
        y = adata[msk, name].X.ravel().copy()
        s, p = scipy.stats.spearmanr(x, y)
        df.append([name, s, p])

    # Process dataframe
    df = pd.DataFrame(df, columns=['name', 'stat', 'pval'])
    df['padj'] = dc.p_adjust_fdr(df['pval'])
    df['cat'] = np.sign(df['stat']).astype(int).astype('U')
    return df.sort_values('stat').reset_index(drop=True)

groupby = 'niches'
order = ['WM', 'PPWM', 'LR', 'LC', 'VI']

msk = np.isin(pdata.obs[groupby], order)
pdata = pdata[msk, :].copy()

df = compute_spearman(pdata, groupby=groupby, order=order)
df = df[(df['padj'] < 0.05) & (df['stat'].abs() > 0.5)]

# Write
df.to_csv(out_path, index=False)

# Plots
fig1, ax = plt.subplots(1, 1, figsize=(2, 1.5), tight_layout=True, dpi=150)
df.groupby('cat').count()[['stat']].plot.bar(ax=ax, width=0.9)
ax.get_legend().remove()
ax.set_xlabel('')
ax.set_ylabel('# genes')
ax.set_axisbelow(True)

mdata = dc.get_pseudobulk(
    adata=pdata,
    mode='mean',
    sample_col='niches',
    groups_col=None,
    min_cells=0,
    min_counts=0
)
mdata.obs[groupby] = pd.Categorical(mdata.obs[groupby].astype(str), categories=order)

n_top = 10
genes = np.array(list(df.head(n_top)['name']) + list(df.tail(n_top)['name']))

fig2 = sc.pl.matrixplot(
    mdata[:, genes].copy(),
    genes,
    'niches',
    dendrogram=False,
    colorbar_title='mean z-score',
    standard_scale='var',
    cmap='viridis',
    return_fig=True,
)
fig2.show()
fig2 = fig2.fig

def get_expr_df(adata, gene, groupby, order):
    import seaborn as sns
    gex = adata[:, gene].X.ravel().copy()
    df = pd.DataFrame()
    df['group'] = pd.Categorical(adata.obs.loc[:, groupby].values, categories=order)
    df['sample'] = adata.obs.loc[:, 'Sample id'].values
    df['value'] = gex
    return df.reset_index(drop=True)

gene = 'HDAC11'
gene_df = get_expr_df(pdata, gene, groupby, order)
fig3, ax = plt.subplots(1, 1, figsize=(2.5, 2.5), tight_layout=True, dpi=150)
sns.boxplot(data=gene_df, x='group', y='value', ax=ax, fliersize=0)
sns.stripplot(data=gene_df, x='group', y='value', ax=ax, linewidth=1)
ax.set_ylabel(gene)
ax.set_xlabel('')
ax.tick_params(axis='x', labelrotation=90)

gene = 'LAMA5'
gene_df = get_expr_df(pdata, gene, groupby, order)
fig4, ax = plt.subplots(1, 1, figsize=(2.5, 2.5), tight_layout=True, dpi=150)
sns.boxplot(data=gene_df, x='group', y='value', ax=ax, fliersize=0)
sns.stripplot(data=gene_df, x='group', y='value', ax=ax, linewidth=1)
ax.set_ylabel(gene)
ax.set_xlabel('')
ax.tick_params(axis='x', labelrotation=90)

mat = df.assign(index=lambda x: 'trajectory').pivot(index='index', columns='name', values='stat')
net = pd.read_csv('config/progeny.csv')
net = net.groupby('source').head(100)

pw_acts, pw_pvals = dc.run_ulm(mat=mat, net=net)
fig5 = dc.plot_barplot(
    pw_acts,
    'trajectory',
    top=25,
    vertical=False,
    figsize=(4, 2),
    return_fig=True,
    dpi=150
)
pw_pvals = (
    pw_pvals.T
    .rename(columns={'trajectory': 'pval'})
    .assign(padj=lambda x: dc.p_adjust_fdr(x['pval']))
)
pws = pw_pvals[pw_pvals['padj'] < 0.05].index

fig6, ax = plt.subplots(1, 1, figsize=(3, 3), dpi=150)
dc.plot_targets(df.set_index('name'), stat='stat', source_name='TNFa', net=net, top=5, ax=ax)
ax.set_xlim(-15, 15)
plt.show()

dc.run_ulm(mat=pdata, net=net, use_raw=False)
pwdata = dc.get_acts(pdata, 'ulm_estimate')
pwdata = dc.get_pseudobulk(
    adata=pwdata,
    mode='mean',
    sample_col='niches',
    groups_col=None,
    min_cells=0,
    min_counts=0,
    skip_checks=True
)
pwdata.obs[groupby] = pd.Categorical(pwdata.obs[groupby].astype(str), categories=order)

fig7 = sc.pl.matrixplot(
    pwdata,
    pws,
    'niches',
    colorbar_title='mean activity',
    cmap='coolwarm',
    return_fig=True
)
fig7.show()
fig7 = fig7.fig

net = dc.read_gmt('config/c2.cp.reactome.v2023.1.Hs.symbols.gmt')
net['source'] = [s.split('REACTOME')[1].replace('_', ' ').lstrip() for s in net['source']]
pw = pd.read_csv(pres_path)
pw = pw[pw['pvals'] < 0.15]['source'].unique().astype(str)
net = net[np.isin(net['source'], pw)]

genes = pd.DataFrame(index=df[df['cat'] == '1']['name'])

# Run ora
enr_pvals = dc.get_ora_df(
    df=genes,
    net=net,
    source='source',
    target='target'
)

enr_pvals['-log10(FDR)'] = -np.log10(enr_pvals['FDR p-value'])

fig8, ax = plt.subplots(1, 1, figsize=(3, 3), dpi=150)
dc.plot_barplot_df(
    enr_pvals.sort_values('-log10(FDR)', ascending=False).head(10),
    x='-log10(FDR)',
    y='Term',
    thr=-np.log10(0.05),
    color='#d62728',
    ax=ax
)
ax.set_axisbelow(True)
ax.set_xlim(0, 10)

genes = pd.DataFrame(index=df[df['cat'] == '-1']['name'])

# Run ora
enr_pvals = dc.get_ora_df(
    df=genes,
    net=net,
    source='source',
    target='target'
)

enr_pvals['-log10(FDR)'] = -np.log10(enr_pvals['FDR p-value'])

fig9, ax = plt.subplots(1, 1, figsize=(3, 3), dpi=150)
dc.plot_barplot_df(
    enr_pvals.sort_values('-log10(FDR)', ascending=False).head(10),
    x='-log10(FDR)',
    y='Term',
    thr=-np.log10(0.05),
    color='#1f77b4',
    ax=ax
)
ax.set_axisbelow(True)
ax.set_xlim(0, 10)

# Save to pdf
pdf = matplotlib.backends.backend_pdf.PdfPages(plot_path)
for fig in [fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8, fig9]:
    pdf.savefig(fig, bbox_inches='tight')
pdf.close()
