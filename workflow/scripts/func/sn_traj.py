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
from scipy.cluster import hierarchy
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--ann_path', required=True)
parser.add_argument('-n','--pw_path', required=True)
parser.add_argument('-r','--pres_path', required=True)
parser.add_argument('-o','--out_path', required=True)
parser.add_argument('-p','--plot_path', required=True)
args = vars(parser.parse_args())

ann_path = args['ann_path']
pw_path = args['pw_path']
pres_path = args['pres_path']
out_path = args['out_path']
plot_path = args['plot_path']

# Inputs
groupby = 'Lesion type'
order = ['Ctrl', 'CA', 'CI']

# Get psbulk data
sn_adata = sc.read_h5ad(ann_path)
sn_adata = dc.get_pseudobulk(
    adata=sn_adata,
    sample_col='Sample id',
    groups_col='leiden',
    layer='counts',
).copy()

ctypes = sn_adata.obs['leiden'].unique().astype(str)


def compute_spearman(adata, groupby, order):
    # Define ordered cats
    order_dict = {k: i for i, k in enumerate(order)}

    # Subset by cats in order
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


df = []
for ctype in ctypes:
    print(ctype)
    tmp = sn_adata[sn_adata.obs['leiden'] == ctype].copy()
    genes = dc.filter_by_expr(tmp, group=None, min_count=10, min_total_count=15)
    tmp = tmp[:, genes].copy()
    sc.pp.normalize_total(tmp, target_sum=1e4)
    sc.pp.log1p(tmp)

    cdf = compute_spearman(tmp, groupby=groupby, order=order)
    cdf = cdf[(cdf['padj'] < 0.05) & (cdf['stat'].abs() > 0.5)]
    cdf['ctype'] = ctype
    df.append(cdf)
    
df = pd.concat(df)
df.to_csv(out_path, index=False)

# Plots traj
bar = df.groupby(['ctype', 'cat'])[['name']].count().reset_index()
order = bar.groupby('ctype').sum(numeric_only=True)['name'].sort_values(ascending=False).index


fig1, ax = plt.subplots(1, 1, figsize=(4, 2), tight_layout=True, dpi=150)
(
    bar
    .pivot(index='ctype', columns='cat', values='name')
    .loc[order]
    .plot(kind='bar', stacked=False, ax=ax, width=0.9, color=['#1f77b4', '#d62728'], log=False)
)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
ax.set_xlabel('')
ax.set_ylabel('# genes')
ax.set_axisbelow(True)
ax.grid()

bar = (
    df
    .groupby(['name', 'cat'])
    .count()
    .sort_values('stat', ascending=False)
    .head(20)
    .reset_index()
)
rows = bar['name'].values

fig2, ax = plt.subplots(1, 1, figsize=(3.5, 4), tight_layout=True, dpi=150)
(
    bar
    .pivot(index='name', columns='cat', values='stat')
    .fillna(0)
    .loc[np.flip(rows)]
    .plot(kind='barh', width=0.9, stacked=True, ax=ax, color=['#1f77b4', '#d62728'])
    
)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
ax.set_ylabel('')
ax.set_xlabel('# cell types')
ax.set_axisbelow(True)
ax.grid()

table = (
    df
    .pivot(index='name', columns='ctype', values='stat')
    .loc[rows]
    .dropna(how='all', axis=1)
    .fillna(0)
)

fig3, ax = plt.subplots(1, 1, figsize=(3, 4), tight_layout=True, dpi=150)
sns.heatmap(
    table,
    xticklabels=True,
    yticklabels=True,
    cmap='coolwarm',
    ax=ax,
    cbar_kws={"shrink": 0.5},
    robust=True
)
ax.tick_params(axis='x', rotation=90)
ax.set_ylabel('')
ax.set_xlabel('')

# Run enrichment
def run_ora(df, net):
    ctypes = df['ctype'].unique()
    cats = df['cat'].unique()
    enr = []
    for ctype in ctypes:
        for cat in cats:
            genes = pd.DataFrame(index=df[(df['cat'] == cat) & (df['ctype'] == ctype)]['name'])
            # Run ora
            enr_pvals = dc.get_ora_df(
                df=genes,
                net=net,
                source='source',
                target='target'
            )
            enr_pvals['-log10(FDR)'] = -np.log10(enr_pvals['FDR p-value'])
            enr_pvals = enr_pvals[enr_pvals['FDR p-value'] < 0.05]
            enr_pvals = enr_pvals.sort_values('-log10(FDR)', ascending=False)
            if cat == '-1':
                enr_pvals['-log10(FDR)'] = -1 * enr_pvals['-log10(FDR)']
            enr_pvals['ctype'] = ctype
            enr.append(enr_pvals)
    enr = pd.concat(enr)
    return enr

net = dc.read_gmt(pw_path)
net['source'] = [s.split('REACTOME')[1].replace('_', ' ').lstrip() for s in net['source']]
pw = pd.read_csv(pres_path)
pw = pw[pw['pvals'] < 0.15]['source'].unique().astype(str)
net = net[np.isin(net['source'], pw)]
enr = run_ora(df, net=net)

# Plot enrichment
def find_max_abs(x):
    m = np.max(np.abs(x))
    sign = x[np.abs(x) == m]
    return m * sign

table = (
    enr
    .groupby(['Term', 'ctype'])
    [['-log10(FDR)']]
    .agg(lambda x: find_max_abs(x))
    .reset_index()
    .pivot(index='Term', columns='ctype', values='-log10(FDR)')
)
table['nonnans'] = np.sum(~np.isnan(table.values), axis=1)
table = table[table['nonnans'] > 1].drop(columns='nonnans')
table = table.dropna(how='all', axis=1).fillna(0)
table.loc[:, :] = np.sign(table.values)

rows = hierarchy.dendrogram(hierarchy.linkage(table), no_plot=True)['leaves']
cols = hierarchy.dendrogram(hierarchy.linkage(table.T), no_plot=True)['leaves']

fig4, ax = plt.subplots(1, 1, figsize=(1, 6), dpi=150)
sns.heatmap(
    table.iloc[rows, cols],
    xticklabels=True,
    yticklabels=True,
    cmap='coolwarm',
    cbar_kws={"shrink": 0.5}
)
ax.set_ylabel('')
ax.set_xlabel('')

# Save to pdf
pdf = matplotlib.backends.backend_pdf.PdfPages(plot_path)
for fig in [fig1, fig2, fig3, fig4]:
    pdf.savefig(fig, bbox_inches='tight')
pdf.close()
