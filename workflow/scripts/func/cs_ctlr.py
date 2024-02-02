import pandas as pd
import numpy as np
from tqdm import tqdm
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
parser.add_argument('-g','--markers_path', required=True)
parser.add_argument('-t','--table_path', required=True)
parser.add_argument('-s','--snlr_path', required=True)
parser.add_argument('-c','--ctlr_path', required=True)
parser.add_argument('-m','--meta_path', required=True)
parser.add_argument('-o','--out_path', required=True)
parser.add_argument('-q','--colors_dict', required=True)
parser.add_argument('-p','--plot_path', required=True)
args = vars(parser.parse_args())

markers_path = args['markers_path']
table_path = args['table_path']
snlr_path = args['snlr_path']
ctlr_path = args['ctlr_path']
meta_path = args['meta_path']
out_path = args['out_path']
colors_dict = args['colors_dict']
plot_path = args['plot_path']

# Get palette
palette = dict(item.split(':') for item in colors_dict.strip("'").split(';'))

# Read markers
ctypes = pd.read_csv(markers_path)['cell_type'].unique().astype('U')
deg = []
for ctype in ctypes:
    tmp = pd.read_csv('data/prc/ctypes/{0}_deg.csv'.format(ctype))
    if tmp.shape[0] > 0:
        deg.append(tmp)
deg = pd.concat(deg).set_index(['group'])
tc_bc_sc = deg.index[deg.index.str.startswith('BC') | deg.index.str.startswith('TC') | deg.index.str.startswith('SC')].unique().values.astype(str)

# Read composition changes of cell states
df = pd.read_csv(table_path)
df = df.rename(columns={'ctype': 'cell_state'})
df = df[(df['padj'] < 0.10) & (df['type'] == 'cs')]
df['lesion_type'] = [g if np.sign(s) > 0 else r for g, r, s in zip(df['group'], df['ref'], df['stat'])]
df = df[['lesion_type', 'cell_state']].drop_duplicates().sort_values(['lesion_type', 'cell_state'])
df['ctype'] = [x.split('_')[0] for x in df['cell_state']]

# Add TC, BC and SC cell states to CA and CI since no presence in Ctrl
rows = []
for cs in tc_bc_sc:
    ct = cs.split('_')[0]
    rows.append(['CA', ct, cs])
    rows.append(['CI', ct, cs])
rows = pd.DataFrame(rows, columns=['lesion_type', 'ctype', 'cell_state'])
df = pd.concat([df, rows]).sort_values(['lesion_type', 'ctype', 'cell_state']).drop_duplicates(['lesion_type', 'cell_state'])

# Remove cell states with no markers
df = df[np.isin(df['cell_state'], deg.index.unique())]
df = df.set_index(['lesion_type', 'ctype'])

sn_cc = pd.read_csv(snlr_path)
vs_cc = pd.read_csv(ctlr_path)
vs_cc[['source', 'target', 'ligand', 'receptor']] = vs_cc['names'].str.split('^', expand=True, regex=False)[[0, 1, 2, 3]]
vs_cc['lesion_type'] = [x.split('vs')[0] if s > 0 else x.split('vs')[1] for x, s in zip(vs_cc['contrast'], vs_cc['sign'])]
sn_cc['lesion_type'] = [x.split('vs')[0] if s > 0 else x.split('vs')[1] for x, s in zip(sn_cc['contrast'], sn_cc['sign'])]

# Check which interactions have cstates
col_source_cs = []
col_target_cs = []
col_hascs = []
iter = zip(vs_cc['lesion_type'], vs_cc['ligand'], vs_cc['receptor'], vs_cc['source'], vs_cc['target'])

for lesion_type, ligand, receptor, source, target in tqdm(iter):
    
    source_cs = df.loc[lesion_type]
    s_msk = False
    source_list = []
    if source in source_cs.index:
        source_cs = source_cs.loc[[source]]['cell_state'].values.astype(str)
        if np.any(np.isin(source_cs, deg.index)):
            s_msk = deg.loc[source_cs].reset_index().groupby('group')['names'].apply(lambda x: np.any(x.str.contains(ligand)))
            source_list = source_cs[s_msk]
    col_source_cs.append(';'.join(source_list))
    
    target_cs = df.loc[lesion_type]
    t_msk = False
    target_list = []
    if target in target_cs.index:
        target_cs = target_cs.loc[[target]]['cell_state'].values.astype(str)
        if np.any(np.isin(target_cs, deg.index)):
            t_msk = deg.loc[target_cs].reset_index().groupby('group')['names'].apply(lambda x: np.any(x.str.contains(receptor)))
            target_list = target_cs[t_msk]
    col_target_cs.append(';'.join(target_list))

    hascs = False
    if np.any(s_msk) or np.any(t_msk):
        hascs = True
    col_hascs.append(hascs)
vs_cc['source_cs'] = col_source_cs
vs_cc['target_cs'] = col_target_cs
vs_cc['hascs'] = col_hascs
res = vs_cc[vs_cc['hascs']].drop(columns='hascs')
res.to_csv(out_path, index=False)

# Format for plotting
bar = pd.concat([
    sn_cc[['contrast', 'interaction']].assign(type='sn'),
    vs_cc[['contrast', 'names']].rename(columns={'names': 'interaction'}).assign(type='vs'),
    res.assign(interaction=lambda x: x['ligand'] + '^' + x['receptor'])[['contrast', 'interaction']].assign(type='cs')
])
bar = bar.groupby(['contrast', 'type']).count().reset_index()
bar['contrast'] = pd.Categorical(bar['contrast'].values.astype(str), categories=['CAvsCtrl', 'CIvsCtrl', 'CAvsCI'])

fig1, axes = plt.subplots(1, 2, figsize=(4.5, 2.5), tight_layout=True, dpi=150, sharey=True)
ax = axes[0]
(
    bar
    .pivot(index='contrast', columns='type', values='interaction')
    [['cs', 'vs', 'sn']]
    .plot(kind='bar', stacked=True, ax=ax, width=0.9, log=False)
)
ax.get_legend().remove()
ax.set_xlabel('')
ax.set_ylabel('# interactions')
ax.set_yscale('log')
ax.set_yticks([1, 10, 100, 1000])
ax.grid(axis='y')
ax.set_axisbelow(True)

bar = pd.concat([
    sn_cc[['lesion_type', 'interaction']].assign(type='sn'),
    vs_cc[['lesion_type', 'names']].rename(columns={'names': 'interaction'}).assign(type='vs'),
    res.assign(interaction=lambda x: x['ligand'] + '^' + x['receptor'])[['lesion_type', 'interaction']].assign(type='cs')
])
bar = bar.groupby(['lesion_type', 'type']).count().reset_index()
bar['lesion_type'] = pd.Categorical(bar['lesion_type'].values.astype(str), categories=['Ctrl', 'CA', 'CI'])
ax = axes[1]
(
    bar
    .pivot(index='lesion_type', columns='type', values='interaction')
    [['cs', 'vs', 'sn']]
    .plot(kind='bar', stacked=True, ax=ax, width=0.9, log=False)
)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
ax.set_xlabel('')
ax.set_ylabel('# interactions')
ax.set_yscale('log')
ax.set_yticks([1, 10, 100, 1000])
ax.grid(axis='y')
ax.set_axisbelow(True)

bar = (
    res[['lesion_type', 'names', 'source', 'target']]
    .melt(['lesion_type', 'names'])
    .drop_duplicates(['lesion_type', 'names', 'value'])
    .groupby(['lesion_type', 'value'])
    .count()
    .reset_index()
    [['lesion_type', 'value', 'variable']]
)
bar['lesion_type'] = pd.Categorical(bar['lesion_type'].values.astype(str), categories=['Ctrl', 'CA', 'CI'])
order = bar.groupby(['value']).sum(numeric_only=True).sort_values('variable', ascending=False).index
bar['value'] = pd.Categorical(bar['value'].values.astype(str), categories=order)

fig2, ax = plt.subplots(1, 1, figsize=(4, 2), tight_layout=True, dpi=150)
(
    bar
    .pivot(index='value', columns='lesion_type', values='variable')
    .plot(kind='bar', ax=ax, width=0.9, stacked=True, color=palette)
)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
ax.set_xlabel('')
ax.set_ylabel('# interactions')
ax.grid(axis='y')
ax.set_axisbelow(True)

bar = (
    res[['lesion_type', 'ligand', 'receptor']]
    .melt(['lesion_type'])
    .groupby(['lesion_type', 'value'])
    .count()
    .reset_index()
    .sort_values('variable')
)
bar['lesion_type'] = pd.Categorical(bar['lesion_type'].values.astype(str), categories=['Ctrl', 'CA', 'CI'])
order = bar.groupby(['value']).sum(numeric_only=True).sort_values('variable', ascending=True).index

fig3, ax = plt.subplots(1, 1, figsize=(4, 4), tight_layout=True, dpi=150)
(
    bar
    .pivot(index='value', columns='lesion_type', values='variable')
    .loc[order[-20:]]
    .plot(kind='barh', ax=ax, width=0.9, stacked=True, color=palette)
)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
ax.set_xlabel('')
ax.set_ylabel('# interactions')
ax.grid(axis='x')
ax.set_axisbelow(True)

g = res[['source', 'target', 'names', 'lesion_type']].groupby(['source', 'target', 'lesion_type']).count().reset_index().rename(columns={'names': 'n'})
g['edge'] = ['^'.join(sorted([x, y])) for x, y in zip(g['source'], g['target'])]
g = g.groupby(['edge', 'lesion_type']).sum(numeric_only=True).reset_index()
g[['source', 'target']] = g['edge'].str.split('^', expand=True, regex=False)[[0, 1]]
g = g[['source', 'target', 'lesion_type', 'n']]
g.loc[g['lesion_type'] == 'CA', 'n'] = -(g.loc[g['lesion_type'] == 'CA', 'n'])
g.loc[g['lesion_type'] == 'CI', 'n'] = -(g.loc[g['lesion_type'] == 'CI', 'n'])
g_df = g.groupby(['source', 'target']).sum(numeric_only=True).reset_index()
g_df = g_df[abs(g_df['n']) > 1]

cmap = plt.get_cmap('coolwarm')

g = nx.from_pandas_edgelist(g_df, source='source', target='target', edge_attr=['n']) # create_using=nx.DiGraph

c_dict = {'CA': 'red', 'Ctrl': 'green'}
edges = g.edges()
colors = ['red' if g[u][v]['n'] < 0 else 'green' for u,v in edges]
weights = [abs(g[u][v]['n']) for u,v in edges]
weights = np.log(np.array(weights) + 1) + 1

layout = {
    'TC': np.array([1.00000000e+00, 1.98682148e-08]),
    'BC': np.array([0.76604444, 0.64278759]),
    'EC': np.array([0.17364823, 0.98480774]),
    'AS': np.array([-0.50000005,  0.8660254 ]),
    'OL': np.array([-0.9396926 ,  0.34202023]),
    'OPC': np.array([-0.9396926 , -0.34202013]),
    'MG': np.array([-0.4999999 , -0.86602542]),
    'NEU': np.array([ 0.17364817, -0.9848077 ]),
    'SC': np.array([ 0.76604432, -0.64278773])
}

fig4, ax = plt.subplots(1, 1, figsize=(3, 3), tight_layout=True, dpi=150)
nx.draw(
    g,
    pos=layout,
    with_labels=True,
    edge_color=colors,
    width=weights,
    ax=ax,
    alpha=1,
    node_size=700,
    node_color='white',
    edgecolors='gray'
)

bar_lt = (
    res.drop_duplicates(['names', 'lesion_type'])[['contrast', 'names', 'ligand', 'receptor', 'lesion_type']]
    .assign(inter=lambda x: x['ligand'] + '^' + x['receptor'])
    .groupby(['inter', 'lesion_type'])
    .count()
    .reset_index()
    [['inter', 'lesion_type', 'names']]
)

bar_lt = bar_lt.pivot(index='inter', columns='lesion_type', values='names').loc[:, ['Ctrl', 'CA', 'CI']] # 'CI' np.flip(order)

source_tab = (
    res.drop_duplicates(['names', 'lesion_type'])[['contrast', 'names', 'source', 'target', 'ligand', 'receptor']]
    .assign(inter=lambda x: x['ligand'] + '^' + x['receptor'])
    .groupby(['inter', 'source'])
    .count()
    .reset_index()
    [['inter', 'source', 'target']]
    .set_index('inter')
)
source_tab = source_tab.reset_index().pivot(index='inter', columns='source', values='target')#.loc[order]

target_tab = (
    res.drop_duplicates(['names', 'lesion_type'])[['contrast', 'names', 'source', 'target', 'ligand', 'receptor']]
    .assign(inter=lambda x: x['ligand'] + '^' + x['receptor'])
    .groupby(['inter', 'target'])
    .count()
    .reset_index()
    [['inter', 'target', 'source']]
    .set_index('inter')
)
target_tab = target_tab.reset_index().pivot(index='inter', columns='target', values='source')#.loc[order]

msk = target_tab.index[target_tab.sum(1) > 1]
bar_lt = bar_lt.loc[np.flip(msk)]
source_tab = source_tab.loc[msk]
target_tab = target_tab.loc[msk]

# Sort
msk = target_tab.sum(1).sort_values(ascending=False).index
bar_lt = bar_lt.loc[np.flip(msk)]
source_tab = source_tab.loc[msk]
target_tab = target_tab.loc[msk]

fig5, axes = plt.subplots(1, 3, figsize=(6.5, 6.5), tight_layout=True, dpi=150) # gridspec_kw={'width_ratios': [1.5, 1.5, 2, 1]}
ax = axes[0]
sns.heatmap(
    source_tab,
    cmap='viridis',
    annot=True,
    ax=ax,
    cbar=False,
    yticklabels=True,
    xticklabels=True,
    vmax=4.0
)
ax.set_facecolor('xkcd:gray')
ax.xaxis.set_tick_params(rotation=90)
ax.set_ylabel('')
ax = axes[1]
sns.heatmap(
    target_tab,
    cmap='viridis',
    annot=True,
    ax=ax,
    cbar=False,
    xticklabels=True,
    vmax=4.0,
)
ax.set_ylabel('')
ax.set_yticks([])
ax.xaxis.set_tick_params(rotation=90)
ax.set_facecolor('xkcd:gray')

ax = axes[2]
(
    bar_lt
    .plot(kind='barh', ax=ax, width=0.9, stacked=True, align='center', color=palette)
)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
ax.set_xlabel('# cell type pairs')
ax.set_ylabel('')
ax.set_yticks([])
ax.margins(y=0)
ax.grid(axis='x')
ax.set_axisbelow(True)
ax.set_xticks([0.0, 2.0, 4.0, 6.0])

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

# Compute mean scores
mean_scores = dc.get_pseudobulk(
    adata=scores,
    sample_col='niches',
    groups_col=None,
    mode='mean',
    min_cells=0,
    min_counts=0,
)
mean_scores = mean_scores[(mean_scores.obs['niches'] != 'Ependym') & (mean_scores.obs['niches'] != 'GM')].copy()
mean_scores.var['inter'] = ['^'.join(x.split('^')[2:]) for x in mean_scores.var.index]
mean_scores = mean_scores.T.copy()
mean_scores.obs = mean_scores.obs.reset_index(names='names').set_index('inter')
mean_scores = mean_scores.to_df().reset_index().groupby('inter').mean(numeric_only=True)
mean_scores.loc[:, :] = sc.pp.scale(mean_scores.values.T).T

cg = sns.clustermap(
    mean_scores.loc[source_tab.index, :], #res['names'].unique()
    xticklabels=True,
    yticklabels=True,
    cmap='coolwarm',
    vmin=-1,
    vmax=+1,
    figsize=(3, 6),
    cbar_pos=(1, .4, .03, .2),
    dendrogram_ratio=(0.25, 0.05)
)
cg.ax_heatmap.set_xticklabels(cg.ax_heatmap.get_xticklabels(), rotation=90)
cg.fig.set_dpi(150)
fig6 = cg.fig

# Save to pdf
pdf = matplotlib.backends.backend_pdf.PdfPages(plot_path)
for fig in [fig1, fig2, fig3, fig4, fig5, fig6]:
    pdf.savefig(fig, bbox_inches='tight')
pdf.close()
