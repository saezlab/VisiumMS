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
parser.add_argument('-b','--colors_dict', required=True)
parser.add_argument('-c','--plot_path', required=True)
args = vars(parser.parse_args())

ctlr_path = args['ctlr_path']
colors_dict = args['colors_dict']
plot_path = args['plot_path']

# Get palette
palette = dict(item.split(':') for item in colors_dict.strip("'").split(';'))

# Read
res = pd.read_csv(ctlr_path)

# Number of diff interactions per cell type
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

fig1, ax = plt.subplots(1, 1, figsize=(4, 2), tight_layout=True, dpi=150)
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

# Venn diagram
tab = res.groupby(['lesion_type'])['names'].apply(lambda x: set(x))
labels = tab.index
sets = tab.values

fig2, ax = plt.subplots(1, 1, figsize=(4, 4), dpi=150)
v = venn3(sets, labels, set_colors=[palette[l] for l in labels], alpha=0.75, ax=ax)
c = venn3_circles(sets, linewidth=1, ax=ax)

# Networks
g = res[['source', 'target', 'names', 'lesion_type']].groupby(['source', 'target', 'lesion_type']).count().reset_index().rename(columns={'names': 'n'})
g['edge'] = ['^'.join(sorted([x, y])) for x, y in zip(g['source'], g['target'])]
g = g.groupby(['edge', 'lesion_type']).sum(numeric_only=True).reset_index()
g[['source', 'target']] = g['edge'].str.split('^', expand=True, regex=False)[[0, 1]]
g = g[['source', 'target', 'lesion_type', 'n']]

def plot_graph(df, lesion_type, layout, ax, palette):
    if lesion_type is not None:
        df = df[df['lesion_type'] == lesion_type]
    
    G = nx.from_pandas_edgelist(
        df,
        source='source',
        target='target',
        edge_attr=['lesion_type', 'n']
    )
    G.add_nodes_from([k for k in layout if k not in G.nodes()])
    edges = G.edges()
    nodes = G.nodes()
    weights = [G[u][v]['n'] for u,v in edges]
    weights = np.log2(np.array(weights) + 1) + 1
    if lesion_type is not None:
        edge_color = [palette[lesion_type] for i in range(len(weights))]
    else:
        edge_color = [palette[G[u][v]['lesion_type']] for u, v in edges]
    nx.draw(
        G,
        pos=layout,
        with_labels=True,
        edge_color=edge_color,
        width=weights,
        ax=ax,
        alpha=1,
        node_size=700,
        node_color='white',
        edgecolors='gray',
    )

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

fig3, ax = plt.subplots(1, 1, figsize=(4, 4), dpi=150)
plot_graph(g, 'Ctrl', layout, ax, palette)

fig4, ax = plt.subplots(1, 1, figsize=(4, 4), dpi=150)
plot_graph(g, 'CA', layout, ax, palette)

fig5, ax = plt.subplots(1, 1, figsize=(4, 4), dpi=150)
plot_graph(g, 'CI', layout, ax, palette)

g.loc[g['lesion_type'] == 'CA', 'n'] = -(g.loc[g['lesion_type'] == 'CA', 'n'])
g.loc[g['lesion_type'] == 'CI', 'n'] = -(g.loc[g['lesion_type'] == 'CI', 'n'])
g_df = g.groupby(['source', 'target']).sum(numeric_only=True).reset_index()
g_df = g_df[abs(g_df['n']) > 1]

g = nx.from_pandas_edgelist(g_df, source='source', target='target', edge_attr=['n']) # create_using=nx.DiGraph

c_dict = {'CA': 'red', 'Ctrl': 'green'}
edges = g.edges()
colors = ['red' if g[u][v]['n'] < 0 else 'green' for u,v in edges]
weights = [abs(g[u][v]['n']) for u,v in edges]
weights = np.log2(np.array(weights) + 1) + 1

fig6, ax = plt.subplots(1, 1, figsize=(3, 3), tight_layout=True, dpi=150)
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

# Biomarkers
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

fig7, ax = plt.subplots(1, 1, figsize=(3, 5), tight_layout=True, dpi=150)
(
    bar
    .pivot(index='value', columns='lesion_type', values='variable')
    .loc[order[-30:]]
    .plot(kind='barh', ax=ax, width=0.9, stacked=True, color=palette)
)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
ax.set_ylabel('')
ax.set_xlabel('# interactions')
ax.grid(axis='x')
ax.set_axisbelow(True)

# Save to pdf
pdf = matplotlib.backends.backend_pdf.PdfPages(plot_path)
for fig in [fig1, fig2, fig3, fig4, fig5, fig6, fig7]:
    pdf.savefig(fig, bbox_inches='tight')
pdf.close()
