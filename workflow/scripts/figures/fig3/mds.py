import pandas as pd
import numpy as np
import scanpy as sc
import decoupler as dc
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import seaborn as sns
from composition_stats import closure
from sklearn.manifold import MDS
import adjustText as at
import anndata as ad
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-b','--data_path', required=True)
parser.add_argument('-c','--meta_path', required=True)
parser.add_argument('-d','--sn_mofa_path', required=True)
parser.add_argument('-e','--vs_mofa_path', required=True)
parser.add_argument('-f','--plot_path', required=True)
args = vars(parser.parse_args())

data_path = args['data_path']
meta_path = args['meta_path']
sn_mofa_path = args['sn_mofa_path']
vs_mofa_path = args['vs_mofa_path']
plot_path = args['plot_path']

# Read
adata = sc.read_h5ad(data_path)
meta = pd.read_csv(meta_path)
vs_samples = meta[~meta['Batch vs'].isnull()]['Sample id'].values.astype('U')

def clr(mat):
    msk = mat.astype(bool)
    lmat = np.log(mat, where=msk)
    gm = lmat.mean(axis=-1, keepdims=True, where=msk)
    return np.where(msk, (lmat - gm).squeeze(), 0.)

# Compute clr props sn
sn_prop = (
    adata
    .obs
    .groupby(['Sample id', 'leiden'])
    .count()
    [['total_counts']]
    .reset_index()
    .pivot(index='Sample id', columns='leiden', values='total_counts')
)
sn_prop.loc[:, :] = clr(closure(sn_prop.values))

# Compute clr props cell states
ctypes = adata.obs['leiden'].unique().astype(str)
cs_prop = []
for ctype in ctypes:
    tmp = pd.read_csv('data/prc/ctypes/{ctype}_ann.csv'.format(ctype=ctype), index_col=0)
    tmp['Sample id'] = [i.split('-')[0] for i in tmp.index]
    tmp['total_counts'] = 1
    cs_prop.append(tmp)

cs_prop = (
    pd.concat(cs_prop)
    .groupby(['Sample id', 'cell_states'])
    .count()
    [['total_counts']]
    .reset_index()
    .pivot(index='Sample id', columns='cell_states', values='total_counts')
    .fillna(0)
)
cs_prop.loc[:, :] = clr(closure(cs_prop.values))

# Compute clr props visium
vs_prop = []
for sample_id in vs_samples:
    tmp = (
        pd.read_csv('data/prc/vs/{0}/abunds.csv'.format(sample_id), index_col=0)
        .sum(0)
        .reset_index(name=sample_id)
        .set_index('index')
        .T
        .fillna(0)
    )
    vs_prop.append(tmp)
vs_prop = pd.concat(vs_prop)
vs_prop.loc[:, :] = clr(closure(vs_prop.values))

# Read mofa results
sn_mofa = pd.read_csv(sn_mofa_path, index_col=0).pivot(
    index='Sample id', columns='Factor', values='value')
vs_mofa = pd.read_csv(vs_mofa_path, index_col=0).pivot(
    index='Sample id', columns='Factor', values='value')

def get_distances(df):
    cors = 1 - df.T.corr()
    return cors

def get_coords(cors, meta):
    mds = MDS(n_components=2, normalized_stress='auto', n_jobs=-1, dissimilarity='precomputed', random_state=0)
    coords = mds.fit_transform(cors.values)
    coords = pd.DataFrame(coords, columns=['MDS1', 'MDS2'])
    coords['Sample id'] = cors.index
    coords = pd.merge(coords, meta[['Sample id', 'Lesion type', 'Sex']])
    return coords


def plot_mds(coords, title, figsize=(4, 4)):
    fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=150)
    sns.scatterplot(data=coords, x="MDS1", y="MDS2", hue="Lesion type", ax=ax)
    texts = []
    for x, y, t in zip(coords['MDS1'], coords['MDS2'], coords['Sample id']):
            texts.append(ax.text(x, y, t))
    at.adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black'), ax=ax)
    # Reorder legens
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, bbox_to_anchor=(1,0.5), loc='center left', frameon=False)
    ax.set_title(title)
    return fig

# Get distances
sn_prop_ds = get_distances(sn_prop)
cs_prop_ds = get_distances(cs_prop)
vs_prop_ds = get_distances(vs_prop)
sn_mofa_ds = get_distances(sn_mofa)
vs_mofa_ds = get_distances(vs_mofa)

# Get final distances
total_ds = (
    sn_prop_ds
    .add(cs_prop_ds, fill_value=0)
    .add(vs_prop_ds, fill_value=0)
    .add(sn_mofa_ds, fill_value=0)
    .add(vs_mofa_ds, fill_value=0)
)
msk = np.isnan(total_ds).sum() == 0
total_ds = total_ds.loc[msk, msk]

# Get MDS coords
sn_prop_coo = get_coords(sn_prop_ds, meta)
cs_prop_coo = get_coords(cs_prop_ds, meta)
vs_prop_coo = get_coords(vs_prop_ds, meta)
sn_mofa_coo = get_coords(sn_mofa_ds, meta)
vs_mofa_coo = get_coords(vs_mofa_ds, meta)
total_coo = get_coords(total_ds, meta)

fg1 = plot_mds(sn_prop_coo, title='sn_prop')
fg2 = plot_mds(cs_prop_coo, title='cs_prop')
fg3 = plot_mds(vs_prop_coo, title='vs_prop')
fg4 = plot_mds(sn_mofa_coo, title='sn_mofa')
fg5 = plot_mds(vs_mofa_coo, title='vs_mofa')
fg6 = plot_mds(total_coo, title='total')

# Save to pdf
pdf = matplotlib.backends.backend_pdf.PdfPages(plot_path)
for fig in [fg1, fg2, fg3, fg4, fg5, fg6]:
    pdf.savefig(fig, bbox_inches='tight')
pdf.close()