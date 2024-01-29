import pandas as pd
import numpy as np
import scanpy as sc
import decoupler as dc
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import seaborn as sns
from composition_stats import closure
from sklearn.manifold import MDS
from sklearn.metrics import silhouette_score, silhouette_samples
import adjustText as at
import anndata as ad
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--colors_dict', required=True)
parser.add_argument('-b','--data_path', required=True)
parser.add_argument('-c','--meta_path', required=True)
parser.add_argument('-d','--sn_mofa_path', required=True)
parser.add_argument('-e','--vs_mofa_path', required=True)
parser.add_argument('-f','--plot_path', required=True)
args = vars(parser.parse_args())

colors_dict = args['colors_dict']
data_path = args['data_path']
meta_path = args['meta_path']
sn_mofa_path = args['sn_mofa_path']
vs_mofa_path = args['vs_mofa_path']
plot_path = args['plot_path']

# Get palette
palette = dict(item.split(':') for item in colors_dict.strip("'").split(';'))

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


def plot_mds(coords, title, figsize=(4, 4), palette=None):
    fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=150)
    sns.scatterplot(data=coords, x="MDS1", y="MDS2", hue="Lesion type", ax=ax, palette=palette)
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

# Compute Sil. score
def get_sil_scores(df, name, meta):
    labels = meta.set_index('Sample id').loc[df.index, 'Lesion type'].values.astype('U')
    scores = silhouette_samples(df, labels, metric='precomputed')
    scores = pd.DataFrame(scores.reshape(-1, 1), columns=['score'])
    scores['name'] = name
    return scores

sn_prop_sil = get_sil_scores(sn_prop_ds, 'sn_prop', meta)
cs_prop_sil = get_sil_scores(cs_prop_ds, 'cs_prop', meta)
vs_prop_sil = get_sil_scores(vs_prop_ds, 'vs_prop', meta)
sn_mofa_sil = get_sil_scores(sn_mofa_ds, 'sn_mofa', meta)
vs_mofa_sil = get_sil_scores(vs_mofa_ds, 'vs_mofa', meta)
total_sil = get_sil_scores(total_ds, 'total', meta)
sil_df = pd.concat([sn_prop_sil, cs_prop_sil, vs_prop_sil, sn_mofa_sil, vs_mofa_sil, total_sil])

# Get MDS coords
sn_prop_coo = get_coords(sn_prop_ds, meta)
cs_prop_coo = get_coords(cs_prop_ds, meta)
vs_prop_coo = get_coords(vs_prop_ds, meta)
sn_mofa_coo = get_coords(sn_mofa_ds, meta)
vs_mofa_coo = get_coords(vs_mofa_ds, meta)
total_coo = get_coords(total_ds, meta)

# Supps
fg1 = plot_mds(sn_prop_coo, title='sn_prop', palette=palette)
fg2 = plot_mds(cs_prop_coo, title='cs_prop', palette=palette)
fg3 = plot_mds(vs_prop_coo, title='vs_prop', palette=palette)
fg4 = plot_mds(sn_mofa_coo, title='sn_mofa', palette=palette)
fg5 = plot_mds(vs_mofa_coo, title='vs_mofa', palette=palette)
fg6 = plot_mds(total_coo, title='total', palette=palette)

# Main figure
fg0, ax = plt.subplots(1, 1, figsize=(4., 2.5), facecolor='white', dpi=150, tight_layout=True)
sns.scatterplot(data=total_coo, x="MDS1", y="MDS2", hue="Lesion type", style="Sex", ax=ax, palette=palette)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, bbox_to_anchor=(1,0.5), loc='center left', frameon=False)

# Distirbution of sil scores
fg7, ax = plt.subplots(1, 1, figsize=(3, 3), facecolor='white', dpi=150, tight_layout=True)
sns.boxplot(data=sil_df, x='name', y='score', fliersize=5, ax=ax)
ax.set_xlabel('')
ax.set_ylabel('Silhouette Coeff')
ax.tick_params(axis='x', labelrotation=90)

# Save to pdf
pdf = matplotlib.backends.backend_pdf.PdfPages(plot_path)
for fig in [fg0, fg1, fg2, fg3, fg4, fg5, fg6, fg7]:
    pdf.savefig(fig, bbox_inches='tight')
pdf.close()
