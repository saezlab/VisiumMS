import pandas as pd
import numpy as np
import decoupler as dc
import scanpy as sc
from sklearn.metrics import adjusted_rand_score
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import seaborn as sns
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--meta_path', required=True)
parser.add_argument('-b','--n_colors_dict', required=True)
parser.add_argument('-c','--l_colors_dict', required=True)
parser.add_argument('-d','--plot_path', required=True)
args = vars(parser.parse_args())

meta_path = args['meta_path']
n_colors_dict = args['n_colors_dict']
l_colors_dict = args['l_colors_dict']
plot_path = args['plot_path']

# Get palette
n_palette = dict(item.split('_') for item in n_colors_dict.strip("'").split(';'))
l_palette = dict(item.split(':') for item in l_colors_dict.strip("'").split(';'))

# Read meta
meta = pd.read_csv(meta_path)
vs_samples = meta[~meta['Batch vs'].isnull()]['Sample id'].values.astype('U')

def read_annots(sample_id):
    slide = sc.read_h5ad('data/prc/vs/{0}/adata.h5ad'.format(sample_id))
    slide.obs['niches'] = pd.read_csv('data/prc/vs/{0}/niches.csv'.format(sample_id), index_col=0)
    df = slide.obs[['areas', 'niches']].copy()
    return df
    
def jaccard(a, b):
    inter = np.sum(a * b)
    union = np.sum(a | b)
    return inter / union

ari_df = []
res = []
for sample_id in vs_samples:
    print(sample_id)
    df = read_annots(sample_id)
    cats = np.intersect1d(df['areas'].values.astype('U'), df['niches'].values.astype('U'))
    df = df[df['niches'].isin(cats) & df['areas'].isin(cats)].copy()  # Filter by shared cats
    for cat in cats:
        a = (df['niches'] == cat).values
        b = (df['areas'] == cat).values
        j = jaccard(a, b)
        res.append([sample_id, cat, j])
    ari = adjusted_rand_score(df['areas'].values.astype('U'), df['niches'].values.astype('U'))
    ari_df.append([sample_id, ari])
res = pd.DataFrame(res, columns=['sample_id', 'cat', 'jacc'])
ari_df = pd.DataFrame(ari_df, columns=['sample_id', 'ari'])

s, p = scipy.stats.ttest_1samp(
    a=ari_df['ari'].values,
    popmean=0,
    alternative='greater'
)
print('ARI:', p)

# Plot
fig1, ax = plt.subplots(1, 1, figsize=(2, 2), dpi=150)
order = res.groupby('cat')['jacc'].median().sort_values(ascending=False).index
sns.boxplot(data=res, x='cat', y='jacc', ax=ax, fliersize=5, order=order, hue='cat', legend=False, palette=n_palette)
ax.tick_params(axis='x', labelrotation=90)
ax.set_ylabel('Jaccard index')
ax.set_ylim(-0.05, 1.05)
ax.set_xlabel('')

fig2, ax = plt.subplots(1, 1, figsize=(2, 2), dpi=150)
data = pd.merge(ari_df, meta[['Sample id', 'Lesion type']].rename(columns={'Sample id': 'sample_id'}))
sns.boxplot(data=data, x='Lesion type', hue='Lesion type', legend=False, y='ari', ax=ax, fliersize=5, palette=l_palette)
ax.tick_params(axis='x', labelrotation=90)
ax.set_ylim(-0.05, 1.05)
ax.set_xlabel('')
ax.set_ylabel('Adjusted Rand Index')

# Save to pdf
pdf = matplotlib.backends.backend_pdf.PdfPages(plot_path)
for fig in [fig1, fig2]:
    pdf.savefig(fig, bbox_inches='tight')
pdf.close()
