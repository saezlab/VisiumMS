import pandas as pd
import numpy as np
import decoupler as dc
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--meta_path', required=True)
parser.add_argument('-b','--colors_dict', required=True)
parser.add_argument('-c','--plot_path', required=True)
args = vars(parser.parse_args())

meta_path = args['meta_path']
colors_dict = args['colors_dict']
plot_path = args['plot_path']

# Get palette
palette = dict(item.split('_') for item in colors_dict.strip("'").split(';'))

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

res = []
for sample_id in vs_samples:
    print(sample_id)
    df = read_annots(sample_id)
    cats = np.intersect1d(df['areas'].values.astype('U'), df['niches'].values.astype('U'))
    for cat in cats:
        a = (df['niches'] == cat).values
        b = (df['areas'] == cat).values
        j = jaccard(a, b)
        res.append([sample_id, cat, j])
res = pd.DataFrame(res, columns=['sample_id', 'cat', 'jacc'])

# Plot
fig, ax = plt.subplots(1, 1, figsize=(2, 2), dpi=150)
order = res.groupby('cat')['jacc'].median().sort_values(ascending=False).index
sns.boxplot(data=res, x='cat', y='jacc', ax=ax, fliersize=5, order=order, hue='cat', legend=False, palette=palette)
ax.tick_params(axis='x', labelrotation=90)
ax.set_ylabel('Jaccard index')
ax.set_ylim(0, 1)
ax.set_xlabel('')

# Write
fig.savefig(plot_path, bbox_inches='tight')
