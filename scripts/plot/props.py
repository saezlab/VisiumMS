import pandas as pd
import numpy as np
import scanpy as sc
import os
import seaborn as sns
import matplotlib.pyplot as plt


# Read raw sc object
adata = sc.read('data/prc/sc/raw.h5ad')

# Compute props for sc
sc = (adata.obs
 .groupby(['sample_id', 'leiden'])
 .count()
 .reset_index()
 [['sample_id', 'leiden', 'total_counts']]
)
sc['props'] = (sc['total_counts'] /
               sc.groupby(['sample_id'])['total_counts']
               .transform('sum'))

# Add lesion type labels 
meta = pd.read_csv('data/metadata.csv').set_index('sample_id')
sc['lesion_type'] = [meta.loc[s, 'lesion_type'] for s in sc['sample_id']]

# Filter by lesion type
meta = meta[np.isin(meta['lesion_type'], ['Chronic Active', 'Control'])]
sc = sc[np.isin(sc['lesion_type'], ['Chronic Active', 'Control'])]

# Plot
g = (sns.catplot(col='leiden', y='props', x='lesion_type', data=sc, col_wrap=5, kind='box', sharey=False, height=2)
 .set_titles(template='{col_name}').set_xlabels('').set_ylabels('Proportions').set_xticklabels(rotation=45)
)
g.fig.subplots_adjust(top=0.8)
g.fig.suptitle('Cell type proportions in single-cell')
res_path = 'figures/props'
os.makedirs(res_path, exist_ok=True)
plt.savefig('{0}/sc.pdf'.format(res_path), bbox_inches='tight')

# Extract proportions from visium slides
vs = []
for sample_id in meta.index:
    prop = pd.read_csv('data/prc/visium/{0}/cell_props.csv'.format(sample_id), index_col=0)
    prop = pd.DataFrame(prop.mean(axis=0), columns=['props']).reset_index().rename({'index': 'leiden'}, axis=1)
    prop['sample_id'] = sample_id
    prop['lesion_type'] = meta.loc[sample_id, 'lesion_type']
    vs.append(prop)
vs = pd.concat(vs)

# Plot
g = (sns.catplot(col='leiden', y='props', x='lesion_type', data=vs, col_wrap=5, kind='box', sharey=False, height=2)
 .set_titles(template='{col_name}').set_xlabels('').set_ylabels('Proportions').set_xticklabels(rotation=45)
)
g.fig.subplots_adjust(top=0.8)
g.fig.suptitle('Cell type proportions in visium')
plt.savefig('{0}/vm.pdf'.format(res_path), bbox_inches='tight')
