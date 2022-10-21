import scanpy as sc
import pandas as pd
import numpy as np
import decoupler as dc
import matplotlib.pyplot as plt
import seaborn as sns
import os
from sklearn.decomposition import PCA
import adjustText as at


# Defina path
fig_path = 'figures/manuscript/supp_fig2/'
fig_name = 'pca_sc.pdf'
plt.rcParams['font.sans-serif'] = 'Arial'


# Read adata
meta = pd.read_csv(os.path.join('data','metadata.csv')).set_index('sample_id')
adata = sc.read_h5ad('data/prc/sc/raw.h5ad')
meta = meta[np.isin(meta.index.values, np.unique(adata.obs['sample_id'].values))]

# Pseudobulk
padata = dc.get_pseudobulk(adata, sample_col='sample_id', groups_col=None, min_prop=0.2, min_smpls=3, use_raw=False)

# Normalize
sc.pp.normalize_total(padata, target_sum=1e4)
sc.pp.log1p(padata)

# PCA
sc.pp.scale(padata, max_value=10)
pca = PCA(n_components=2)
pca_x = pca.fit_transform(padata.X)
explained_variance = pca.explained_variance_ratio_

# Plot
palette = pd.read_csv('data/cond_palette.csv')
palette['rgb'] = [(r, g, b) for r, g, b in zip(palette['r'], palette['g'], palette['b'])]
palette = {k: v for k,v in zip(palette['cond'], palette['rgb'])}

# Extract PCs
df = pd.DataFrame(pca_x, columns=['PC {0}'.format(i + 1) for i in range(pca_x.shape[1])])
df['sample_id'] = padata.obs['sample_id'].values.copy()
df['lesion_type'] = padata.obs['lesion_type'].values.copy()
df['sex'] = [meta.loc[i, 'sex'] for i in meta.index]

# Plot
fig, ax = plt.subplots(1, 1, figsize=(4.5, 2.5), facecolor='white', dpi=150, tight_layout=True)

sns.scatterplot(data=df, x="PC 1", y="PC 2", hue="lesion_type", style="sex", palette=palette, ax=ax)
ax.set_xlabel('PC 1 | Exp. var. {0:.2f}'.format(explained_variance[0]))
ax.set_ylabel('PC 2 | Exp. var. {0:.2f}'.format(explained_variance[1]))

# Reorder legens
handles, labels = ax.get_legend_handles_labels()
#labels = ['lesion_type', 'Control', 'Acute', 'Chronic Active', 'Chronic Inactive', 'sex', 'F', 'M']
#idx = np.array([labels.index(i) for i in labels])
#handles, labels = np.array(handles)[idx], np.array(labels)[idx]
ax.legend(handles, labels, bbox_to_anchor=(1,0.5), loc='center left', frameon=False)
os.makedirs(fig_path, exist_ok=True)
fig.savefig(os.path.join(fig_path, fig_name), bbox_inches='tight')
