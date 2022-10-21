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
fig_path = 'figures/manuscript/fig2/'
sfig_path = 'figures/manuscript/supp_fig2/'
fig_name = 'pca_vm.pdf'
plt.rcParams['font.sans-serif'] = 'Arial'


def read_slide(sample_id):

    # Read rna-seq
    slide = sc.read_visium('data/raw/visium/{0}/outs/'.format(sample_id))
    slide.var_names_make_unique()
    sc.pp.filter_genes(slide, min_cells=3)
    sc.pp.filter_cells(slide, min_genes=200)

    # Store raw counts
    slide.layers['counts'] = slide.X.copy()

    # Normalize
    sc.pp.normalize_total(slide, target_sum=1e4)
    sc.pp.log1p(slide)

    # Read props
    props = pd.read_csv('data/prc/visium/{0}/cell_props.csv'.format(sample_id), index_col=0)
    inter = slide.obs.index.intersection(props.index)
    slide.obsm['props'] = props.loc[inter]

    return slide


# Read slides
meta = pd.read_csv(os.path.join('data','metadata.csv')).set_index('sample_id')
slides = []
for sample_id in meta.index.values:
    slide = read_slide(sample_id)
    slide.obs['sample_id'] = sample_id
    slide.obs['lesion_type'] = meta.loc[sample_id, 'lesion_type']
    slides.append(slide)
    
slides = slides[0].concatenate(slides[1:], join='outer')

# Pseudobulk
padata = dc.get_pseudobulk(slides, sample_col='sample_id', layer='counts', groups_col=None, min_prop=0.2, min_smpls=3, use_raw=False)

# Normalize
sc.pp.normalize_total(padata, target_sum=1e4)
sc.pp.log1p(padata)

# PCA
sc.pp.scale(padata, max_value=10)
pca = PCA(n_components=padata.X.shape[0] - 1)
pca_x = pca.fit_transform(padata.X)
explained_variance = pca.explained_variance_ratio_

# Read palette
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

sns.scatterplot(data=df, x="PC 2", y="PC 3", hue="lesion_type", style="sex", palette=palette, ax=ax)
ax.set_xlabel('PC 2 | Exp. var. {0:.2f}'.format(explained_variance[1]))
ax.set_ylabel('PC 3 | Exp. var. {0:.2f}'.format(explained_variance[2]))

# Reorder legens
handles, labels = ax.get_legend_handles_labels()
#labels = ['lesion_type', 'Control', 'Acute', 'Chronic Active', 'Chronic Inactive', 'sex', 'F', 'M']
#idx = np.array([labels.index(i) for i in labels])
#handles, labels = np.array(handles)[idx], np.array(labels)[idx]
ax.legend(handles, labels, bbox_to_anchor=(1,0.5), loc='center left', frameon=False)
os.makedirs(fig_path, exist_ok=True)
fig.savefig(os.path.join(fig_path, fig_name), bbox_inches='tight')

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
os.makedirs(sfig_path, exist_ok=True)
fig.savefig(os.path.join(sfig_path, fig_name), bbox_inches='tight')
