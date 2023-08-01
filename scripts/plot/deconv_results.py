import pandas as pd
import numpy as np
import scanpy as sc
import os
from anndata import AnnData
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp


path_slides = 'data/raw/visium/'
path_deconv = 'data/prc/visium/'
output_path = 'figures/deconv/'

# Read annotated sc atlas
meta = pd.read_csv('data/metadata.csv').set_index('sample_id')

# Read deconv results
prop = []
sample_ids = []
for sample_id in meta.index:
    sample_ids.append(sample_id)
    tmp = AnnData(pd.read_csv(os.path.join(path_deconv, '{0}/cell_props.csv'.format(sample_id)), index_col=0))
    tmp.obs['sample_id'] = sample_id
    prop.append(tmp)
prop = prop[0].concatenate(prop[1:])

# Define params
n_celltypes = prop.var_names.shape[0]
n_rows = int(np.ceil(n_celltypes / 4))

# Plot proportions
for sample_id in sample_ids:
    # Read original slide
    slide = sc.read_visium(os.path.join(path_slides, '{0}/outs/'.format(sample_id)))
    tmp = prop[prop.obs['sample_id'] == sample_id]
    tmp.uns['spatial'] = slide.uns['spatial'].copy()
    tmp.obsm['spatial'] = slide.obsm['spatial'].copy()
    del slide
    
    # Plot props in slide
    fig, axes = plt.subplots(n_rows, 4, tight_layout=False, facecolor='white', figsize=(4*4, 3*n_rows))
    axes = axes.flatten()
    for celltype, ax in zip(prop.var_names, axes):
        sc.pl.spatial(tmp, color=celltype, size=2, frameon=False, cmap='Reds', return_fig=False, show=False, ax=ax)
    fig.suptitle('{0} proportions'.format(sample_id))
    os.makedirs(output_path, exist_ok=True)
    fig.savefig(os.path.join(output_path, '{0}.pdf'.format(sample_id)))
    del tmp
