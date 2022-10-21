import pandas as pd
import numpy as np
import scanpy as sc
import os
from anndata import AnnData
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp


path_slides = 'data/raw/visium/'
path_ann = 'data/prc/sc/annotated.h5ad'
path_deconv = 'data/prc/visium/'
fig_path = 'figures/manuscript/supp_fig1/'
fig_name_samples = 'corr_deconvsamples.pdf'
fig_name_celltypes = 'corr_deconvcelltypes.pdf'
plt.rcParams['font.sans-serif'] = 'Arial'

# Read annotated sc atlas
meta = pd.read_csv('data/metadata.csv').set_index('sample_id')

# Compute proportions for sc
sc_df = (sc.read(path_ann).obs
 .groupby(['sample_id', 'cell_type'])[['sample_id']]
 .count()
 .rename({'sample_id': 'n'}, axis=1)
 .reset_index()
 .assign(sc=lambda x: x['n'] / x.groupby('sample_id')['n'].transform('sum'))
 .drop('n', axis=1)
)

# Read deconv results
prop = []
sample_ids = []
for sample_id in meta.index:
    sample_ids.append(sample_id)
    tmp = AnnData(pd.read_csv(os.path.join(path_deconv, '{0}/cell_props.csv'.format(sample_id)), index_col=0))
    tmp.obs['sample_id'] = sample_id
    prop.append(tmp)
prop = prop[0].concatenate(prop[1:])

# Compute average props for visium
vm_df = pd.DataFrame(prop.X, index=prop.obs.index, columns=prop.var.index)
vm_df['sample_id'] = prop.obs['sample_id']
vm_df = vm_df.groupby('sample_id').mean(1)
vm_df = vm_df.melt(value_vars=vm_df.columns, ignore_index=False).reset_index()
vm_df = vm_df.rename({'variable': 'cell_type', 'value': 'vm'}, axis=1)

# Merge dataframes
df = pd.merge(sc_df, vm_df, how='inner')
df['sc'] = np.log10(df['sc'].replace({0: np.nan}))
df['vm'] = np.log10(df['vm'].replace({0: np.nan}))
df = df.fillna(-4)

# Plot correlation samples
def annotate(data, **kws):
    x = data['sc'].values
    y = data['vm'].values
    msk = np.isfinite(x) & np.isfinite(y)
    r, p = sp.stats.pearsonr(x[msk], y[msk])
    ax = plt.gca()
    ax.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(r, p), transform=ax.transAxes)
    
# Plot corrs per sample
g = sns.lmplot(x='sc', y='vm', data=df, col='sample_id', col_wrap=4, height=3)
g.map_dataframe(annotate)
g.fig.set_facecolor('white')
g.fig.suptitle('Correlation per sample')
g.fig.set_tight_layout(tight='pad')
os.makedirs(fig_path, exist_ok=True)
g.fig.savefig(os.path.join(fig_path, fig_name_samples), bbox_inches='tight')

# Plot corrs per cell type
g = sns.lmplot(x='sc', y='vm', data=df, col='cell_type', col_wrap=4, height=3)
g.map_dataframe(annotate)
g.fig.set_facecolor('white')
g.fig.suptitle('Correlation per sample')
g.fig.set_tight_layout(tight='pad')
os.makedirs(fig_path, exist_ok=True)
g.fig.savefig(os.path.join(fig_path, fig_name_celltypes), bbox_inches='tight')
