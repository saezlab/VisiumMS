import pandas as pd
import numpy as np
import scanpy as sc
import os
from anndata import AnnData
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp
import argparse


# Read command line and set args
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--path_slides', help='Input path to raw slides', required=True)
parser.add_argument('-a', '--path_ann', help='Annotated AnnData object', required=True)
parser.add_argument('-d', '--path_deconv', help='Path to deconv results', required=True)
parser.add_argument('-o', '--output_path', help='Path were to save plots', required=True)
args = vars(parser.parse_args())

path_slides = args['path_slides']
path_ann = args['path_ann']
path_deconv = args['path_deconv']
output_path = args['output_path']
###############################

# Read annotated sc atlas
meta = sc.read(path_ann).obs

# Compute proportions for sc
sc_df = (meta
 .groupby(['sample_id', 'leiden'])[['sample_id']]
 .count()
 .rename({'sample_id': 'n'}, axis=1)
 .reset_index()
 .assign(sc=lambda x: x['n'] / x.groupby('sample_id')['n'].transform('sum'))
 .drop('n', axis=1)
)

# Read deconv results
prop = []
sample_ids = []
for sample_id in os.listdir(path_deconv):
    if os.path.isdir(os.path.join(path_deconv, sample_id)) and not sample_id.startswith('.') and sample_id != 'reg_model':
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
vm_df = vm_df.rename({'variable': 'leiden', 'value': 'vm'}, axis=1)

# Merge dataframes
df = pd.merge(sc_df, vm_df, how='inner')
df['sc'] = np.log10(df['sc'].replace({0: np.nan}))
df['vm'] = np.log10(df['vm'].replace({0: np.nan}))
df = df.fillna(-4)

# Plot correlations:
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
g.fig.savefig(os.path.join(output_path, 'deconv_corr_sample.pdf'), dpi=300)

# Plot corrs per cell type
g = sns.lmplot(x='sc', y='vm', data=df, col='leiden', col_wrap=4, height=3)
g.map_dataframe(annotate)
g.fig.set_facecolor('white')
g.fig.suptitle('Correlation per sample')
g.fig.set_tight_layout(tight='pad')
g.fig.savefig(os.path.join(output_path, 'deconv_corr_celltype.pdf'), dpi=300)

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
    fig.savefig(os.path.join(output_path, 'deconv_props_{0}.pdf'.format(sample_id)))
    del tmp
