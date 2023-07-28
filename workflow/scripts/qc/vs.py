import pandas as pd
import numpy as np
import scanpy as sc
import os
import matplotlib.pyplot as plt
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--sample_path', required=True)
parser.add_argument('-g','--min_genes', required=True)
parser.add_argument('-c','--min_cells', required=True)
parser.add_argument('-d','--colors_dict', required=True)
args = vars(parser.parse_args())

path_sample = args['sample_path']
min_genes = int(args['min_genes'])
min_cells = int(args['min_cells'])
colors_dict = args['colors_dict']
colors_dict = dict(item.split('_') for item in colors_dict.strip("'").split(';'))
sample_id = os.path.normpath(path_sample).split(os.path.sep)[-1]
plot_path = os.path.join('results', 'qc', 'vs_{0}.pdf'.format(sample_id))
out_path = os.path.join(path_sample, 'adata.h5ad')

# Read slide
adata = sc.read_visium(path_sample, count_file='raw_feature_bc_matrix.h5')
adata.var_names_make_unique()
adata.obsm['spatial'] = adata.obsm['spatial'].astype(float)

# Read areas
path_areas = os.path.join(path_sample, 'areas.csv')
adata.obs['areas'] = areas = pd.read_csv(path_areas, index_col=0)['Pathology_annotation']

# Remove spots not in tissue
msk = adata.obs['in_tissue'].astype(int) > 0
adata = adata[msk].copy()

# Basic filtering
sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_genes(adata, min_cells=min_cells)

# Compute QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Add colors
adata.uns['areas_colors'] = [colors_dict[c] for c in np.unique(adata.obs['areas'].dropna().values.astype('U'))]

# Plot
fig, axes = plt.subplots(1, 4, figsize=(12, 2.5), tight_layout=True, dpi=150)
axes = axes.ravel()

ax = axes[0]
sc.pl.spatial(adata, ax=ax, return_fig=False, show=False, frameon=False)
ax.set_title('H&E')

ax = axes[1]
adata.obs['log_total_counts'] = np.log10(adata.obs['total_counts'])
sc.pl.spatial(adata, color='log_total_counts', ax=ax, return_fig=False, show=False, frameon=False, size=1.5)
del adata.obs['log_total_counts']
ax.set_title('Number of UMIs (log10)')

ax = axes[2]
adata.obs['log_n_genes_by_counts'] = np.log10(adata.obs['n_genes_by_counts'])
sc.pl.spatial(adata, color='log_n_genes_by_counts', ax=ax, return_fig=False, show=False, frameon=False, size=1.5)
del adata.obs['log_n_genes_by_counts']
ax.set_title('Number of genes (log10)')

ax = axes[3]
sc.pl.spatial(adata, color='areas', ax=ax, return_fig=False, show=False, frameon=False, size=1.5)
ax.set_title('Areas')

fig.suptitle(sample_id)
fig.savefig(plot_path, bbox_inches='tight')

# Delet extra atribs
del adata.obs['in_tissue']
del adata.var
del adata.uns['areas_colors']

# Store counts and normalize
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# Save
adata.write(out_path)
