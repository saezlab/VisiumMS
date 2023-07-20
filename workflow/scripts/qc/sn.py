import scanpy as sc
import scanpy.external as sce
import numpy as np
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--sample_path', required=True)
parser.add_argument('-g','--min_genes', required=True)
parser.add_argument('-c','--min_cells', required=True)
parser.add_argument('-m','--mt_thr', required=True)
parser.add_argument('-n','--ct_thr', required=True)
parser.add_argument('-d','--db_thr', required=True)
args = vars(parser.parse_args())

path_sample = args['sample_path']
min_genes = int(args['min_genes'])
min_cells = int(args['min_cells'])
mt_thr = float(args['mt_thr'])
ct_thr = float(args['ct_thr'])
db_thr = float(args['db_thr'])
sample_id = os.path.normpath(path_sample).split(os.path.sep)[-2]
plot_path = os.path.join('results', 'qc', 'sn_{0}.pdf'.format(sample_id))
out_path = os.path.join(os.path.dirname(path_sample), 'adata.h5ad')

# Read
adata = sc.read_10x_h5(path_sample)
adata.var_names_make_unique()

# Basic filtering
sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_genes(adata, min_cells=min_cells)

# COmpute QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Compute doublet scores
sce.pp.scrublet(adata, verbose=False)

# Find msk
mt_msk = adata.obs['pct_counts_mt'] < mt_thr
ct_thr = np.quantile(adata.obs['n_genes_by_counts'], ct_thr)
ct_msk = adata.obs['n_genes_by_counts'] < ct_thr
db_msk = adata.obs['doublet_score'] < db_thr
msk = mt_msk & ct_msk & db_msk

# Count removed cells
n_mt = np.sum(~mt_msk)
n_ct = np.sum(~ct_msk)
n_db = np.sum(~db_msk)
n_total = np.sum(~msk)

# Plot
fig, ax = plt.subplots(2, 2, figsize=(7, 6), tight_layout=True, dpi=150)
ax = ax.ravel()

sns.histplot(
    data=adata.obs, x='total_counts', y='pct_counts_mt', cbar=True, cbar_kws=dict(shrink=.75),
    bins=100, log_scale=[True, False], ax=ax[0]
)
ax[0].set_xlabel('Total counts')
ax[0].set_ylabel('Fraction MT counts')
ax[0].axhline(y=mt_thr, linestyle='--', color='black')

sns.histplot(
    data=adata.obs, x='doublet_score', cbar=True, cbar_kws=dict(shrink=.75),
    bins=100, element="step", ax=ax[1]
)
ax[1].set_xlabel('Doublet score')
ax[1].set_ylabel('Number of cells')
ax[1].axvline(x=db_thr, linestyle='--', color='black')

sns.histplot(
    data=adata.obs, x='total_counts', y='n_genes_by_counts', cbar=True, cbar_kws=dict(shrink=.75),
    bins=100, log_scale=[True, True], ax=ax[2]
)
ax[2].set_xlabel('Total counts')
ax[2].set_ylabel('Number of genes expr')
ax[2].axhline(y=ct_thr, linestyle='--', color='black')


sns.barplot(
    x=['MT', 'Gene', 'Doublet', 'Total'], y=[n_mt, n_ct, n_db, n_total], ax=ax[3], color='#589cc4'
)
ax[3].set_ylabel('Number of cells lost')
fig.suptitle(sample_id)
fig.savefig(plot_path, bbox_inches='tight')

# Filter
adata = adata[msk, :].copy()
adata.write(out_path)
