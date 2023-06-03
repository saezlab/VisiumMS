import pandas as pd
import numpy as np
import scanpy as sc
import decoupler as dc
import os


# Load raw object
path_raw = 'data/prc/sc/raw.h5ad'
adata = sc.read(path_raw)

# Get pseudo-bulk profile
padata = dc.get_pseudobulk(adata, sample_col='sample_id', groups_col='cell_type',
                           min_prop=0.2, min_smpls=3, use_raw=False)

# Normalize
sc.pp.normalize_total(padata, target_sum=1e4)
sc.pp.log1p(padata)

# Run contrast
logFCs, pvals = dc.get_contrast(padata,
                                group_col='cell_type',
                                condition_col='lesion_type',
                                condition='Chronic Active',
                                reference='Control',
                                method='t-test'
                               )

# Extract deg
deg = dc.format_contrast_results(logFCs, pvals)

# Adjust pvals
adj_pvals = []
for contrast in np.unique(deg['contrast'].values):
    adj_pvals.extend(list(dc.p_adjust_fdr(deg[deg['contrast'] == contrast]['pvals'].values)))
deg['adj_pvals'] = adj_pvals

# Filter by basic thrs
deg = deg[(np.abs(deg['logFCs']) > 0.5) & (deg['pvals'] < 0.05)]

# Filter MT genes
deg = deg[[not g.startswith('MT-') for g in deg['name']]]

# Filter genes that are sign in other cell_types
counts = deg.groupby('name').count()
genes = counts.index[counts['contrast'] == 1].values
deg = deg[np.isin(deg['name'], genes)]

# Filter genes that have too much change
deg = deg[np.abs(deg['logFCs']) < 10]

# Split between pos and neg DEG
pos = deg[deg['logFCs'] > 0]
neg = deg[deg['logFCs'] < 0]

# Switch neg to pos values
neg = neg.assign(logFCs=lambda x: np.abs(x['logFCs']))

# Sort
pos = pos.sort_values(['contrast', 'pvals'])
neg = neg.sort_values(['contrast', 'pvals'])

# Write
res_path = 'data/prc/sign/deg'
os.makedirs(res_path, exist_ok=True)
pos.to_csv('{0}/pos.csv'.format(res_path), index=False)
neg.to_csv('{0}/neg.csv'.format(res_path), index=False)
deg.to_csv('{0}/deg.csv'.format(res_path), index=False)
