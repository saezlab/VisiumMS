import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import os
import decoupler as dc

# Defina path
fig_path = 'figures/manuscript/supp_fig3/'
fig_name = 'microglia_msigdb.pdf'
plt.rcParams['font.sans-serif'] = 'Arial'

# Read
adata = sc.read_h5ad('data/prc/sc/microglia.h5ad')
adata = adata[~np.isin(adata.obs['leiden'], ['5'])]
adata.uns['log1p']["base"] = None

# Load resource
msigdb = dc.get_resource('MSigDB')

# Filter by hallmark
msigdb = msigdb[msigdb['collection']=='hallmark']

# Remove duplicated entries
msigdb = msigdb[~msigdb.duplicated(['geneset', 'genesymbol'])]
msigdb['geneset'] = [x.split('HALLMARK_')[1] for x in msigdb['geneset']]

# Filter relevant genesets
genesets = [
    'COMPLEMENT',
    'INFLAMMATORY_RESPONSE',
    'INTERFERON_ALPHA_RESPONSE',
    'INTERFERON_GAMMA_RESPONSE'
]
msigdb = msigdb[np.isin(msigdb['geneset'], genesets)]

# Estimate acts
dc.run_ulm(adata, msigdb, source='geneset', target='genesymbol', weight=None)

# Plot
acts = dc.get_acts(adata, 'ulm_estimate')
fig = sc.pl.matrixplot(acts, acts.var_names, 'leiden', dendrogram=False, cmap='coolwarm', standard_scale='var', colorbar_title='Mean scaled activity', swap_axes=True, return_fig=True)
fig.savefig(os.path.join(fig_path, fig_name), bbox_inches='tight')
