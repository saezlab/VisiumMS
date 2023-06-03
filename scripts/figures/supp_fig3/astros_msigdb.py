import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import os
import decoupler as dc

# Defina path
fig_path = 'figures/manuscript/supp_fig3/'
fig_name = 'astros_msigdb.pdf'
plt.rcParams['font.sans-serif'] = 'Arial'

# Read
adata = sc.read_h5ad('data/prc/sc/astros.h5ad')
adata = adata[~np.isin(adata.obs['leiden'], ['3', '6', '7'])]
adata.uns['log1p']["base"] = None

# Load resource
msigdb = dc.get_resource('MSigDB')

# Gsets to test
genesets = ['REACTOME_NR1H2_NR1H3_REGULATE_GENE_EXPRESSION_LINKED_TO_LIPOGENESIS',
'REACTOME_PLASMA_LIPOPROTEIN_ASSEMBLY',
'REACTOME_PLASMA_LIPOPROTEIN_ASSEMBLY_REMODELING_AND_CLEARANCE',
'REACTOME_PLASMA_LIPOPROTEIN_CLEARANCE',
'REACTOME_PLASMA_LIPOPROTEIN_REMODELING',
'REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR',
'REACTOME_CILIUM_ASSEMBLY']

# Filter
netdb = msigdb[np.isin(msigdb.geneset, genesets)]
netdb = netdb[~netdb.duplicated(['geneset', 'genesymbol'])]
netdb['geneset'] = [x.split('REACTOME_')[1] for x in netdb['geneset']]

# Estimate acts
dc.run_ulm(adata, netdb, source='geneset', target='genesymbol', weight=None)

# Plot
acts = dc.get_acts(adata, 'ulm_estimate')
fig = sc.pl.matrixplot(acts, acts.var_names, 'leiden', dendrogram=False, cmap='coolwarm', standard_scale='var', colorbar_title='Mean scaled activity', swap_axes=True, return_fig=True)
fig.savefig(os.path.join(fig_path, fig_name), bbox_inches='tight')
