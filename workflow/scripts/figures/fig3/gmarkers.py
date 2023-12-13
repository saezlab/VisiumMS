import pandas as pd
import numpy as np
import scanpy as sc
import decoupler as dc
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--ann_path', required=True)
parser.add_argument('-b','--deg_path', required=True)
parser.add_argument('-c','--plot_path', required=True)
args = vars(parser.parse_args())

ann_path = args['ann_path']
deg_path = args['deg_path']
plot_path = args['plot_path']

# Read
adata = sc.read_h5ad(ann_path)
deg = pd.read_csv(deg_path)

# Make psbulk (takes a while)
pdata = dc.get_pseudobulk(
    adata=adata,
    sample_col='Sample id',
    groups_col='leiden',
    layer='counts'
)

# Subset by cell type
ctype = 'AS'
cdata = pdata[pdata.obs['leiden'] == ctype].copy()
genes = dc.filter_by_expr(cdata, group='Lesion type', min_count=10, min_total_count=15)
cdata = cdata[:, genes].copy()
sc.pp.normalize_total(cdata, target_sum=1e4)
sc.pp.log1p(cdata)

# AS genes
genes = [
    'ACSF2', 'ALDH1L1','PDK1','TANGO2',
    'AQP1','ARL17B','SMAD6','HMGB1','HLA-F-AS1', 'LRP1', 'PRKG2',
    'NFE2L2','SERPINA3', 'APP', 'C3', 'TNFRSF1A','DST', 'EEA1','DOCK11','RNF7','CCND2', 
    'CNDP2','HSPBP1', 'FOXJ1','ANXA1','SLITRK2','ITGA2','SCRG1','CLU','ITGB1','CFAP97',
    'EBLN3P','PYGO1', 'COL6A1','CTSD','IFITM3',
    'TAGLN','AKT3','CPNE3','KCNK2','PIK3IP1','ARL5B','IFITM2','NNT','CUL4B'
]

# Plot
fig1 = sc.pl.matrixplot(
    cdata,
    genes,
    'Lesion type',
    dendrogram=False,
    standard_scale='var',
    colorbar_title='Z-scaled expression',
    cmap='viridis',
    categories_order=['Ctrl', 'CA', 'CI'],
    return_fig=True,
    show=False,
    title=ctype,
)
fig1.show()
fig1 = fig1.fig

# Save to pdf
pdf = matplotlib.backends.backend_pdf.PdfPages(plot_path)
for fig in [fig1]:
    pdf.savefig(fig, bbox_inches='tight')
pdf.close()
