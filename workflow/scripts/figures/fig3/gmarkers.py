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

# AS genes
ctype = 'AS'
cdata = pdata[pdata.obs['leiden'] == ctype].copy()
genes = dc.filter_by_expr(cdata, group='Lesion type', min_count=10, min_total_count=15)
cdata = cdata[:, genes].copy()
sc.pp.normalize_total(cdata, target_sum=1e4)
sc.pp.log1p(cdata)
genes = [
    'PDK1','ACSF2','TANGO2','AQP1','ARL17B','SMAD6','HMGB1','HLA-F-AS1', 'LRP1', 'PRKG2','OSMR',
    'SERPINA3','NFE2L2', 'APP', 'C3', 'DST','TNFRSF1A', 'EEA1','DOCK11','RNF7','CCND2',
    'CNDP2','HSPBP1', 'FOXJ1','ANXA1','SLITRK2','ITGA2','SCRG1','CLU','ITGB1','CFAP97',
    'EBLN3P', 'COL6A1','CTSD','IFITM3','TAGLN','AKT3','CPNE3','KCNK2','PIK3IP1',
    'IFITM2', 'ARL5B', 'NNT', 'CUL4B', 'ALDH1L1'
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

# MG genes
ctype = 'MG'
cdata = pdata[pdata.obs['leiden'] == ctype].copy()
genes = dc.filter_by_expr(cdata, group='Lesion type', min_count=10, min_total_count=15)
cdata = cdata[:, genes].copy()
sc.pp.normalize_total(cdata, target_sum=1e4)
sc.pp.log1p(cdata)
genes = [
    'APPL2','PRAM1','SNHG5','ARHGAP12','RGS10','TMEM119',
    'HAMP','P2RY12','HMOX1', 'NOTCH1','SERPINB6','FAM20C','TNS3',
    'SIRPA','SPP1','SIGLEC11','TRAF3','CCDC40', 'FCGBP', 'PARP9',
    'CLN6','RBM47', 'FLT1','KATNAL2','SND1','TLN1','RALGAPA2','RAF1',
    'LACC1','PDIA4','PLIN2','PARVG','RHOBTB3','MS4A6A','CD14','PPARG',
    'SURF1','COPB2','C1QB','APOC1','ARMC2','SPCS1','FTL','C1QC',
    'GPNMB', 'APP','VPS4A','AHI1','FPR3','CD163','SERF2','EIF4G2',
    'CTSD','ZC3H12C','IQGAP1','KIF13A','S100A11','PLXNC1','APOE','SPCS3',
    'ADRM1','CD68', 'ITGA4','CLIP4', 'LRP10','ANXA2','TPP1','LGALS1', 
    'GPX3','TGFBI','MS4A7','ALDH1A1','EIF3F','ITGB1','PDCD10','PIK3CG', 
    'AGPS','LRRK2','PBX3','CLEC12A','ROR1', 'CLECL1','HAVCR1','SLC4A7',
    'CLEC7A','TNFRSF13C','FRMD4A','IGF1','CA8','CX3CR1','FOXP2','MYO1C'
]
genes = [
    'PRAM1','ARHGAP12','TMEM119','HAMP','P2RY12',
    'HMOX1','SIRPA','SPP1','TRAF3','PARP9','FLT1','PARVG','MS4A6A','CD14','PPARG','COPB2','C1QB',
    'APOC1','FTL','C1QC','GPNMB','APP','AHI1','FPR3','CD163','CTSS','IQGAP1',
    'PLXNC1','APOE','CD68','ITGA4','ANXA2','LGALS1','TGFBI','MS4A7','ITGB1',
    'PIK3CG','PDCD10','AGPS','LRRK2','CLEC12A','CLEC7A','FRMD4A','CX3CR1', 'FOXP2'
]

# Plot
fig2 = sc.pl.matrixplot(
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
fig2.show()
fig2 = fig2.fig

# OL genes
ctype = 'OL'
cdata = pdata[pdata.obs['leiden'] == ctype].copy()
genes = dc.filter_by_expr(cdata, group='Lesion type', min_count=10, min_total_count=15)
cdata = cdata[:, genes].copy()
sc.pp.normalize_total(cdata, target_sum=1e4)
sc.pp.log1p(cdata)
genes = [
    'ACKR2','CAMK2A','ERBB2','CD226','PGM1','FAM13C','HMX1','SOX13','CDH1','TCFL5','ADAMTS4',
    'CHST6','BDH1','PSEN1','TRIM2','DNM2','DOCK3','NDE1','TRAF3','EPHB1','ANP32A','CPEB1', 
    'DSCAML1','CDK18','CHKA','CORO2B','DNMT1','AXIN1','E2F3','CPXM2','AHRR','COX19','CD82', 
    'COX10', 'SRA1','CANT1','BIN3','HSPB1','LRRC8A','HVCN1','IRF1','IL17RB','EIF5','PDIA5',
    'NFKB2','CALM1','THAP1','TSR2','HSP90B1','PGM3','CAMK2D','HCN3','USP1','IGHGP','CDCA8',
    'NGFR','GPRC5A','ATF4','BAG6','CD274', 'TGFBR3','PCBP1','SLC22A17','EIF2S1','ANGPT2',
    'BMI1','RAB3D','DCC','MAF1','CTSD','AUP1','BRCA2','AHR','ITGB1','RHEB','PRRC2A','LRRN3',
    'SOX4','LGALS3','ANXA2','AMD1','LPL','MYRIP','MPZ','GMFB','H1FX','CLGN','ROBO2','NECTIN3',
    'RNF141','CNTNAP3','CITED4','TGFBR2','VWA8','COL28A1','CD2AP','ARNTL','PHACTR2','MAGI1',
    'ELOVL6','SVEP1','RORA','TMEFF1','SLC9B1','SCD'
]
genes = [
    'CAMK2A','ERBB2','CD226','PGM1','SOX13','ADAMTS4','PSEN1','DOCK3','NDE1','EPHB1',
    'CDK18','AXIN1','SRA1','HSPB1','IRF1','EIF5','NFKB2','CALM1','THAP1','HSP90B1',
    'CAMK2D','USP1','NGFR','ATF4','CD274','SLC22A17','EIF2S1','DCC','AUP1','BRCA2',
    'ITGB1','SOX4','LGALS3','ANXA2','MYRIP','MPZ','GMFB','TGFBR2','VWA8','CD2AP',
    'ARNTL','ELOVL6','SVEP1','RORA','SCD'
]
genes = [
    'CAMK2A','ERBB2','CD226','PGM1','SOX13','ADAMTS4','PSEN1','DOCK3','NDE1','EPHB1',
     'CDK18','AXIN1','SRA1','HSPB1','IRF1','EIF5','NFKB2','CALM1','THAP1','HSP90B1',
     'CAMK2D','USP1','NGFR','ATF4','CD274','SLC22A17','EIF2S1','DCC','BRCA2',
     'ITGB1','SOX4','LGALS3','OSMR','ANXA2','MYRIP','MPZ','GMFB','TGFBR2','VWA8','CD2AP',
     'ARNTL','ELOVL6','SVEP1','RORA','SCD'
]

# Plot
fig3 = sc.pl.matrixplot(
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
fig3.show()
fig3 = fig3.fig

# Save to pdf
pdf = matplotlib.backends.backend_pdf.PdfPages(plot_path)
for fig in [fig1, fig2, fig3]:
    pdf.savefig(fig, bbox_inches='tight')
pdf.close()
