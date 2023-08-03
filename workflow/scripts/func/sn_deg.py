import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import decoupler as dc
import seaborn as sns
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--inp_path', required=True)
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

inp_path = args['inp_path']
out_path = args['out_path']

# Read
adata = sc.read_h5ad(inp_path)

# Get filtered pseudo-bulk profile
pdata = dc.get_pseudobulk(
    adata,
    sample_col='Sample id',
    groups_col='leiden',
    layer='counts',
    mode='sum',
    min_cells=10,
    min_counts=1000
)

def do_contrast(cdata, contrast):
    try:
        dds = DeseqDataSet(
            adata=cdata,
            design_factors=['Lesion type','Sex'],
            refit_cooks=True,
            n_cpus=64
        )
        # Compute LFCs
        dds.deseq2()
        # Extract contrast between one condition vs rest
        stat_res = DeseqStats(
            dds,
            contrast=["Lesion type"] + contrast,
            n_cpus=64,
            independent_filter=False
        )
        # Compute Wald test
        stat_res.summary()
        # Extract results
        results_df = stat_res.results_df

        return results_df
    except:
        return None

def try_lesion_contrast(cdata, df, cell_type, sample_count, lesions):
    sub_count = sample_counts.loc[(sample_counts['leiden'] == cell_type) & (np.isin(sample_counts['Lesion type'], lesions)), 'count']
    if (sub_count.min() > 2) and (sub_count.size > 1):
        c_df = do_contrast(cdata, contrast=lesions)
        c_df['contrast'] = '{0}vs{1}'.format(lesions[0], lesions[1])
        c_df['cell_type'] = cell_type
        df.append(c_df)


cell_types = adata.obs['leiden'].unique()
df = []
sample_counts = pdata.obs.groupby(["Lesion type", "leiden"]).size().reset_index(name='count')

for cell_type in cell_types:
    cdata = pdata[pdata.obs['leiden'] == cell_type].copy()
    genes = dc.filter_by_expr(cdata, group='Lesion type', min_count=10, min_total_count=15)
    print(cell_type, genes.size)

    cdata = cdata[:, genes].copy()

    try_lesion_contrast(cdata, df, cell_type, sample_counts, lesions=["CA", "Ctrl"])
    try_lesion_contrast(cdata, df, cell_type, sample_counts, lesions=["CI", "Ctrl"])
    try_lesion_contrast(cdata, df, cell_type, sample_counts, lesions=["CA", "CI"])
        
df = pd.concat(df)
df.to_csv(out_path)

