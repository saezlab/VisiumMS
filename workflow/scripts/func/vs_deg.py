import scanpy as sc
import numpy as np
import pandas as pd
import os
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import decoupler as dc
import argparse
import anndata as ad


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-m','--meta_path', required=True)
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

meta_path = args['meta_path']
out_path = args['out_path']

# Read meta
meta = pd.read_csv(meta_path)
vs_samples = meta[~meta['Batch vs'].isnull()]['Sample id'].values.astype('U')

# Read
adata = []
for sample_id in vs_samples:
    print(sample_id)
    slide = sc.read_h5ad('data/prc/vs/{0}/adata.h5ad'.format(sample_id))
    slide.obs['niches'] = pd.read_csv('data/prc/vs/{0}/niches.csv'.format(sample_id), index_col=0)
    slide.obs['Sample id'] = sample_id
    slide.obs_names = [sample_id + '|' + i for i in slide.obs_names]
    adata.append(slide)
adata = ad.concat(adata, join='outer')

# Get filtered pseudo-bulk profile
pdata = dc.get_pseudobulk(
    adata,
    sample_col='Sample id',
    groups_col='niches',
    layer='counts',
    mode='sum',
    min_cells=10,
    min_counts=1000
)

def do_contrast(cdata, contrast):
    cdata = cdata[np.isin(cdata.obs['niches'], contrast)].copy()
    genes = dc.filter_by_expr(cdata, group='niches', min_count=10, min_total_count=15)
    print(contrast, genes.size)
    cdata = cdata[:, genes].copy()
    try:
        dds = DeseqDataSet(
            adata=cdata,
            design_factors=['niches','Sex'],
            refit_cooks=True,
            n_cpus=64
        )
        # Compute LFCs
        dds.deseq2()
        # Extract contrast between one niche vs another
        stat_res = DeseqStats(
            dds,
            contrast=["niches"] + contrast,
            n_cpus=64,
            cooks_filter=False,
            independent_filter=False
        )
        # Compute Wald test
        stat_res.summary()
        # Extract results
        results_df = stat_res.results_df

        results_df['contrast'] = '{0}vs{1}'.format(contrast[0], contrast[1])

        return results_df
    except:
        return None

contrasts = [['PPWM', 'WM'], ['LR', 'PPWM'], ['LC', 'LR'], ['VI', 'LC']]
df = []
for contrast in contrasts:
    tmp = do_contrast(pdata, contrast=contrast).sort_values('padj', ascending=False)
    df.append(tmp)
df = pd.concat(df)
df.to_csv(out_path)
