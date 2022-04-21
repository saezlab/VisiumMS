import scanpy as sc
import scanpy.external as sce

import numpy as np
import pandas as pd

import pickle
import os
import argparse
import urllib.request

"""
Script to run basic QC filtering. Stores data for future plotting and writes new AnnData objects.
"""

# Read command line and set args
parser = argparse.ArgumentParser(prog='qc', description='Run QC per sample')
parser.add_argument('-i', '--input_dir', help='Input directory containing all sample directories', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the processed objects', required=True)
args = vars(parser.parse_args())

base_path = args['input_dir']
output_path = args['output_dir']
###############################

# Load diss genes
diss_fname = 'coregene_df-FALSE-v3.csv'
diss_path = os.path.join(base_path, diss_fname)

if not os.path.isfile(diss_path):
    url = 'https://raw.githubusercontent.com/kieranrcampbell/' \
    'scrnaseq-digestion-paper/master/data/deliverables/coregene_df-FALSE-v3.csv'
    urllib.request.urlretrieve(url, diss_path)
    
diss_genes = pd.read_csv(diss_path).sort_values('PValue').head(200).gene_symbol.tolist()
os.remove(diss_path)

# If needed set up specific thresholds for doublet detection
# The value should divide the bimodal distirbution in two
doublet_thresholds = {
}

# Perform QC for each sample independently
for sample in os.listdir(base_path):
    if sample.startswith('.'):
        continue
    print(sample)
    # Read raw data
    adata = sc.read_10x_mtx(os.path.join(base_path, sample, 'filtered_feature_bc_matrix'),
                            var_names='gene_symbols', cache=True)
    adata.var_names_make_unique()

    # Basic filtering
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    # Compute QC metrics
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    # Compute doublets score
    sce.pp.scrublet(adata, verbose=False)
    
    # Compute dissociation score and normalize it (0-1)
    intersect_genes = adata.var.index.intersection(diss_genes)
    if len(intersect_genes) > 1:
        sc.tl.score_genes(adata, gene_list=intersect_genes, ctrl_size=len(intersect_genes), 
                          score_name='diss_score')
        adata.obs['diss_score'] = (adata.obs.diss_score-np.min(adata.obs.diss_score)) / \
        (np.max(adata.obs.diss_score) - np.min(adata.obs.diss_score))
    else:
        adata.obs['diss_score'] = 0

    # Set filter values (can be changed)
    mt_thr = 5.0
    gene_qnt = 0.99
    diss_qnt = 0.99

    if sample in doublet_thresholds:
        doublet_thr = doublet_thresholds[sample]
    else:
        doublet_thr = 0.2

    # Save cell meta data
    df = adata.obs[['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'doublet_score', 'diss_score']]
    plot_data = {'mt_thr' : mt_thr,
                 'gene_qnt' : gene_qnt,
                 'doublet_thr' : doublet_thr,
                 'diss_qnt' : diss_qnt,
                 'df' : df
                }

    # Filter
    gene_thr = np.quantile(adata.obs.n_genes_by_counts, gene_qnt)
    diss_thr = np.quantile(adata.obs.diss_score, diss_qnt)
    msk = (adata.obs.n_genes_by_counts < gene_thr) & \
          (adata.obs.pct_counts_mt < mt_thr) &  \
          (adata.obs.doublet_score < doublet_thr) & \
          (adata.obs.diss_score < diss_thr)
    adata = adata[msk, :]

    # Save results
    os.makedirs(os.path.join(output_path, sample), exist_ok=True)
    pickle.dump(plot_data, open(os.path.join(output_path, sample, sample+'.pkl'), "wb"))
    adata.write(os.path.join(output_path, sample,sample+'.h5ad'))
