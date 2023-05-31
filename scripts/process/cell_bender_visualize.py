
import scanpy as sc
import anndata as ad

import numpy as np
import pandas as pd

import os
import re
from pathlib import Path

# tmp
# scp -r philipp@172.27.41.194:/home/philipp/Work/VisiumMS/out/* /Users/pschafer/Python_projects/VisiumMS/out

# gene list for dotplot
gene_list = np.array(["GFAP", "ADCY2", "AQP4", "CFAP299", "CFAP43", "SPAG17", "IGHG1", "IGKC", "FCRL5", "FLT1", "CLDN5", "VWF", "APOE", "CD74", 
                      "FTL", "LRMDA", "ARHGAP24", "TBXAS1", "PDGFRA", "PTPRZ1", "PCDH15", "MOG", "ST18", "MOBP", "BCAS1", "SIRT2", "GPR17", 
                      "LAMA2", "CEMIP", "ABCA9", "PARP8", "SKAP1", "PRKCH"])


current_folder = Path(__file__).parent
out_folder = current_folder / ".." / ".." / "out" / "dotplots"
out_folder.mkdir(parents=True, exist_ok=True)
sc._settings.ScanpyConfig.figdir = out_folder  # where scanpy saves plots
output_dir = current_folder / ".." / ".." / "data" / "cellbender_out"

samples = [sample for sample in os.listdir(output_dir) if not sample.startswith(".")]

# load the annoated cellbender output
adata_objects = {sample: sc.read_h5ad(output_dir / sample / "cell_bender_matrix_filtered_qc_annotated.h5ad") for sample in samples}

for sample in adata_objects.keys():
    adata = adata_objects[sample]
    adata.var_names_make_unique()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    genes = gene_list[np.isin(gene_list, adata.var_names)]
    sc.pl.dotplot(adata, var_names=genes, groupby="cell_type", title=sample, save=sample)

adata_global = sc.concat(list(adata_objects.values()), join="outer", label="sample", keys=list(adata_objects.keys()))
sc.pl.dotplot(adata, var_names=gene_list, groupby="cell_type", swap_axes=True, figsize=(6, 12), save="global")
