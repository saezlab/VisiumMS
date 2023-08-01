
# python scripts/process/integrate.py --output cellbender
# python scripts/process/integrate.py --output cellranger

import scanpy as sc
import scanpy.external as sce
import numpy as np
import pandas as pd

import argparse
from pathlib import Path
import os

# add command line flag arguments to specify either "cellbender" or "cellranger" output
parser = argparse.ArgumentParser()
parser.add_argument("--output", type=str, required=True)
args = parser.parse_args()

# set up relative paths within the project
current_folder = Path(__file__).parent
output_dir = current_folder / ".." / ".." / "data" / "prc" / "sc"
if args.output == "cellbender":
    input_path = current_folder / ".." / ".." / "data" / "prc" / "sc" / "cellbender_merged.h5ad"
    out_name = "cellbender_integrated.h5ad"
elif args.output == "cellranger":
    input_path = current_folder / ".." / ".." / "data" / "prc" / "sc" / "cellranger_merged.h5ad"
    out_name = "cellranger_integrated.h5ad"
else:
    raise ValueError("output must be either 'cellbender' or 'cellranger'")

# Read merged object
adata = sc.read_h5ad(input_path)

# Run harmony
sce.pp.harmony_integrate(adata, 'batch', adjusted_basis='X_pca', max_iter_harmony=30)
sc.pp.neighbors(adata)

# Run umap with updated connectivity
sc.tl.umap(adata)

# Write to file
adata.write(output_dir / out_name)
