<<<<<<< HEAD

# python scripts/process/integrate.py --output cellbender
# python scripts/process/integrate.py --output cellranger

=======
>>>>>>> 4059a4c8e2af9c39af787ebee1439fc854d311d6
import scanpy as sc
import scanpy.external as sce
import numpy as np
import pandas as pd

import argparse
<<<<<<< HEAD
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
=======
import os

"""
Script to plot different QC metrics after filtering the data.
"""

# Read command line and set args
parser = argparse.ArgumentParser(prog='qc', description='Run QC per sample')
parser.add_argument('-i', '--input_path', help='Input path to merged object', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
args = vars(parser.parse_args())

input_path = args['input_path']
output_path = args['output_dir']
###############################
>>>>>>> 4059a4c8e2af9c39af787ebee1439fc854d311d6

# Read merged object
adata = sc.read_h5ad(input_path)

# Run harmony
sce.pp.harmony_integrate(adata, 'batch', adjusted_basis='X_pca', max_iter_harmony=30)
sc.pp.neighbors(adata)

# Run umap with updated connectivity
sc.tl.umap(adata)

# Write to file
<<<<<<< HEAD
adata.write(output_dir / out_name)
=======
adata.write(os.path.join(output_path, 'integrated.h5ad'))
>>>>>>> 4059a4c8e2af9c39af787ebee1439fc854d311d6
