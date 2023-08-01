<<<<<<< HEAD

# python scripts/process/merge.py --output cellbender
# python scripts/process/merge.py --output cellranger

=======
>>>>>>> 4059a4c8e2af9c39af787ebee1439fc854d311d6
import scanpy as sc
import scanpy.external as sce

import numpy as np
import pandas as pd

import os
import argparse
<<<<<<< HEAD
from pathlib import Path
=======
>>>>>>> 4059a4c8e2af9c39af787ebee1439fc854d311d6

'''
Open all samples QC processed files, concatenate them and run integration
'''

<<<<<<< HEAD
# add command line flag arguments to specify either "cellbender" or "cellranger" output
parser = argparse.ArgumentParser()
parser.add_argument("--output", type=str, required=True)
args = parser.parse_args()

# set up relative paths within the project
current_folder = Path(__file__).parent
output_dir = current_folder / ".." / ".." / "data" / "prc" / "sc"
if args.output == "cellbender":
    input_dir = current_folder / ".." / ".." / "data" / "prc" / "sc" / "cellbender_qc"
    out_name = "cellbender_merged.h5ad"
elif args.output == "cellranger":
    input_dir = current_folder / ".." / ".." / "data" / "prc" / "sc" / "cellranger_qc"
    out_name = "cellranger_merged.h5ad"
else:
    raise ValueError("output must be either 'cellbender' or 'cellranger'")

# verbose
print("input_dir: ", input_dir)
print("out_name: ", out_name)

# Load meta data
meta = pd.read_excel(current_folder / ".." / ".." / "data" / "Metadata_all.xlsx", sheet_name="snRNA-seq")
samples = np.unique(meta['sample_id'])

files = [f for f in os.listdir(input_dir) if not f.endswith('.h5ad')]
file_names = [f.split('.')[0] for f in files]

# check that all samples are present
for sample in samples:
    if sample not in file_names:
        raise ValueError("Sample {} not present in input directory".format(sample))

adata = []
for sample in samples:
    print(sample)
    
    # Read adata
    tmp = sc.read_h5ad(input_dir / (sample + '.h5ad'))
=======
# Read command line and set args
parser = argparse.ArgumentParser(prog='qc', description='Run QC per sample')
parser.add_argument('-i', '--input_dir', help='Input directory containing all sample directories', required=True)
parser.add_argument('-m', '--metadata', help='Path to metadata', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the processed object', required=True)
args = vars(parser.parse_args())

input_path = args['input_dir']
meta_path = args['metadata']
output_path = args['output_dir']
###############################

# Load meta data
meta = pd.read_csv(meta_path)
samples = np.unique(meta['sample_id'])

adata = []
for sample in os.listdir(input_path):
    if not os.path.isdir(os.path.join(input_path, sample)) or sample.startswith('.'):
        continue
    print(sample)
    
    # Read adata
    tmp = sc.read_h5ad(os.path.join(input_path, sample, sample+'.h5ad'))
>>>>>>> 4059a4c8e2af9c39af787ebee1439fc854d311d6
    
    # Fetch sample metadata
    m = meta[meta['sample_id'] == sample]
    
    # Add metadata to adata
    for col in m.columns:
        tmp.obs[col] = m[col].values[0]
    
    # Append
    adata.append(tmp)
    del tmp
    
# Merge objects and delete list
adata = adata[0].concatenate(adata[1:], join='outer')
<<<<<<< HEAD
print(adata)

# Log-normalize expression
sc.pp.normalize_total(adata, target_sum=1e4)
=======

# Log-normalize expression
sc.pp.normalize_total(adata)
>>>>>>> 4059a4c8e2af9c39af787ebee1439fc854d311d6
sc.pp.log1p(adata)

# Compute HVG
sc.pp.highly_variable_genes(adata, batch_key='batch')
adata.var = adata.var[['highly_variable','highly_variable_nbatches']]
adata.raw = adata

# Filter by HVG
num_hvg_genes = 3000
batch_msk = np.array(adata.var.highly_variable_nbatches > 1)
hvg = adata.var[batch_msk].sort_values('highly_variable_nbatches').tail(num_hvg_genes).index
adata.var['highly_variable'] = [g in hvg for g in adata.var.index]
adata.var = adata.var[['highly_variable','highly_variable_nbatches']]
adata = adata[:,hvg]

# Run PCA
sc.pp.scale(adata)
sc.tl.pca(adata, svd_solver='arpack')

# Run UMAP
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# Write to file
<<<<<<< HEAD
adata.write(output_dir / out_name)
=======
adata.write(os.path.join(output_path, 'merged.h5ad'))
>>>>>>> 4059a4c8e2af9c39af787ebee1439fc854d311d6
