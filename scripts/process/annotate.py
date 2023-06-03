
# python scripts/process/annotate.py --output cellbender --resolution 1
# python scripts/process/annotate.py --output cellranger --resolution 1

import scanpy as sc
import numpy as np
import pandas as pd

import argparse
from pathlib import Path
import os

"""
Cluster cells and optionally annotate them.
"""

# Read command line and set args
parser = argparse.ArgumentParser(prog='ra', description='Reanotates atlas')
parser.add_argument("--output", type=str, required=True)
parser.add_argument('--resolution', help='Resoulution for leiden clustering algorithm', required=True)
parser.add_argument('--dictionary', help='Dictionary of clusters and labels', required=False)
args = vars(parser.parse_args())

# set up relative paths within the project
current_folder = Path(__file__).parent
output_dir = current_folder / ".." / ".." / "data" / "prc" / "sc"
if args['output'] == "cellbender":
    input_path = current_folder / ".." / ".." / "data" / "prc" / "sc" / "cellbender_integrated.h5ad"
    out_name = "cellbender_annotated.h5ad"
elif args.output == "cellranger":
    args['output'] = current_folder / ".." / ".." / "data" / "prc" / "sc" / "cellranger_integrated.h5ad"
    out_name = "cellranger_annotated.h5ad"
else:
    raise ValueError("output must be either 'cellbender' or 'cellranger'")

resolution = float(args['resolution'])
dictionary = args['dictionary']
###############################

# Read merged object
adata = sc.read_h5ad(input_path)

# Cluster cells
sc.tl.leiden(adata, resolution=resolution, key_added='leiden')

# Annotate if dict is available
if dictionary is not None:
    # Extract dict
    dictionary = {pair.split(':')[0]:pair.split(':')[1] for pair in dictionary.split(',')}
    
    # Update names
    adata.obs['leiden'] = [dictionary[clust] for clust in adata.obs['leiden']]
    
# Write to file
adata.write(output_dir / out_name)
