<<<<<<< HEAD

# python scripts/process/annotate.py --output cellbender --resolution 0.75 
# python scripts/process/annotate.py --output cellranger --resolution 0.75

import scanpy as sc

import argparse
from pathlib import Path
from annotate_dicts import annotation_dict
=======
import scanpy as sc
import numpy as np
import pandas as pd

import argparse
import os
>>>>>>> 4059a4c8e2af9c39af787ebee1439fc854d311d6

"""
Cluster cells and optionally annotate them.
"""

# Read command line and set args
parser = argparse.ArgumentParser(prog='ra', description='Reanotates atlas')
<<<<<<< HEAD
parser.add_argument("--output", type=str, required=True)
parser.add_argument('--resolution', help='Resoulution for leiden clustering algorithm', required=True)
args = vars(parser.parse_args())

# set up relative paths within the project
current_folder = Path(__file__).parent
output_dir = current_folder / ".." / ".." / "data" / "prc" / "sc"
if args['output'] == "cellbender":
    input_path = current_folder / ".." / ".." / "data" / "prc" / "sc" / "cellbender_integrated.h5ad"
    out_name = "cellbender_annotated.h5ad"
elif args['output'] == "cellranger":
    input_path = current_folder / ".." / ".." / "data" / "prc" / "sc" / "cellranger_integrated.h5ad"
    out_name = "cellranger_annotated.h5ad"
else:
    raise ValueError("output must be either 'cellbender' or 'cellranger'")

resolution = float(args['resolution'])
dictionary = annotation_dict.get(args['output']).get(args['resolution'])
=======
parser.add_argument('-i', '--input_path', help='Input path to integrated object', required=True)
parser.add_argument('-r', '--resolution', help='Resoulution for leiden clustering algorithm', required=True)
parser.add_argument('-d', '--dictionary', help='Dictionary of clusters and labels', required=False)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
args = vars(parser.parse_args())

input_path = args['input_path']
resolution = float(args['resolution'])
dictionary = args['dictionary']
output_path = args['output_dir']
>>>>>>> 4059a4c8e2af9c39af787ebee1439fc854d311d6
###############################

# Read merged object
adata = sc.read_h5ad(input_path)

# Cluster cells
sc.tl.leiden(adata, resolution=resolution, key_added='leiden')
<<<<<<< HEAD
print(adata.obs.leiden.value_counts())

# Annotate if dict is available
if dictionary is not None:
=======

# Annotate if dict is available
if dictionary is not None:
    # Extract dict
    dictionary = {pair.split(':')[0]:pair.split(':')[1] for pair in dictionary.split(',')}
>>>>>>> 4059a4c8e2af9c39af787ebee1439fc854d311d6
    
    # Update names
    adata.obs['leiden'] = [dictionary[clust] for clust in adata.obs['leiden']]
    
# Write to file
<<<<<<< HEAD
adata.write(output_dir / out_name)
=======
adata.write(os.path.join(output_path, 'annotated.h5ad'))
>>>>>>> 4059a4c8e2af9c39af787ebee1439fc854d311d6
