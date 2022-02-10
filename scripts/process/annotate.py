import scanpy as sc
import numpy as np
import pandas as pd

import argparse
import os

"""
Cluster cells and optionally annotate them.
"""

# Read command line and set args
parser = argparse.ArgumentParser(prog='ra', description='Reanotates atlas')
parser.add_argument('-i', '--input_path', help='Input path to integrated object', required=True)
parser.add_argument('-r', '--resolution', help='Resoulution for leiden clustering algorithm', required=True)
parser.add_argument('-d', '--dictionary', help='Dictionary of clusters and labels', required=False)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
args = vars(parser.parse_args())

input_path = args['input_path']
resolution = float(args['resolution'])
dictionary = args['dictionary']
output_path = args['output_dir']
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
adata.write(os.path.join(output_path, 'annotated.h5ad'))
