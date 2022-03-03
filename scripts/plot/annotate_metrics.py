import scanpy as sc
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import argparse
import os

"""
Script to plot different metrics after annotation.
"""

# Read command line and set args
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_path', help='Input path to object', required=True)
parser.add_argument('-m', '--markers_path', help='Input path to markers csv', required=True)
parser.add_argument('-v', '--version', help='Version of the annotation', required=False)
parser.add_argument('-o', '--output_dir', help='Output directory', required=True)
args = vars(parser.parse_args())

input_path = args['input_path']
markers_path = args['markers_path']
version = args['version']
output_path = args['output_dir']
###############################

# Read merged object
adata = sc.read_h5ad(input_path)

# Read marker genes
markers = pd.read_csv(markers_path)
markers = dict(markers.groupby('cell_type')['gene'].apply(list))

# Plot clusters and marker genes
fig, ax = plt.subplots(1,2, figsize=(14,7), facecolor='white', tight_layout=True, dpi=150, gridspec_kw={'width_ratios': [1, 1.5]})
sc.pl.umap(adata, color='leiden', ax=ax[0], return_fig=False, show=False, add_outline=False, s=10, legend_loc='on data')
sc.pl.dotplot(adata, markers, 'leiden', dendrogram=False, ax=ax[1], return_fig=False, show=False)

# Save
if version is not None:
    fname = 'annotate_summary_{0}.png'.format(version)
else:
    fname = 'annotate_summary.png'
fig.savefig(os.path.join(output_path, fname))