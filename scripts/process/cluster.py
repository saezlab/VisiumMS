<<<<<<< HEAD

# python scripts/process/cluster.py --output cellbender --n_iters 20 --delta 0.05
# python scripts/process/cluster.py --output cellranger --n_iters 20 --delta 0.05

import scanpy as sc
import numpy as np
import pandas as pd
import argparse
from pathlib import Path
=======
import scanpy as sc
import numpy as np
import pandas as pd
>>>>>>> 4059a4c8e2af9c39af787ebee1439fc854d311d6

from sklearn.metrics import silhouette_score
from sklearn.metrics.pairwise import pairwise_distances

"""
Script to cluster using different resolutions to find the most optimal.
"""

<<<<<<< HEAD
# add command line flag arguments to specify either "cellbender" or "cellranger" output
parser = argparse.ArgumentParser()
parser.add_argument("--output", type=str, required=True)
parser.add_argument("--n_iters", type=str, required=True)
parser.add_argument("--delta", type=str, required=True)
args = vars(parser.parse_args())

# set up relative paths within the project
current_folder = Path(__file__).parent
output_dir = current_folder / ".." / ".." / "data" / "prc" / "sc"
if args['output'] == "cellbender":
    input_path = current_folder / ".." / ".." / "data" / "prc" / "sc" / "cellbender_integrated.h5ad"
elif args['output'] == "cellranger":
    input_dir = current_folder / ".." / ".." / "data" / "prc" / "sc" / "cellranger_integrated.h5ad"
else:
    raise ValueError("output must be either 'cellbender' or 'cellranger'")

n_iters = int(args['n_iters'])
delta = float(args['delta'])
if delta > 0.05:
    delta = 0.05

=======
# Read command line and set args
parser = argparse.ArgumentParser(prog='qc', description='Run QC per sample')
parser.add_argument('-i', '--input_path', help='Input path to integrated object', required=True)
parser.add_argument('-n', '--n_iters', help='Number of iterations', required=True)
args = vars(parser.parse_args())

input_path = args['input_path']
n_iters = int(args['n_iters'])
delta = int(args['delta'])
if delta > 0.05:
    delta = 0.05
>>>>>>> 4059a4c8e2af9c39af787ebee1439fc854d311d6
###############################

adata = sc.read_h5ad(input_path)

# Compute pair distances
distances = pairwise_distances(adata.obsm['X_pca'])

# Start with a low resolution
old_res = 0.05
sc.tl.leiden(adata, resolution=old_res, key_added='leiden_{0}'.format(old_res))
old_scr = silhouette_score(distances, np.array(adata.obs['leiden_{0}'.format(old_res)]), metric='precomputed')

# Store res and sscores
res_lst = [old_res]
scr_lst = [old_scr]

for i in range(n_iters):
    res = old_res + delta
<<<<<<< HEAD
    print('Iter {0}: {1} res'.format(i+1, res))
=======
    print('Iter {0}: {1} res'.format(i+1, res)
>>>>>>> 4059a4c8e2af9c39af787ebee1439fc854d311d6
    sc.tl.leiden(adata, resolution=old_res, key_added='leiden_{0}'.format(res))
    scr = silhouette_score(distances, np.array(adata.obs['leiden_{0}'.format(res)]), metric='precomputed')
    diff = scr - old_scr
    if diff > 0:
        delta = np.min(delta * (delta / diff), 0.05)
        res_lst.append(res)
        scr_lst.append(scr)
        old_res = res
        old_scr = scr
    else:
        break
    