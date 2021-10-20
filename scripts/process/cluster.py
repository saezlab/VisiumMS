import scanpy as sc
import numpy as np
import pandas as pd

from sklearn.metrics import silhouette_score
from sklearn.metrics.pairwise import pairwise_distances

"""
Script to cluster using different resolutions to find the most optimal.
"""

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
    print('Iter {0}: {1} res'.format(i+1, res)
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
    