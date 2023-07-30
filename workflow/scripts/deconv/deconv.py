import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
# this line forces theano to use the GPU and should go before importing cell2location
os.environ["THEANO_FLAGS"] = 'device=cuda0,floatX=float32,force_device=True'
import cell2location
import scvi
import argparse


# Read command line and set args
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--slide_path', required=True)
parser.add_argument('-r', '--reg_path', required=True)
parser.add_argument('-n', '--n_cells_spot', required=False, type=float, default=5)
parser.add_argument('-a', '--d_alpha', required=False, type=float, default=20)
parser.add_argument('-e', '--max_epochs', required=False)
parser.add_argument('-p', '--plot_path', required=False)
parser.add_argument('-o', '--out_path', required=True)
args = vars(parser.parse_args())

slide_path = args['slide_path']
reg_path = args['reg_path']
n_cells_spot = int(args['n_cells_spot'])
d_alpha = int(args['d_alpha'])
max_epochs = int(args['max_epochs'])
plot_path = args['plot_path']
out_path = args['out_path']
sample_id = os.path.basename(os.path.dirname(slide_path))

# Read inputs
adata_vis = sc.read_h5ad(slide_path)
del adata_vis.X
adata_vis.X = adata_vis.layers['counts'].copy()
del adata_vis.layers['counts']
inf_aver = pd.read_csv(reg_path, index_col=0)

# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis)

# create and train the model
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    N_cells_per_location=n_cells_spot,
    detection_alpha=d_alpha
)
mod.view_anndata_setup()

# Train
mod.train(max_epochs=max_epochs,
          batch_size=None,
          train_size=1,
          use_gpu=True)

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
)

# Extract abundances df and rename cols to cell types
cell_abunds = adata_vis.obsm['q05_cell_abundance_w_sf'].copy()
cell_abunds.columns = adata_vis.uns['mod']['factor_names']

# Compute proportions
cell_props = cell_abunds / np.sum(cell_abunds, axis=1).values.reshape(-1, 1)

# Store results
cell_abunds.to_csv(os.path.join(out_path, 'abunds.csv'))
cell_props.to_csv(os.path.join(out_path, 'props.csv'))

# Plot

mod.plot_history(1000)
plt.legend(labels=['full data training'])

print('Plot')
fig, ax = plt.subplots(1, 1, figsize=(5, 3), dpi=150)
mod.plot_history(iter_start=1000, ax=ax)
fig.suptitle('Deconv model for sample {0}'.format(sample_id))
fig.savefig(plot_path, bbox_inches='tight')

