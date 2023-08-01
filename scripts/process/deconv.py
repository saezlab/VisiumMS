<<<<<<< HEAD

# usage examples
# python scripts/process/cell_bender_c2l.py --output cellbender --model all --sample MS466
# python scripts/process/cell_bender_c2l.py --output cellbender --model condition --sample MS466
# python scripts/process/cell_bender_c2l.py --output cellbender --model lesion_type --sample MS466
# python scripts/process/cell_bender_c2l.py --output cellranger --model all --sample MS466
# python scripts/process/cell_bender_c2l.py --output cellranger --model condition --sample MS466
# python scripts/process/cell_bender_c2l.py --output cellranger --model lesion_type --sample MS466

import scanpy as sc
import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path
import argparse
=======
"""
Script to deconvolute visium slides into cell types using a pre-trained regression model.
"""

import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
>>>>>>> 4059a4c8e2af9c39af787ebee1439fc854d311d6

# this line forces theano to use the GPU and should go before importing cell2location
os.environ["THEANO_FLAGS"] = 'device=cuda0,floatX=float32,force_device=True'

import cell2location
<<<<<<< HEAD

# get cmd line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--output", type=str, required=True, help="['cellbender', 'cellranger']")
parser.add_argument("--model", type=str, required=True, help="['all', 'condition', 'lesion_type']")
parser.add_argument("--sample", type=str, required=True, help="sample name")
parser.add_argument("--n_cells_spot", type=int, required=False, default=5, help="prior for number of cells per spot")
parser.add_argument("--d_alpha", type=int, required=False, default=20, help="prior for heterogeneity?")
parser.add_argument("--recompute", type=str, required=False, default="False", help="prior for heterogeneity?")
args = parser.parse_args()

# print all the arguments for debugging
print("Arguments:")
for arg in vars(args):
    print(arg, getattr(args, arg))

# check the arguments
if args.output not in ['cellbender', 'cellranger']:
    raise ValueError("Model must be in ['cellbender', 'cellranger']'")
if args.model not in ['all', 'condition', 'lesion_type']:
    raise ValueError("Model must be in ['all', 'condition', 'lesion_type']")
if args.recompute not in ["True", "true", "False", "false"]:
    raise ValueError("Recompute must be in ['True', 'true', 'False', 'false']")

sample = args.sample
n_cells_spot = args.n_cells_spot
d_alpha = args.d_alpha
recompute = args.recompute in ["True", "true"]

# set the paths
current_folder = Path(__file__).parent
visium_dir = current_folder / ".." / ".." / "data" / "raw" / "vis"
model_dir = current_folder / ".." / ".." / "data" / "prc" / "sc" / "c2l_model" / args.output

# get the correct regression model
sample_meta = pd.read_excel(current_folder / ".." / ".." / "data" / "Metadata_all.xlsx", sheet_name="Visium")
sample_entry = sample_meta.loc[sample_meta.sample_id == sample, :].to_dict(orient="records")[0]
print(sample_entry)

if args.model == "all":
    reg_model = "All"
elif args.model == "condition":
    if sample_entry["Condition"] == "Control":
        reg_model = "Control"
    elif sample_entry["Condition"] == "MS":
        reg_model = "MS"
    else:
        raise ValueError("Unknown condition")
elif args.model == "lesion_type":
    if sample_entry["lesion_type"] == "Ctrl":
        reg_model = "Control"
    elif sample_entry["lesion_type"] == "CI":
        reg_model = "CI"
    elif sample_entry["lesion_type"] == "CA":
        reg_model = "CA"
    elif sample_entry["lesion_type"] == "A":
        reg_model = "A"
    else:
        raise ValueError("Unknown lesion type")
else:
    raise ValueError("Unknown model")
reg_path = model_dir / (reg_model + "_reg_model")
inf_aver = pd.read_csv(reg_path / "inf_aver.csv", index_col=0)
print("Using model", reg_model, "for sample", sample, "and model specification", args.model)

c2l_out = current_folder / ".." / ".." / "data" / "prc" / "vis" / "c2l_out" / args.output
c2l_out.mkdir(parents=True, exist_ok=True)

tmp_out = c2l_out / sample
tmp_out.mkdir(parents=True, exist_ok=True)

# TODO: Adjust to new save format
if (not recompute) and (tmp_out / f"cell_abunds_{reg_model}.csv").exists() and (tmp_out / f"cell_props_{reg_model}.csv").exists():
    print("Found existing results for sample", sample, "skipping")
    sys.exit()
print("Running cell2location for sample", sample, "saving in", tmp_out)

adata_vis = sc.read_visium(visium_dir / sample / "outs")
adata_vis.var_names_make_unique()
print(adata_vis.X.data)  # do we have to round something here?

=======
import scvi


import argparse


# Read command line and set args
parser = argparse.ArgumentParser(prog='qc', description='Generates regression model from raw single nuc data')
parser.add_argument('-s', '--slide_path', help='Path to slide adata', required=True)
parser.add_argument('-r', '--reg_path', help='Path to regression model', required=True)
parser.add_argument('-n', '--n_cells_spot', help='Number of cells per spot', required=False, type=float, default=5)
parser.add_argument('-a', '--d_alpha', help='Regularization parammeter for technical effects', required=False, type=float, default=20)
parser.add_argument('-o', '--path_output', help='Path were to save deconvolution', required=True)
args = vars(parser.parse_args())

slide_path = args['slide_path']
reg_path = args['reg_path']
n_cells_spot = int(args['n_cells_spot'])
d_alpha = int(args['d_alpha'])
path_output = args['path_output']

# Read inputs
adata_vis = sc.read_visium(os.path.join(slide_path, 'outs'))
adata_vis.var_names_make_unique()
inf_aver = pd.read_csv(reg_path, index_col=0)

# find shared genes and subset both anndata and reference signatures
>>>>>>> 4059a4c8e2af9c39af787ebee1439fc854d311d6
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

<<<<<<< HEAD
=======
# prepare anndata for cell2location model
>>>>>>> 4059a4c8e2af9c39af787ebee1439fc854d311d6
cell2location.models.Cell2location.setup_anndata(adata=adata_vis)

# create and train the model
mod = cell2location.models.Cell2location(
<<<<<<< HEAD
    adata_vis, 
    cell_state_df=inf_aver,  # marker gene expression averaged over cell types/states
    # the expected average cell abundance: tissue-dependent hyper-prior which can be estimated from paired histology:
    N_cells_per_location=n_cells_spot,
    # hyperparameter controlling normalisation of within-experiment variation in RNA detection:
    detection_alpha=d_alpha
)

=======
    adata_vis, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=n_cells_spot,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=d_alpha
)
>>>>>>> 4059a4c8e2af9c39af787ebee1439fc854d311d6
mod.view_anndata_setup()

# Train
mod.train(max_epochs=30000,
<<<<<<< HEAD
            # train using full data (batch_size=None), why?
            batch_size=None,
            # use all data points in training because we need to estimate cell abundance at all locations
            train_size=1,
            use_gpu=True)
=======
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
          use_gpu=True)
>>>>>>> 4059a4c8e2af9c39af787ebee1439fc854d311d6

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
)

<<<<<<< HEAD
# Save the model for future inspection
mod.save(str(tmp_out / f"c2l_mod_{reg_model}"), overwrite=True)
adata_vis.write(tmp_out / f"sp_{reg_model}.h5ad")

# Save the results for 5% quantile, mean, and 95% quantile
# NOTE: Why do we take the 5% quantile of the posterior?
# NOTE: We are not computing the proportion here anymore, because we can do that downstream
for point_estimates in ["q05_cell_abundance_w_sf", "means_cell_abundance_w_sf", "q95_cell_abundance_w_sf"]:
    cell_abunds = adata_vis.obsm[point_estimates].copy()
    cell_abunds.columns = adata_vis.uns['mod']['factor_names']
    cell_abunds.to_csv(tmp_out / f"cell_abunds_{reg_model}_{point_estimates}.csv")
=======
# Extract abundances df and rename cols to cell types
cell_abunds = adata_vis.obsm['q05_cell_abundance_w_sf'].copy()
cell_abunds.columns = adata_vis.uns['mod']['factor_names']

# Compute proportions
cell_props = cell_abunds / np.sum(cell_abunds, axis=1).values.reshape(-1, 1)

# Store results
os.makedirs(os.path.join(path_output), exist_ok=True)
cell_abunds.to_csv(os.path.join(path_output, 'cell_abunds.csv'))
cell_props.to_csv(os.path.join(path_output, 'cell_props.csv'))

>>>>>>> 4059a4c8e2af9c39af787ebee1439fc854d311d6
