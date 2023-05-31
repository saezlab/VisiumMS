
import scanpy as sc
import pandas as pd
import numpy as np
import os
from pathlib import Path
import re

# this line forces theano to use the GPU and should go before importing cell2location
os.environ["THEANO_FLAGS"] = 'device=cuda0,floatX=float32,force_device=True'

import cell2location

# TODO: harcoded configs
n_cells_spot = 5
d_alpha = 20
recompute = False

current_folder = Path(__file__).parent
output_dir = current_folder / ".." / ".." / "data" / "cellbender_out"
#visium_dir = current_folder / ".." / ".." / "data" / "uscsc_dump"
visium_dir = current_folder / ".." / ".." / "data" / "raw" / "visium"
model_dir = current_folder / ".." / ".." / "data" / "c2l_models"
c2l_out = current_folder / ".." / ".." / "data" / "c2l_out"
c2l_out.mkdir(parents=True, exist_ok=True)

samples = [f for f in os.listdir(visium_dir) if not f.startswith(".")]
#samples = [f for f in os.listdir(visium_dir) if f.startswith("visium")]

for sample in samples:

    tmp_out = c2l_out / sample
    tmp_out.mkdir(parents=True, exist_ok=True)

    if (not recompute) and (tmp_out / "cell_abunds.csv").exists() and (tmp_out / "cell_props.csv").exists():
        print("Found existing results for sample", sample, "skipping")
        continue
    print("Running cell2location for sample", sample, "saving in", tmp_out)

    #base_name = re.sub(r"\.h5ad$", "", sample)
    #base_name = re.sub(r"^visium_", "", base_name)
    base_name = sample
    #adata_vis = sc.read_h5ad(visium_dir / sample)
    adata_vis = sc.read_visium(visium_dir / sample / "outs")
    adata_vis.var_names_make_unique()
    print(adata_vis.X.data)
    adata_vis.X.data = np.round(adata_vis.X.data).astype(int) 
    print(adata_vis.X.data)

    reg_path = model_dir / "MS_reg_model" if base_name.startswith("MS") else model_dir / "Control_reg_model"
    print(f"Using model {reg_path} for sample {sample}")

    inf_aver = pd.read_csv(reg_path / "inf_aver.csv", index_col=0)

    intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
    adata_vis = adata_vis[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()

    cell2location.models.Cell2location.setup_anndata(adata=adata_vis)

    # create and train the model
    mod = cell2location.models.Cell2location(
        adata_vis, 
        cell_state_df=inf_aver,  # marker gene expression averaged over cell types/states
        # the expected average cell abundance: tissue-dependent hyper-prior which can be estimated from paired histology:
        N_cells_per_location=n_cells_spot,
        # hyperparameter controlling normalisation of within-experiment variation in RNA detection:
        detection_alpha=d_alpha
    )

    mod.view_anndata_setup()

    # Train
    mod.train(max_epochs=30000,
              # train using full data (batch_size=None), why?
              #batch_size=None,
              # TODO: I will use a smaller batch size here to save GPU memory
              batch_size=2048,
              # use all data points in training because
              # we need to estimate cell abundance at all locations
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
    cell_abunds.to_csv(tmp_out / "cell_abunds.csv")
    cell_props.to_csv(tmp_out / "cell_props.csv")


