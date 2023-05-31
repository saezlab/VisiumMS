
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from pathlib import Path

# this line forces theano to use the GPU and should go before importing cell2location
os.environ["THEANO_FLAGS"] = 'device=cuda0,floatX=float32,force_device=True'

import cell2location

from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel

# TODO: harcoded config
sample_id = "sample_id"
label_name = "cell_type"
labels_to_remove = ["unannotated"]

# current_folder = globals()['_dh'][0]
current_folder = Path(__file__).parent
output_dir = current_folder / ".." / ".." / "data" / "cellbender_out"
model_out = current_folder / ".." / ".." / "data" / "c2l_models"
model_out.mkdir(parents=True, exist_ok=True)

samples = [sample for sample in os.listdir(output_dir) if not sample.startswith(".")]

# load the cellbender output
adata_objects = {sample: sc.read_h5ad(output_dir / sample / "cell_bender_matrix_filtered_qc_annotated.h5ad") for sample in samples}
adata_raw = sc.concat(list(adata_objects.values()), join="outer", label=sample_id, keys=list(adata_objects.keys()), index_unique="_")
adata_raw.var_names_make_unique()

# remove unannotated cells
adata_raw = adata_raw[~adata_raw.obs[label_name].isin(labels_to_remove), :]

# split MS and control samples
sample_meta = pd.read_excel(current_folder / ".." / ".." / "data" / "Metadata_all.xlsx")

#ms_samples = sample_meta["Brain bank ID"][sample_meta.Condition=="MS"]
#ctrl_samles = sample_meta["Brain bank ID"][sample_meta.Condition=="Control"]

ms_samples = [sample for sample in adata_raw.obs[sample_id].unique() if sample.startswith("MS")]
ctrl_samles = [sample for sample in adata_raw.obs[sample_id].unique() if sample not in ms_samples]

ms_adata_raw = adata_raw[adata_raw.obs[sample_id].isin(ms_samples), :].copy()
# print size of MS dataset
print(f"MS dataset:\n{ms_adata_raw}")
print(f"MS samples:\n{ms_adata_raw.obs[sample_id].unique()}")
ctrl_adata_raw = adata_raw[adata_raw.obs[sample_id].isin(ctrl_samles), :].copy()
print(f"Control dataset:\n{ms_adata_raw}")
print(f"Control samples:\n{ctrl_adata_raw.obs[sample_id].unique()}")
del(adata_raw)

# NOTE: I previously though that the cellbender output was floats
#ms_adata_raw.X = np.round(ms_adata_raw.X).astype(int)
#ctrl_adata_raw.X = np.round(ctrl_adata_raw.X).astype(int)

# Run one model for MS and for healthy controls
for adata, identifier in zip([ms_adata_raw, ctrl_adata_raw], ["MS", "Control"]):

    tmp_out = model_out / (identifier + "_reg_model")
    tmp_out.mkdir(parents=True, exist_ok=True)
    print(f"Running regression model for {identifier}, saving in {tmp_out}")

    # Filter by cell2loc thresholds
    selected = filter_genes(adata, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
    adata = adata[:, selected].copy()

    # use integer encdoing for sample and celltype covariates (scvi utility)
    cell2location.models.RegressionModel.setup_anndata(adata=adata,
                                                       # 10X reaction / sample / batch
                                                       batch_key=sample_id,
                                                       # cell type, covariate used for constructing signatures
                                                       labels_key=label_name
    )

    # Run regression model
    # See https://github.com/BayraktarLab/cell2location/blob/a583a836b3a932ac6b4de54edd56b8dcf235245a/cell2location/models/reference/_reference_module.py#L13
    mod = RegressionModel(adata)
    mod.view_anndata_setup()

    # Training 
    mod.train(max_epochs=250, batch_size=2500, train_size=1, lr=0.002, use_gpu=True)

    # Save training plot
    fig, ax = plt.subplots(1,1, facecolor='white')
    mod.plot_history(20, ax=ax)
    fig.savefig(tmp_out / "training_plot.png", dpi=300, bbox_inches='tight')

    # In this section, we export the estimated cell abundance (summary of the posterior distribution).
    adata = mod.export_posterior(
        adata, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
    )

    # export estimated expression in each cluster
    if 'means_per_cluster_mu_fg' in adata.varm.keys():
        inf_aver = adata.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                        for i in adata.uns['mod']['factor_names']]].copy()
    else:
        inf_aver = adata.var[[f'means_per_cluster_mu_fg_{i}'
                                        for i in adata_raw.uns['mod']['factor_names']]].copy()
    inf_aver.columns = adata.uns['mod']['factor_names']

    inf_aver.to_csv(tmp_out / "inf_aver.csv")
