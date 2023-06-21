
# usage
# python scripts/process/run_mofa.py --model all --ct_metric abunds --image_features None --n_factors 10
# python scripts/process/run_mofa.py --model all --ct_metric props_ilr --image_features None --n_factors 10

# bash loop
#for model in all condition; do
#    for ct_metric in abunds props_ilr; do
#        for image_features in None histogram summary texture; do
#            for n_factors in 5 10; do
#                python scripts/process/run_mofa.py --model $model --ct_metric $ct_metric --image_features $image_features --n_factors $n_factors --recompute True
#            done
#        done
#    done
#done

# resources
# MOFApy: https://github.com/bioFAM/mofapy2/blob/master/mofapy2/notebooks/getting_started_python.ipynb
# mofax: https://github.com/bioFAM/mofax/blob/master/notebooks/getting_started_pbmc10k.ipynb

from mofapy2.run.entry_point import entry_point
import pandas as pd
import numpy as np
import os
import scanpy as sc
import mofax as mfx
from pathlib import Path

from umap import UMAP
import matplotlib.pyplot as plt
import seaborn as sns
import decoupler as dc
import plotnine as pn
import argparse

# get cmd line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--model", type=str, required=True, help="['all', 'condition', 'lesion_type']")
parser.add_argument("--ct_metric", type=str, required=True, help="['abunds', 'props', 'props_ilr']")
parser.add_argument("--image_features", type=str, required=True, help="['None', 'histogram', 'summary', 'texture']")
parser.add_argument("--n_factors", type=int, required=True)
parser.add_argument("--recompute", type=str, required=False, default="False")

# parse arguments and check
args = parser.parse_args()
assert args.model in ["all", "condition", "lesion_type"]
assert args.ct_metric in ["abunds", "props", "props_ilr"]
assert args.image_features in ["None", "histogram", "summary", "texture"]
assert args.recompute in ["True", "true", "False", "false"]
args.recompute = args.recompute in ["True", "true"]

# command line arguments, NOTE: hallmarks are always used
obsm_to_use = ["hallmark_estimates"]

if args.image_features != "None":
    obsm_to_use.append(args.image_features)

if args.ct_metric in ["abunds", "props"]:
    obsm_to_use.append(f"{args.ct_metric}_{args.model}")
elif args.ct_metric == "props_ilr":
    obsm_to_use.append(f"props_{args.model}_ilr")

print(f"Using the following views: {obsm_to_use}")

# initialise the entry point
ent = entry_point()

# set the paths
current_path = Path(__file__).parent
out_dir = current_path / ".." / ".." / "data" / "prc" / "vis" / "mofa_tests" / f"{args.model}__{args.ct_metric}__{args.image_features}__{str(args.n_factors)}"
plot_dir = out_dir / "plots"
sc.settings.figdir = str(plot_dir)
out_dir.mkdir(parents=True, exist_ok=True)
plot_dir.mkdir(parents=True, exist_ok=True) 
visium_path = current_path / ".." / ".." / "data" / "prc" / "vis" / "processed"

# check for recompute
if (out_dir / "mofa_model.hdf5").exists() and (not args.recompute):
    print(f"MOFA model for {out_dir} already exists, skipping...")
    exit()

# get visium samples
visium_samples = [f.split(".")[0] for f in os.listdir(visium_path) if not f.startswith(".")]
print(np.array(visium_samples))

# get the metadata
sample_meta = pd.read_excel(current_path / ".." / ".." / "data" / "Metadata_all.xlsx", sheet_name="Visium")

# load the visium data
vis_dict = {smp: sc.read_h5ad(visium_path / f"{smp}.h5ad") for smp in visium_samples}

# TODO: temporary fix to put the ilr transformed data into a dataframe
# NOTE: It is important that the column names do not resemble integers otherwise there are hd5 errors when saving the MOFA model
for values in vis_dict.values():
    values.obsm["props_all_ilr"] = pd.DataFrame(values.obsm["props_all_ilr"], index=values.obs.index, 
                                                columns=["irl_" + str(s) for s in list(range(values.obsm["props_all_ilr"].shape[1]))])
    values.obsm["props_condition_ilr"] = pd.DataFrame(values.obsm["props_condition_ilr"], index=values.obs.index, 
                                                      columns=["irl_" + str(s) for s in list(range(values.obsm["props_condition_ilr"].shape[1]))])
    values.obsm["props_lesion_type_ilr"] = pd.DataFrame(values.obsm["props_lesion_type_ilr"], index=values.obs.index, 
                                                        columns=["irl_" + str(s) for s in list(range(values.obsm["props_lesion_type_ilr"].shape[1]))])

# checks
print(visium_samples[0])
print(f"Adata object for first visium sample:\n{vis_dict[visium_samples[0]]}")
print(f"Available obsm layers:\n{vis_dict[visium_samples[0]].obsm_keys()}")
obsm_features = {obsm_key: vis_dict[visium_samples[0]].obsm[obsm_key].columns.to_list() for obsm_key in obsm_to_use}

# get cell metadata
meta_list = []
for sample in visium_samples:
    df = vis_dict[sample].obs.copy()
    df.index = [sample + "_" + s for s in df.index]
    df["sample_id"] = sample
    df["condition"] = sample_meta.loc[sample_meta.sample_id == sample, "Condition"].values[0]
    df["lesion_type"] = sample_meta.loc[sample_meta.sample_id == sample, "lesion_type"].values[0]
    meta_list.append(df)
meta_df = pd.concat(meta_list, axis=0)

# create mofa dataframe
df_list = []
for obsm_key in obsm_to_use:
    for sample in visium_samples:
        df = vis_dict[sample].obsm[obsm_key].copy()
        df.index = [sample + "_" + s for s in df.index] # unique barcodes are required!
        df = df.reset_index().melt(id_vars="index", var_name="feature", value_name="value")
        df = df.rename(columns={"index": "sample"})
        df["group"] = sample
        df["view"] = obsm_key
        df = df[["sample", "group", "feature", "value", "view"]]
        df_list.append(df)
data_dt = pd.concat(df_list)
print("Mofa dataframe:\n")
print(data_dt.head())
print(data_dt.tail())

# scale each view to unit variance
ent.set_data_options(
    scale_views = True
)

# set the likelihoods for all views
ent.set_data_df(data_dt, likelihoods = ["gaussian" for _ in range(len(obsm_to_use))])

# set the model options
ent.set_model_options(
    factors = args.n_factors,
    spikeslab_weights = True, 
    ard_weights = True
)

# set the training options
ent.set_train_options(
    convergence_mode = "fast", 
    dropR2 = 0.001, 
    gpu_mode = True, 
    seed = 1
)

ent.build()

ent.run()

ent.save(outfile=str(out_dir / "mofa_model.hdf5"))

# downstream part
model = mfx.mofa_model(out_dir / "mofa_model.hdf5")

# verbosity
print(f"""\
Cells: {model.shape[0]}
Features: {model.shape[1]}
Groups of cells: {', '.join(model.groups)}
Views: {', '.join(model.views)}
""")

# create adata object based on MOFA results
adata = sc.AnnData(X=model.get_factors())
adata.obs = model.get_cells()
adata.obs.set_index("cell", inplace=True)
sc.pp.neighbors(adata, n_neighbors=15, use_rep="X")
resolutions = [0.2, 0.3, 0.4, 0.5, 0.6]
for res in resolutions:
    sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res}")
adata.obs.rename(columns={"group": "sample_id"}, inplace=True)
adata.obs = adata.obs.join(sample_meta.set_index("sample_id"), on="sample_id")
sc.tl.umap(adata)

for obsm_key in obsm_to_use:
    df_list = []
    for sample in visium_samples:
        df = vis_dict[sample].obsm[obsm_key].copy()
        df.index = [sample + "_" + s for s in df.index] # unique barcodes are required!
        df_list.append(df)
    df = pd.concat(df_list)
    adata.obsm[obsm_key] = df.loc[adata.obs_names]

adata.write(out_dir / "adata.h5ad")

# check umaps
sc.pl.umap(adata, color=["sample_id", "Condition", "leiden_0.5"], ncols=1, save=".pdf", show=False)

# count the cluster fractions per sample, divide by the rowsums
for res in resolutions:
    df = adata.obs.groupby(["sample_id", f"leiden_{res}"]).size().unstack()
    df = df.div(df.sum(axis=1), axis=0)
    fig, ax = plt.subplots(figsize=(12, 8))
    sns.heatmap(df, cmap="viridis", annot=True, fmt=".2f", ax=ax)
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Sample ID")
    fig.tight_layout()
    fig.suptitle(f"Cluster fractions per sample (resolution {res})")
    fig.savefig(plot_dir / f"cluster_fractions_per_sample_heatmap_res_{res}.pdf")

# count the cluster fractions per lesion type, divide by the rowsums
for res in resolutions:
    df = adata.obs.groupby(["lesion_type", f"leiden_{res}"]).size().unstack()
    df = df.div(df.sum(axis=1), axis=0)
    fig, ax = plt.subplots(figsize=(12, 4))
    sns.heatmap(df, cmap="viridis", annot=True, fmt=".2f", ax=ax)
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Lesion Type")
    fig.tight_layout()
    fig.suptitle(f"Cluster fractions per lesion type (resolution {res})")
    fig.savefig(plot_dir / f"cluster_fractions_per_lesion_type_heatmap_res_{res}.pdf")

# count the cluster fractions per condition, divide by the rowsums
for res in resolutions:
    df = adata.obs.groupby(["Condition", f"leiden_{res}"]).size().unstack()
    df = df.div(df.sum(axis=1), axis=0)
    fig, ax = plt.subplots(figsize=(12, 2))
    sns.heatmap(df, cmap="viridis", annot=True, fmt=".2f", ax=ax)
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Condition")
    fig.tight_layout()
    fig.suptitle(f"Cluster fractions per condition (resolution {res})")
    fig.savefig(plot_dir / f"cluster_fractions_per_condition_heatmap_res_{res}.pdf")

# compute average features per cluster
for obsm_key in obsm_to_use:
    df = adata.obsm[obsm_key].copy()
    meta_keys = ["sample_id", "Condition", "lesion_type"] + ["leiden_" + str(res) for res in resolutions]
    df = df.join(adata.obs[meta_keys], on="cell")
    df.reset_index(inplace=True)
    df = df.melt(id_vars=["cell"]+meta_keys, var_name="feature", value_name="value")
    for res in resolutions:
        df_tmp = df.groupby(["feature", f"leiden_{res}"]).mean().reset_index()
        df_tmp["min_max_value"] = df_tmp.groupby(["feature"]).value.apply(lambda x: (x - x.min()) / (x.max() - x.min()))
        fig, ax = plt.subplots(1, 2, figsize=(24, 12))
        sns.heatmap(data=df_tmp.pivot(index="feature", columns=f"leiden_{res}", values="value"), 
                    cmap="viridis", ax=ax[0])
        sns.heatmap(data=df_tmp.pivot(index="feature", columns=f"leiden_{res}", values="min_max_value"), 
                cmap="viridis", ax=ax[1])
        ax[0].set_title("Raw values")
        ax[1].set_title("Min-max scaled values per feature")
        fig.savefig(plot_dir / f"feature_heatmap_{obsm_key}_res_{res}.pdf")

# compute feature loadings per factor
for obsm_key in obsm_to_use:
    features = obsm_features[obsm_key]
    fig, ax = plt.subplots(figsize=(6, 10))
    sns.heatmap(model.get_weights(df=True).loc[features, :], ax=ax, cmap="coolwarm", center=0)
    fig.tight_layout()
    fig.savefig(plot_dir / f"feature_loadings_{obsm_key}.pdf")

# check the fraction of variance explained per group per factor per view
fig = mfx.plot.plot_r2(model, y="Group", x="Factor")
fig.tight_layout()
fig.savefig(plot_dir / "r2_per_group_per_factor_per_view.pdf")

