
# to run the script from the project directory:
# python scripts/process/deep_image_features_eval.py --quality hires
# python scripts/process/deep_image_features_eval.py --quality fullres

from pathlib import Path
import numpy as np
import pandas as pd
import scanpy
import matplotlib.pyplot as plt
import seaborn as sns
import os
import scanpy as sc
from sklearn.decomposition import PCA
import re
import argparse

# command line arguments for quality to use
parser = argparse.ArgumentParser()
parser.add_argument("--quality", type=str, required=True)
args = parser.parse_args()

if args.quality not in ["lowres", "hires", "fullres"]:
    raise ValueError("quality must be one of 'lowres', 'hires' or 'fullres'")
QUALITY = args.quality

# set up relative paths within the project
current_folder = Path(__file__).parent
input_dir = current_folder / ".." / ".." / "data" / "prc" / "images" / ("image_features_" + QUALITY)
output_dir = current_folder / ".." / ".." / "out" / "deep_image_features_eval"
output_dir.mkdir(parents=True, exist_ok=True)

# load metadata
metadata_path = current_folder / ".." / ".." / "data" / "Metadata_all.xlsx"
metadata = pd.read_excel(metadata_path, index_col=0, sheet_name="Visium")
sample_lesion_dict = dict(zip(metadata.sample_id, metadata["lesion_type"]))

# load the samples
samples = [re.sub(".h5ad", "", f) for f in os.listdir(input_dir) if not f.startswith(".") and f.endswith(".h5ad")]
adata_dict = {sample: sc.read_h5ad(input_dir / (sample + ".h5ad")) for sample in samples}

# get condition vectors for plotting
cond_vec = np.array(["MS" if s.startswith("MS") else "HC" for s in adata_dict.keys()])
stage_vec = np.array([sample_lesion_dict[sample] if sample in sample_lesion_dict.keys() else "NA" for sample in adata_dict.keys()])

### 1. mean image features per sample ###
mean_image_features = np.array([adata.obsm["image_features"].mean(axis=0) for adata in adata_dict.values()])
pca_coord = PCA(n_components=min(mean_image_features.shape)).fit_transform(mean_image_features)

fig, ax = plt.subplots(figsize=(10, 10))
sns.scatterplot(x=pca_coord[:, 0], y=pca_coord[:, 1], hue=stage_vec)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
fig.savefig(output_dir / ("pca_mean_image_features_" + QUALITY + ".png"), dpi=300)
plt.close()
### ------------------------------- ###

### 2. cluster cells based on image feature space and compute fractions of clusters per image ###
image_features_mtx = np.vstack([adata.obsm["image_features"] for adata in adata_dict.values()])
image_features_adata = scanpy.AnnData(image_features_mtx)
image_features_adata.obs = pd.concat([adata.obs for adata in adata_dict.values()])
image_features_adata.obs["sample_id"] = np.hstack([np.repeat(sample, len(adata_dict[sample])) for sample in samples])
image_features_adata.obs["lesion_type"] = image_features_adata.obs["sample_id"].map(sample_lesion_dict)

sc.pp.pca(image_features_adata)
sc.external.pp.bbknn(image_features_adata, batch_key="sample_id")
sc.tl.leiden(image_features_adata, resolution=2)

cluster_counts = image_features_adata.obs.groupby(["sample_id", "leiden"]).size().reset_index()
cluster_counts = cluster_counts.pivot(index="sample_id", columns="leiden", values=0)
cluster_counts = cluster_counts.fillna(0)
cluster_counts = cluster_counts.div(cluster_counts.sum(axis=1), axis=0)
pca_coord = PCA(n_components=min(cluster_counts.shape)).fit_transform(cluster_counts)

fig, ax = plt.subplots(figsize=(10, 10))
sns.clustermap(cluster_counts, cmap="Blues")
fig.savefig(output_dir / ("ht_clustered_image_features_" + QUALITY + ".png"), dpi=300)
plt.close()

fig, ax = plt.subplots(figsize=(10, 10))
sns.scatterplot(x=pca_coord[:, 0], y=pca_coord[:, 1], hue=stage_vec)
for i, txt in enumerate(cluster_counts.index):
    plt.annotate(txt, (pca_coord[i, 0], pca_coord[i, 1]))
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
fig.savefig(output_dir / ("pca_clustered_image_features_" + QUALITY + ".png"), dpi=300)
plt.close()
### ------------------------------- ###