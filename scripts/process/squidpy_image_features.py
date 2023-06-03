
# usage
# python scripts/process/squidpy_image_features.py

from pathlib import Path
import numpy as np
import pandas as pd
import scanpy
import matplotlib.pyplot as plt
import seaborn as sns
import os
import scanpy as sc
from sklearn.decomposition import PCA
from umap import UMAP
import squidpy as sq

plt.rcParams['figure.figsize'] = [10, 10]

current_folder = Path(__file__).parent
image_dir = current_folder / ".." / ".." / "data" / "raw" / "images"
visium_dir = current_folder / ".." / ".." / "data" / "raw" / "vis"

image_features_out = current_folder / ".." / ".." / "data" / "prc" / "images" / "squdipy_features"
image_features_out.mkdir(parents=True, exist_ok=True)

samples = [f for f in os.listdir(image_dir) if not f.startswith(".")]

for sample in samples:

    print(sample)

    adata = sc.read_visium(visium_dir / sample / "outs")
    adata.var_names_make_unique()
    adata.obsm["spatial"] = adata.obsm["spatial"].astype(int) # convert str to int
    print(adata)

    # load image
    img = sq.im.ImageContainer(image_dir / sample / (sample + "_pic.tif"))
    print(img.shape)

    # compute features
    for feature in ["texture", "summary", "histogram"]:
        print(feature)
        sq.im.calculate_image_features(adata, img, features=feature, key_added=feature, n_jobs=8, show_progress_bar=True)
    
    # compute QC metrics for each spot
    adata.obs["detected_genes"] = np.sum(adata.X > 0, axis=1)
    adata.obs["total_umis"] = np.sum(adata.X, axis=1)

    # save object
    adata.write(image_features_out / (sample + ".h5ad"))

    # make QC plots
    for feature in ["texture", "summary", "histogram"]:
        X = adata.obsm[feature].to_numpy()
        pca_coord = PCA(n_components=10).fit_transform(X)
        umap_coord = UMAP(n_components=2).fit_transform(pca_coord)

        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        sns.scatterplot(
            x=umap_coord[:, 0],
            y=umap_coord[:, 1],
            hue=adata.obs["detected_genes"],
            s=20)
        plt.suptitle(sample)
        fig.savefig(image_features_out / (sample + f"_{feature}_umap_detected_genes.png"), dpi=300)
        plt.close()

        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        sns.scatterplot(
            x=umap_coord[:, 0],
            y=umap_coord[:, 1],
            hue=adata.obs["total_umis"],
            s=20)
        plt.suptitle(sample)
        fig.savefig(image_features_out / (sample + f"_{feature}_umap_total_umis.png"), dpi=300)
        plt.close()

