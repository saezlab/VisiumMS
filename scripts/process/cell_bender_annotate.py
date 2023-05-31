
import scanpy as sc
import anndata as ad

from sklearn.neighbors import NearestNeighbors
import numpy as np
import pandas as pd

import os
from pathlib import Path
import re

# TODO: This hypterparameter might be too high for lowly abundant cell types
NN = 10

current_folder = Path(__file__).parent
cellbender_out = current_folder / ".." / ".." / "data" / "cellbender_out"
raw_out = current_folder / ".." / ".." / "data" / "raw" / "sc"

samples = [sample for sample in os.listdir(cellbender_out) if not sample.startswith(".")]

# load the annoated single-nuc data
annotated = sc.read_h5ad(current_folder / ".." / ".." / "data" / "uscsc_dump" / "annotated.h5ad")
annotated.obs.index = [re.sub("-[0-9]+$", "", barcode) for barcode in annotated.obs.index]
hvg = annotated.var_names.to_numpy()

# load the cellbender single-nuc output
adata_no_annot_list = {sample: sc.read_h5ad(cellbender_out / sample / "cell_bender_matrix_filtered_qc.h5") for sample in samples}

# load raw single-nuc data
raw_adata_objects = {sample: sc.read_10x_mtx(raw_out / sample / "filtered_feature_bc_matrix") for sample in samples}

# TODO: for now I write lousy logging files
for sample in samples:

    with open(cellbender_out / sample / "cell_bender_annotation.log", "w") as f:

        sc._settings.ScanpyConfig.figdir = cellbender_out / sample  # where scanpy saves plots

        print(sample)

        adata_no_annot = adata_no_annot_list[sample]
        #sc.pp.normalize_total(adata_no_annot, target_sum=1e4) # done below
        #sc.pp.log1p(adata_no_annot)                           # done below

        #adata_annot = annoated[annoated.obs.sample_id==sample, ]  
        # NOTE: We use the raw data here!
        annotated_ids = annotated.obs.index[annotated.obs.sample_id==sample]
        adata_annot = raw_adata_objects[sample]
        annotated_ids = np.intersect1d(annotated_ids, adata_annot.obs.index)
        adata_annot = adata_annot[annotated_ids, hvg]
        obs = annotated.obs[annotated.obs.sample_id==sample]
        adata_annot.obs = obs.loc[annotated_ids, :]

        # get the set of shared var_names
        shared_var_names = list(set(adata_annot.var_names) & set(adata_no_annot.var_names))
        adata_annot = adata_annot[:, shared_var_names]
        adata_no_annot_tmp = adata_no_annot[:, shared_var_names].copy()  # copy, because I want to keep all genes for later on

        # how many cells are in adata_no_annot that are not in adata_annot?
        print("Number of cells in adata_no_annot: " + str(len(adata_no_annot_tmp.obs_names)))
        print("Number of cells in adata_no_annot: " + str(len(adata_no_annot_tmp.obs_names)), file=f)
        print("Number of cells in adata_annot: " + str(len(adata_annot.obs_names)))
        print("Number of cells in adata_annot: " + str(len(adata_annot.obs_names)), file=f)

        cells_oi = set(adata_no_annot_tmp.obs_names) - set(adata_annot.obs_names)
        print("Number of cells in adata_no_annot that are not in adata_annot: " + str(len(cells_oi)))
        print("Number of cells in adata_no_annot that are not in adata_annot: " + str(len(cells_oi)), file=f)

        # how many cells are in adata_no_annot that are not in adata_annot?
        cells_oi = set(adata_annot.obs_names) - set(adata_no_annot_tmp.obs_names)
        print("Number of cells in adata_annot that are not in adata_no_annot: " + str(len(cells_oi)))
        print("Number of cells in adata_annot that are not in adata_no_annot: " + str(len(cells_oi)), file=f)

        # add prefixes to cell barcodes and concat
        adata_annot.obs_names = ["annot_" + name for name in adata_annot.obs_names]
        adata_no_annot_tmp.obs_names = ["not_annot_" + name for name in adata_no_annot_tmp.obs_names]
        adata_concat = ad.concat([adata_annot, adata_no_annot_tmp], join="outer", label="annotation", keys=["annotated", "unannotated"])

        # NOTE: integration should not be neccessary
        #sc.external.pp.bbknn(adata_concat, batch_key="annotation")

        # normalize and joint pca
        sc.pp.normalize_total(adata_concat, target_sum=1e4)
        sc.pp.log1p(adata_concat)
        sc.pp.pca(adata_concat)

        # umap for QC, one should see high overlap between annotated and non-annotated because they are mostly the same cells
        sc.pp.neighbors(adata_concat)
        sc.tl.umap(adata_concat)
        sc.pl.umap(adata_concat, color="annotation", save="_joint_umap.png")

        # split back into annotated and unannotated and remove prefixes
        adata_annot_tmp = adata_concat[adata_concat.obs.index.str.startswith("annot_"), ]
        adata_annot_tmp.obs_names = [re.sub("^annot_", "", name) for name in adata_annot_tmp.obs_names]
        adata_no_annot_tmp = adata_concat[adata_concat.obs.index.str.startswith("not_annot_"), ]
        adata_no_annot_tmp.obs_names = [re.sub("^not_annot_", "", name) for name in adata_no_annot_tmp.obs_names]

        # for each row in adata_not_annot_pca search for the 20 nearest neighbors in adata_annot_pca
        nbrs = NearestNeighbors(n_neighbors=NN, algorithm='ball_tree').fit(adata_annot_tmp.obsm["X_pca"])
        distances, indices = nbrs.kneighbors(adata_no_annot_tmp.obsm["X_pca"])
        
        # NOTE: we will add a cautionary annotation and a forced annotation, where we don't care for the NN/2 threshold
        annotation, annotation_forced = [], []
        for cell_i in range(len(adata_no_annot)):
            barcode = adata_no_annot.obs.index[cell_i]
            # if the barcode is present in the annotated object, then use the annotation
            if barcode in adata_annot.obs.index:
                annotation.append(adata_annot.obs.loc[barcode, "cell_type"])
                annotation_forced.append(adata_annot.obs.loc[barcode, "cell_type"])
            # else, use the annotation of the nearest neighbors
            else:
                nearest_neighbors = indices[cell_i, :]
                nn_annot = adata_annot.obs.iloc[nearest_neighbors, :]["cell_type"].value_counts()
                # if max count below 50%, then set to "unannotated"
                if nn_annot.iloc[0] < NN/2:
                    annotation.append("unannotated")
                else:
                    # get label with max count
                    annotation.append(nn_annot.index[0])
                annotation_forced.append(nn_annot.index[0])
        adata_no_annot.obs["cell_type"] = annotation

        print(adata_no_annot.obs["cell_type"].value_counts())
        print(adata_no_annot.obs["cell_type"].value_counts(), file=f)

        print(adata_no_annot)

        # save the annotated object
        adata_no_annot.write_h5ad(cellbender_out / sample / "cell_bender_matrix_filtered_qc_annotated.h5ad")

        # save the forced annotation object
        adata_no_annot.obs["cell_type"] = annotation_forced
        adata_no_annot.write_h5ad(cellbender_out / sample / "cell_bender_matrix_filtered_qc_annotated_forced.h5ad")