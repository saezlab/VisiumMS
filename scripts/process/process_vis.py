
# usage
# python scripts/process/process_vis.py --output cellbender
# python scripts/process/process_vis.py --output cellranger

import pandas as pd
import numpy as np
import os
import scanpy as sc
from pathlib import Path
import re

from umap import UMAP
import matplotlib.pyplot as plt
import seaborn as sns
import decoupler as dc
from composition_stats import closure, ilr
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--output", type=str, required=True)
args = parser.parse_args()

# see deconv.py
def lookup(model_call, sample_entry):
    if model_call == "all":
        reg_model = "All"
    elif model_call == "condition":
        if sample_entry["Condition"] == "Control":
            reg_model = "Control"
        elif sample_entry["Condition"] == "MS":
            reg_model = "MS"
        else:
            raise ValueError("Unknown condition")
    elif model_call == "lesion_type":
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
    return reg_model

current_path = Path(__file__).parent
visium_path = current_path / ".." / ".." / "data" / "raw" / "vis"
img_features = current_path / ".." / ".." / "data" / "prc" / "images" / "squdipy_features"
visium_samples = [f for f in os.listdir(visium_path) if not f.startswith(".")]
print(np.array(visium_samples))

if args.output == "cellbender":
    c2l_path = current_path / ".." / ".." / "data" / "prc" / "vis" / "c2l_out" / "cellbender"
    out_file = current_path / ".." / ".." / "data" / "prc" / "vis" / "processed" / "cellbender"
    out_file.mkdir(parents=True, exist_ok=True)
elif args.output == "cellranger":
    c2l_path = current_path / ".." / ".." / "data" / "prc" / "vis" / "c2l_out" / "cellranger"
    out_file = current_path / ".." / ".." / "data" / "prc" / "vis" / "processed" / "cellranger"
    out_file.mkdir(parents=True, exist_ok=True)
else:
    raise ValueError("Unknown output")

sample_meta = pd.read_excel(current_path / ".." / ".." / "data" / "Metadata_all.xlsx", sheet_name="Visium")

def read_slide(sample_id, visium_path, c2l_path):

    # get sample metadata
    sample_entry = sample_meta.loc[sample_meta.sample_id == sample_id, :].to_dict(orient="records")[0]

    # Read rna-seq
    slide = sc.read_visium(visium_path / sample_id / "outs")
    slide.var_names_make_unique()
    
    sc.pp.filter_genes(slide, min_cells=3)
    sc.pp.filter_cells(slide, min_genes=200)

    # Store raw counts
    slide.raw = slide
    slide.layers["counts"] = slide.X.copy()

    # Normalize
    sc.pp.normalize_total(slide, target_sum=1e4)
    sc.pp.log1p(slide)

    # Read props and abunds
    for model_call in ["all", "condition"]:
        suffix = lookup(model_call, sample_entry)
        m = pd.read_csv(c2l_path / sample_id / f"cell_abunds_{suffix}_q05_cell_abundance_w_sf.csv", index_col=0)
        inter = slide.obs.index.intersection(m.index)
        slide.obsm[f"abunds_{model_call}"] = m.loc[inter]
        # check this gives us proportions
        slide.obsm[f"props_{model_call}"] = slide.obsm[f"abunds_{model_call}"].div(slide.obsm[f"abunds_{model_call}"].sum(axis=1), axis=0)
        slide.obsm[f"props_{model_call}"].loc[:, :] = closure(slide.obsm[f"props_{model_call}"].values)
        slide.obsm[f"props_ilr_{model_call}"] = ilr(slide.obsm[f"props_{model_call}"].values)

    # Read image features
    adata_img = sc.read_h5ad(img_features / f"{sample_id}.h5ad")
    for feature in ["summary", "histogram", "texture"]:
        m = adata_img.obsm[feature]
        inter = slide.obs.index.intersection(m.index)
        slide.obsm[feature] = m.loc[inter]

    return slide

# load the visium data
vis_dict = {s: read_slide(s, visium_path, c2l_path) for s in visium_samples}

# also common leiden clustering to evaluate the image results
vis_all = sc.AnnData.concatenate(*vis_dict.values(), batch_key="sample_id", batch_categories=visium_samples)
sc.pp.highly_variable_genes(vis_all, n_top_genes=2000)
sc.pp.scale(vis_all, max_value=10)
sc.pp.pca(vis_all, n_comps=40)
sc.external.pp.bbknn(vis_all, batch_key='sample_id') 
sc.pp.neighbors(vis_all, n_neighbors=10, n_pcs=40)
sc.tl.leiden(vis_all, resolution=0.5, key_added="leiden_0.50")
sc.tl.leiden(vis_all, resolution=0.25, key_added="leiden_0.25")
sc.tl.leiden(vis_all, resolution=0.1, key_added="leiden_0.10")

for sample_id in vis_dict.keys():
    df = vis_all.obs[vis_all.obs.sample_id == sample_id].copy()
    df.index = [re.sub(f"-{sample_id}$", "", i) for i in df.index]
    for res in ["leiden_0.50", "leiden_0.25", "leiden_0.10"]:
        vis_dict[sample_id].obs[res] = df.loc[vis_dict[sample_id].obs.index][res]
del vis_all

# get pathways
msigdb = dc.get_resource('MSigDB')

# get hallmark db
hallmark = msigdb[msigdb['collection']=='hallmark'] # filter by hallmark
hallmark = hallmark[~hallmark.duplicated(['geneset', 'genesymbol'])] # remove duplicates
hallmark.loc[:, 'geneset'] = [name.split('HALLMARK_')[1] for name in hallmark['geneset']] # rename for consistency
hallmark = hallmark.loc[:, ['geneset', 'genesymbol']] # reorder columns

# get progeny db
progeny = dc.get_progeny(top=300) # only possible for progeny because we weights are not available for hallmark or reactome
progeny = progeny.rename(columns={'source': 'geneset', 'target': 'genesymbol'})
progeny = progeny.loc[:, ['geneset', 'genesymbol', 'weight']] # reorder columns

# get reactome db
reactome = msigdb[msigdb['collection'] == 'reactome_pathways']
reactome = reactome.loc[:, ['geneset', 'genesymbol']] # reorder columns
reactome = reactome[~reactome.duplicated(['geneset', 'genesymbol'])]

# infer pathway activity per spot using ulm and wmean (for now remove reactome)
for key, adata in vis_dict.items():
    print(key)

    #for pkn, pkn_name in zip([hallmark, progeny, reactome], ["hallmark", "progeny", "reactome"]):  
    for pkn, pkn_name in zip([hallmark, progeny], ["hallmark", "progeny"]):
        print(pkn_name)

        # run ulm
        dc.run_ulm(
            mat=adata,
            net=pkn,
            source="geneset",
            target="genesymbol",
            weight="weight" if pkn_name in ["progeny"] else None,
            verbose=True,
            use_raw=True)
        adata.obsm[f"{pkn_name}ulm_estimate"] = adata.obsm["ulm_estimate"]
        adata.obsm[f"{pkn_name}_ulm_pvals"] = adata.obsm["ulm_pvals"]
        del adata.obsm["ulm_estimate"], adata.obsm["ulm_pvals"]
        
        # run wmean
        dc.run_wmean(
            mat=adata,
            net=pkn,
            source="geneset",
            target="genesymbol",
            weight="weight" if pkn_name in ["progeny"] else None,
            verbose=True,
            use_raw=True)
        adata.obsm[f"{pkn_name}_wmean_estimate"] = adata.obsm["wmean_estimate"]
        adata.obsm[f"{pkn_name}_wmean_norm"] = adata.obsm["wmean_norm"]
        adata.obsm[f"{pkn_name}_wmean_corr"] = adata.obsm["wmean_corr"]
        adata.obsm[f"{pkn_name}_wmean_pvals"] = adata.obsm["wmean_pvals"]
        del adata.obsm["wmean_estimate"], adata.obsm["wmean_norm"], adata.obsm["wmean_corr"], adata.obsm["wmean_pvals"]

# save the adata objects
for key, adata in vis_dict.items():
    adata.write(out_file / f"{key}.h5ad")
