
# usage
# python scripts/process/process_vis.py

import pandas as pd
import numpy as np
import os
import scanpy as sc
from pathlib import Path

from umap import UMAP
import matplotlib.pyplot as plt
import seaborn as sns
import decoupler as dc

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
out_file = current_path / ".." / ".." / "data" / "prc" / "vis" / "processed"
out_file.mkdir(parents=True, exist_ok=True)
visium_path = current_path / ".." / ".." / "data" / "raw" / "vis"
c2l_path = current_path / ".." / ".." / "data" / "prc" / "vis" / "c2l_out" / "cellranger"
img_features = current_path / ".." / ".." / "data" / "prc" / "images" / "squdipy_features"
visium_samples = [f for f in os.listdir(visium_path) if not f.startswith(".")]
print(np.array(visium_samples))

sample_meta = pd.read_excel(current_path / ".." / ".." / "data" / "Metadata_all.xlsx", sheet_name="Visium")
sample_meta


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
    for model_call in ["all", "condition", "lesion_type"]:
        suffix = lookup(model_call, sample_entry)
        for output in ["abunds", "props"]:
            m = pd.read_csv(c2l_path / sample_id / f"cell_{output}_{suffix}.csv", index_col=0)
            inter = slide.obs.index.intersection(m.index)
            slide.obsm[f"{output}_{model_call}"] = m.loc[inter]

    # Read image features
    adata_img = sc.read_h5ad(img_features / f"{sample_id}.h5ad")
    for feature in ["summary", "histogram", "texture"]:
        m = adata_img.obsm[feature]
        inter = slide.obs.index.intersection(m.index)
        slide.obsm[feature] = m.loc[inter]

    return slide

# load the visium data
vis_dict = {s: read_slide(s, visium_path, c2l_path) for s in visium_samples}

# get pathways
msigdb = dc.get_resource('MSigDB')

# get hallmark db
hallmark = msigdb[msigdb['collection']=='hallmark'] # filter by hallmark
hallmark = hallmark[~hallmark.duplicated(['geneset', 'genesymbol'])] # remove duplicates
hallmark.loc[:, 'geneset'] = [name.split('HALLMARK_')[1] for name in hallmark['geneset']] # rename for consistency
hallmark = hallmark.loc[:, ['geneset', 'genesymbol']] # reorder columns

# get progeny db
progeny = dc.get_progeny(top=300)
progeny = progeny.rename(columns={'source': 'geneset', 'target': 'genesymbol'})
progeny = progeny.loc[:, ['geneset', 'genesymbol', 'weight']] # reorder columns

# get reactome db
reactome = msigdb[msigdb['collection'] == 'reactome_pathways']
reactome = reactome.loc[:, ['geneset', 'genesymbol']] # reorder columns
reactome = reactome[~reactome.duplicated(['geneset', 'genesymbol'])]

# infer pathway activity per slide
for key, adata in vis_dict.items():
    print(key)
    for pkn, pkn_name in zip([hallmark, progeny, reactome], ["hallmark", "progeny", "reactome"]):
        print(pkn_name)
        dc.run_ulm(
            mat=adata,
            net=pkn,
            source="geneset",
            target="genesymbol",
            weight="weight" if pkn_name in ["progeny"] else None,
            verbose=True,
            use_raw=True)
        adata.obsm[f"{pkn_name}_estimates"] = adata.obsm["ulm_estimate"]
        adata.obsm[f"{pkn_name}_pvals"] = adata.obsm["ulm_pvals"]
        del adata.obsm["ulm_estimate"], adata.obsm["ulm_pvals"]

# save the adata objects
for key, adata in vis_dict.items():
    adata.write(out_file / f"{key}.h5ad")
