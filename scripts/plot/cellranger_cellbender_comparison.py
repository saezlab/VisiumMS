
import scanpy as sc
import scanpy.external as sce

import numpy as np
import pandas as pd

import os
from pathlib import Path

import seaborn as sns
import matplotlib.pyplot as plt

import plotnine as pn
import re

current_folder = Path(__file__).parent
raw_input = current_folder / ".." / ".." / "data" / "raw" / "sc"
cellbender_input = current_folder / ".." / ".." / "data" / "prc" / "sc" / "cellbender_qc"
cellranger_input = current_folder / ".." / ".." / "data" / "prc" / "sc" / "cellranger_qc"
plot_out = current_folder / ".." / ".." / "out" / "cellbender_cellranger_comparison"
plot_out.mkdir(parents=True, exist_ok=True)
sc._settings.ScanpyConfig.figdir = plot_out

samples = [s.split(".")[0] for s in os.listdir(cellbender_input) if s.endswith(".h5ad")]

cellbender_adata_list = [sc.read_h5ad(cellbender_input / (sample + ".h5ad")) for sample in samples]
adata_cellbender = sc.concat(cellbender_adata_list, join="outer", label="sample_id", keys=samples)
del cellbender_adata_list
adata_cellbender.obs.index = [smp + "_" + bc for smp, bc in zip(adata_cellbender.obs["sample_id"], adata_cellbender.obs.index)]

cellranger_adata_list = [sc.read_h5ad(cellranger_input / (sample + ".h5ad")) for sample in samples]
adata_cellranger = sc.concat(cellranger_adata_list, join="outer", label="sample_id", keys=samples)
del cellranger_adata_list
adata_cellranger.obs.index = [smp + "_" + bc for smp, bc in zip(adata_cellranger.obs["sample_id"], adata_cellranger.obs.index)]

# add comparison column to obs
adata_cellranger.obs["in_cellbender_out"] = np.isin(adata_cellranger.obs.index.to_numpy(), adata_cellbender.obs.index.to_numpy())
adata_cellbender.obs["in_cellranger_out"] = np.isin(adata_cellbender.obs.index.to_numpy(), adata_cellranger.obs.index.to_numpy())

# reading in the unfiltered data
adata_raw = []
for sample in samples:
    print(sample)
    barcodes_oi = set([re.sub("^" + sample + "_", "", s) for s in adata_cellbender.obs_names[adata_cellbender.obs["sample_id"] == sample]]) | \
        set([re.sub("^" + sample + "_", "", s) for s in adata_cellranger.obs_names[adata_cellranger.obs["sample_id"] == sample]])
    barcodes_oi = np.array(list(barcodes_oi))
    adata = sc.read_10x_h5(raw_input / sample / "raw_feature_bc_matrix.h5")
    adata = adata[barcodes_oi, :]
    adata.var_names_make_unique()
    adata.obs.index = [sample + "_" + bc for bc in adata.obs.index]
    adata_raw.append(adata)
adata_raw = sc.concat(adata_raw, join="outer", label="sample_id", keys=samples)

# compute QC metrics
adata_raw.var['mt'] = adata_raw.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata_raw, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Compute doublets score
sce.pp.scrublet(adata_raw, verbose=False)

# Compute dissociation score and normalize it (0-1)
diss_path = "https://raw.githubusercontent.com/kieranrcampbell/scrnaseq-digestion-paper/master/data/deliverables/coregene_df-FALSE-v3.csv"
diss_genes = pd.read_csv(diss_path).sort_values('PValue').head(200).gene_symbol.tolist()
intersect_genes = adata_raw.var.index.intersection(diss_genes)
if len(intersect_genes) > 1:
    sc.tl.score_genes(adata_raw, gene_list=intersect_genes, ctrl_size=len(intersect_genes), 
                        score_name='diss_score')
    adata_raw.obs['diss_score'] = (adata_raw.obs.diss_score-np.min(adata_raw.obs.diss_score)) / \
    (np.max(adata_raw.obs.diss_score) - np.min(adata_raw.obs.diss_score))
else:
    adata_raw.obs['diss_score'] = 0

# check whether a barcode/cell is present in cellranger, cellbender or both
adata_raw.obs["in_cellranger_out"] = np.isin(adata_raw.obs.index.to_numpy(), adata_cellranger.obs.index.to_numpy())
adata_raw.obs["in_cellbender_out"] = np.isin(adata_raw.obs.index.to_numpy(), adata_cellbender.obs.index.to_numpy())
adata_raw.obs["in_both"] = np.logical_and(adata_raw.obs["in_cellranger_out"], adata_raw.obs["in_cellbender_out"])

cell_labels = []
# loop through the rows of the dataframe
for i, row in adata_raw.obs.iterrows():
    if row["in_both"]:
        cell_labels.append("both")
    elif row["in_cellranger_out"]:
        cell_labels.append("only_cellranger")
    elif row["in_cellbender_out"]:
        cell_labels.append("only_cellbender")
    else:
        raise ValueError("This should not happen")
adata_raw.obs["label"] = cell_labels

# verbosity
print(adata_raw.obs.value_counts("label"))

# global comparison
fig, axis = plt.subplots(2, 3, figsize=(20, 10))
sns.set_theme(style="whitegrid")
sns.countplot(data=adata_raw.obs, x="label", ax=axis[0, 0])
sns.violinplot(data=adata_raw.obs, x="label", y="n_genes_by_counts", ax=axis[0, 1])
sns.violinplot(data=adata_raw.obs, x="label", y="total_counts", ax=axis[0, 2])
sns.violinplot(data=adata_raw.obs, x="label", y="pct_counts_mt", ax=axis[1, 0])
sns.violinplot(data=adata_raw.obs, x="label", y="doublet_score", ax=axis[1, 1])
sns.violinplot(data=adata_raw.obs, x="label", y="diss_score", ax=axis[1, 2])
plt.savefig(plot_out / "global_comparison.png", dpi=300)
plt.close()

# comparison per sample post QC
df = adata_cellranger.obs.sample_id.value_counts().to_frame().reset_index().rename(columns={"count": "n_cellranger"}).merge(
    adata_cellbender.obs.sample_id.value_counts().to_frame().reset_index().rename(columns={"count": "n_cellbender"}),
    on="sample_id", how="outer"
)

plt.figure(figsize=(10, 10))
sns.set_theme(style="whitegrid")
plt.plot([0, 15000], [0, 15000], linewidth=2, color="black")
sns.scatterplot(data=df, x="n_cellranger", y="n_cellbender")
for i, row in df.iterrows():
    plt.text(row["n_cellranger"], row["n_cellbender"], row["sample_id"])
plt.xlabel("Number of cells in cellranger analysis")
plt.ylabel("Number of cells in CellBender analysis")
plt.xlim(1, 15000)
plt.ylim(1, 15000)
plt.savefig(plot_out / "sample_scatter.png", dpi=300)
plt.close()

df_list = []
for smp in np.unique(adata_cellbender.obs.sample_id):
    cellbender_set = set(adata_cellbender.obs_names[adata_cellbender.obs.sample_id == smp])
    prev_set = set(adata_cellranger.obs_names[adata_cellranger.obs.sample_id == smp])
    df_list.append(pd.DataFrame({"sample_id": smp, "both": len(cellbender_set.intersection(prev_set)), "cellbender_only": len(cellbender_set.difference(prev_set)), "cellranger_only": len(prev_set.difference(cellbender_set))}, index=[0]))
df = pd.concat(df_list)
df = df.melt(id_vars="sample_id", var_name="comparison", value_name="count")

plot = pn.ggplot(df,
          pn.aes(x= "comparison", y="count", fill="comparison")) + \
          pn.geom_bar(stat="identity") + \
          pn.facet_wrap("sample_id") + \
          pn.theme(axis_text_x=pn.element_text(rotation=60, hjust=1)) + \
          pn.theme(figure_size=(10, 10)) 
plot.save(plot_out / "sample_boxplot.png", dpi=300)

# total UMI counts per sample per method
sns.displot(data=adata_raw.obs, x="total_counts", hue="label", kind="kde", common_norm = False, fill=True, col="sample_id", col_wrap=3, facet_kws={"sharey": False}, log_scale=True)
plot.save(plot_out / "total_UMI_per_sample_per_method.png", dpi=300)
plt.close()