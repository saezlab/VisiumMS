<<<<<<< HEAD

# usage
# python scripts/process/eval_deconv.py

# TODO:
# repeat evaluation using the posterior means instead of the 5% quantile

=======
>>>>>>> 4059a4c8e2af9c39af787ebee1439fc854d311d6
import pandas as pd
import numpy as np
import scanpy as sc
import os
from anndata import AnnData
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp
<<<<<<< HEAD
from pathlib import Path

# function to annotate seaborn scatterplots with pearson correlation coefficient
=======
import argparse


# Read command line and set args
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--path_slides', help='Input path to raw slides', required=True)
parser.add_argument('-a', '--path_ann', help='Annotated AnnData object', required=True)
parser.add_argument('-d', '--path_deconv', help='Path to deconv results', required=True)
parser.add_argument('-o', '--output_path', help='Path were to save plots', required=True)
args = vars(parser.parse_args())

path_slides = args['path_slides']
path_ann = args['path_ann']
path_deconv = args['path_deconv']
output_path = args['output_path']
###############################

# Read annotated sc atlas
meta = sc.read(path_ann).obs

# Compute proportions for sc
sc_df = (meta
 .groupby(['sample_id', 'leiden'])[['sample_id']]
 .count()
 .rename({'sample_id': 'n'}, axis=1)
 .reset_index()
 .assign(sc=lambda x: x['n'] / x.groupby('sample_id')['n'].transform('sum'))
 .drop('n', axis=1)
)

# Read deconv results
prop = []
sample_ids = []
for sample_id in os.listdir(path_deconv):
    if os.path.isdir(os.path.join(path_deconv, sample_id)) and not sample_id.startswith('.') and sample_id != 'reg_model':
        sample_ids.append(sample_id)
        tmp = AnnData(pd.read_csv(os.path.join(path_deconv, '{0}/cell_props.csv'.format(sample_id)), index_col=0))
        tmp.obs['sample_id'] = sample_id
        prop.append(tmp)
prop = prop[0].concatenate(prop[1:])

# Compute average props for visium
vm_df = pd.DataFrame(prop.X, index=prop.obs.index, columns=prop.var.index)
vm_df['sample_id'] = prop.obs['sample_id']
vm_df = vm_df.groupby('sample_id').mean(1)
vm_df = vm_df.melt(value_vars=vm_df.columns, ignore_index=False).reset_index()
vm_df = vm_df.rename({'variable': 'leiden', 'value': 'vm'}, axis=1)

# Merge dataframes
df = pd.merge(sc_df, vm_df, how='inner')
df['sc'] = np.log10(df['sc'].replace({0: np.nan}))
df['vm'] = np.log10(df['vm'].replace({0: np.nan}))
df = df.fillna(-4)

# Plot correlations:
>>>>>>> 4059a4c8e2af9c39af787ebee1439fc854d311d6
def annotate(data, **kws):
    x = data['sc'].values
    y = data['vm'].values
    msk = np.isfinite(x) & np.isfinite(y)
    r, p = sp.stats.pearsonr(x[msk], y[msk])
    ax = plt.gca()
    ax.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(r, p), transform=ax.transAxes)
<<<<<<< HEAD


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


def get_sc_proportions(meta):
    return (meta
     .groupby(['sample_id', 'cell_types'])[['sample_id']]
     .count()
     .rename({'sample_id': 'n'}, axis=1)
     .reset_index()
     .assign(sc=lambda x: x['n'] / x.groupby('sample_id')['n'].transform('sum'))
     .drop('n', axis=1)
    )


def get_vis_proportions(model_call, method, sample_meta, measure):
    assert model_call in ["all", "condition", "lesion_type"]
    assert method in ["cellbender", "cellranger"]
    assert measure in ["abunds", "props"]
    prop = []
    for sample_id in os.listdir(path_deconv / method):
        sample_entry = sample_meta.loc[sample_meta.sample_id == sample_id, :].to_dict(orient="records")[0]
        suffix = lookup(model_call, sample_entry)
        tmp = AnnData(X=pd.read_csv(path_deconv / method / sample_id / f"cell_{measure}_{suffix}.csv", index_col=0), dtype=np.float32)
        tmp.obs['sample_id'] = sample_id
        prop.append(tmp)

    prop = prop[0].concatenate(prop[1:])

    # Compute average props for visium
    vm_df = pd.DataFrame(prop.X, index=prop.obs.index, columns=prop.var.index)
    vm_df['sample_id'] = prop.obs['sample_id']
    vm_df = vm_df.groupby('sample_id').mean(1)
    vm_df = vm_df.melt(value_vars=vm_df.columns, ignore_index=False).reset_index()
    vm_df = vm_df.rename({'variable': 'cell_types', 'value': 'vm'}, axis=1)
    vm_df['vm'] /= vm_df.groupby('sample_id')['vm'].transform('sum') # normalize the vm column per sample_id
    return(vm_df)


def get_correlation(model_call, method, samples_oi, group, replace_zero=0.001):
    assert group in ["cell_types", "sample_id"]
    assert model_call in ["all", "condition", "lesion_type"]
    assert method in ["cellbender", "cellranger"]

    sc_prop = sc_proportions[method]
    vis_prop = vis_proportions[method][model_call]

    # only keep rows where sample_id is in samples_oi
    sc_prop = sc_prop.loc[sc_prop.sample_id.isin(samples_oi), :]
    vis_prop = vis_prop.loc[vis_prop.sample_id.isin(samples_oi), :]

    df = pd.merge(sc_prop, vis_prop, how='outer', on=['sample_id', 'cell_types'])
    df['sc'] = np.log10(df['sc'].replace({0: np.nan}))
    df['vm'] = np.log10(df['vm'].replace({0: np.nan}))
    df = df.fillna(np.log10(replace_zero))
    grouped_df = df.groupby(group)

    def pearson_corr(x):
        return x["sc"].corr(x["vm"], method="pearson")
    def spearman_corr(x):
        return x["sc"].corr(x["vm"], method="spearman")
    
    df_out = grouped_df.apply(pearson_corr).to_frame()
    df_out = df_out.rename({0: "Pearson_R"}, axis=1)
    df_out['Spearman_R'] = grouped_df.apply(spearman_corr)
    df_out["method"] = method
    df_out["model_call"] = model_call

    return df_out


def plot_correlation(model_call, method, samples_oi, replace_zero=0.001, save_path=None):
    sc_prop = sc_proportions[method]
    vis_prop = vis_proportions[method][model_call]

    # only keep rows where sample_id is in samples_oi
    sc_prop = sc_prop.loc[sc_prop.sample_id.isin(samples_oi), :]
    vis_prop = vis_prop.loc[vis_prop.sample_id.isin(samples_oi), :]

    df = pd.merge(sc_prop, vis_prop, how='outer', on=['sample_id', 'cell_types'])
    df['sc'] = np.log10(df['sc'].replace({0: np.nan}))
    df['vm'] = np.log10(df['vm'].replace({0: np.nan}))
    df = df.fillna(np.log10(replace_zero))

    # Plot corrs per cell type
    g = sns.lmplot(x='sc', y='vm', data=df, col='cell_types', col_wrap=4, height=3, facet_kws={"sharex": False, "sharey": False})
    g.map_dataframe(annotate)
    g.fig.set_facecolor('white')
    g.fig.suptitle(f'{method}_{model_call}_deconv_corr_celltype')
    g.fig.set_tight_layout(tight='pad')
    if save_path:
        g.fig.savefig(os.path.join(save_path, f'{method}_{model_call}_deconv_corr_celltype.pdf'), dpi=300)

    # Plot corrs per sample id
    g = sns.lmplot(x='sc', y='vm', data=df, col='sample_id', col_wrap=4, height=3, facet_kws={"sharex": False, "sharey": False})
    g.map_dataframe(annotate)
    g.fig.set_facecolor('white')
    g.fig.suptitle(f'{method}_{model_call}_deconv_corr_sample')
    g.fig.set_tight_layout(tight='pad')
    if save_path:
        g.fig.savefig(os.path.join(save_path, f'{method}_{model_call}_deconv_corr_sample.pdf'), dpi=300)

# set up the paths
current_folder = Path(__file__).parent
path_ann = current_folder / ".." / ".." / "data" / "prc" / "sc"
path_deconv = current_folder / ".." / ".." / "data" / "prc" / "vis" / "c2l_out"
output_path = current_folder / ".." / ".." / "out" / "deconvolution_eval"
output_path.mkdir(parents=True, exist_ok=True)

sample_meta = pd.read_excel(current_folder / ".." / ".." / "data" / "Metadata_all.xlsx", sheet_name="Visium")

sc_proportions = {
    "cellbender": get_sc_proportions(sc.read_h5ad(path_ann / f"annotated_cellbender_mod.h5ad").obs),
    "cellranger": get_sc_proportions(sc.read_h5ad(path_ann / f"annotated_cellranger.h5ad").obs)
}
#sc_proportions = {method: get_sc_proportions(sc.read_h5ad(path_ann / f"annotated_{method}_mod.h5ad").obs) for method in ["cellbender", "cellranger"]}

vis_proportions = {method: {model_call: get_vis_proportions(model_call, method, sample_meta, measure="abunds") for model_call in ["all", "condition", "lesion_type"]} for method in ["cellbender", "cellranger"]}

# check how well the c2l proportions and normalized abundances correlate
#props = get_vis_proportions("all", "cellranger", sample_meta, measure="props")
#abunds = get_vis_proportions("all", "cellranger", sample_meta, measure="abunds")
#fig, axs = plt.subplots(4, 6, figsize=(15, 10))
#axs = axs.flatten()
#for i, sample_id in enumerate(os.listdir(path_deconv / "cellranger")):
#    sns.scatterplot(x=props.loc[props.sample_id==sample_id].vm, y=abunds.loc[abunds.sample_id==sample_id].vm, ax=axs[i])
#    axs[i].plot([0, 1], [0, 1], transform=axs[i].transAxes, ls="--", c=".3")
#    axs[i].set_xlabel("Proportions")
#    axs[i].set_ylabel("Normalized Abundances")
#    axs[i].set_title(sample_id)
#    plt.title(sample_id)

# get sample_ids for which we have both sc and visium data
samples_oi = set(list(sc_proportions.values())[0].sample_id) & set(list(list(vis_proportions.values())[0].values())[0].sample_id)


for group in ["cell_types", "sample_id"]:
    for metric in ["Pearson_R", "Spearman_R"]:
        df_list = []
        for model_call in ["all", "condition", "lesion_type"]:
            for method in ["cellbender", "cellranger"]:
                df_list.append(get_correlation(model_call, method, samples_oi, group=group))
        df = pd.concat(df_list)
        df["identifer"] = df["method"] + "_" + df["model_call"]
        df.reset_index(inplace=True)

        g = sns.catplot(
            data=df,
            x=group, y=metric, hue="identifer",
            kind="bar", palette="dark", alpha=.6, height=6
        )
        g.despine(left=True)
        g.legend.set_title("")
        g.set_xticklabels(rotation=90)
        g.savefig(output_path / f"barplot_{group}_{metric}.pdf")


for method in ["cellbender", "cellranger"]:
    for model_call in ["all", "condition", "lesion_type"]:
        plot_correlation(model_call, method, samples_oi, save_path=output_path)
=======
    
# Plot corrs per sample
g = sns.lmplot(x='sc', y='vm', data=df, col='sample_id', col_wrap=4, height=3)
g.map_dataframe(annotate)
g.fig.set_facecolor('white')
g.fig.suptitle('Correlation per sample')
g.fig.set_tight_layout(tight='pad')
g.fig.savefig(os.path.join(output_path, 'deconv_corr_sample.pdf'), dpi=300)

# Plot corrs per cell type
g = sns.lmplot(x='sc', y='vm', data=df, col='leiden', col_wrap=4, height=3)
g.map_dataframe(annotate)
g.fig.set_facecolor('white')
g.fig.suptitle('Correlation per sample')
g.fig.set_tight_layout(tight='pad')
g.fig.savefig(os.path.join(output_path, 'deconv_corr_celltype.pdf'), dpi=300)

# Define params
n_celltypes = prop.var_names.shape[0]
n_rows = int(np.ceil(n_celltypes / 4))

# Plot proportions
for sample_id in sample_ids:
    # Read original slide
    slide = sc.read_visium(os.path.join(path_slides, '{0}/outs/'.format(sample_id)))
    tmp = prop[prop.obs['sample_id'] == sample_id]
    tmp.uns['spatial'] = slide.uns['spatial'].copy()
    tmp.obsm['spatial'] = slide.obsm['spatial'].copy()
    del slide
    
    # Plot props in slide
    fig, axes = plt.subplots(n_rows, 4, tight_layout=False, facecolor='white', figsize=(4*4, 3*n_rows))
    axes = axes.flatten()
    for celltype, ax in zip(prop.var_names, axes):
        sc.pl.spatial(tmp, color=celltype, size=2, frameon=False, cmap='Reds', return_fig=False, show=False, ax=ax)
    fig.suptitle('{0} proportions'.format(sample_id))
    fig.savefig(os.path.join(output_path, 'deconv_props_{0}.pdf'.format(sample_id)))
    del tmp
>>>>>>> 4059a4c8e2af9c39af787ebee1439fc854d311d6
