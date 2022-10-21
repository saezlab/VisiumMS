import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns


def read_slide(sample_id):

    # Read rna-seq
    slide = sc.read_visium('data/raw/visium/{0}/outs/'.format(sample_id))
    slide.var_names_make_unique()
    sc.pp.filter_genes(slide, min_cells=3)
    sc.pp.filter_cells(slide, min_genes=200)

    # QC
    slide.var['mt'] = slide.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(slide, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # Store raw counts
    slide.layers['counts'] = slide.X.copy()

    # Normalize
    sc.pp.normalize_total(slide, target_sum=1e4)
    sc.pp.log1p(slide)

    # Read props
    props = pd.read_csv('data/prc/visium/{0}/cell_props.csv'.format(sample_id), index_col=0)
    inter = slide.obs.index.intersection(props.index)
    slide.obsm['props'] = props.loc[inter]

    return slide


# Defina path
fig_path = 'figures/manuscript/supp_fig1/'
fig_name = 'spatial_{0}_{1}_ngenes.pdf'

# Read data
cols = [None, 'log_total_counts', 'log_n_genes_by_counts']
meta = pd.read_csv(os.path.join('data','metadata.csv')).set_index('sample_id')
meta = meta.sort_values('lesion_type')

for i, sample_id in enumerate(meta.index.values):
    slide = read_slide(sample_id)
    slide.obs['sample_id'] = sample_id
    lesion = meta.loc[sample_id, 'lesion_type']
    slide.obs['lesion_type'] = lesion
    slide.obs = slide.obs.assign(log_n_genes_by_counts=lambda x: np.log10(x['n_genes_by_counts']))
    slide.obs = slide.obs.assign(log_total_counts=lambda x: np.log10(x['total_counts']))
    fig, axes = plt.subplots(1, 3, figsize=(10, 3), facecolor='white', dpi=150)
    for j, col in enumerate(cols):
        ax = axes[j]
        sc.pl.spatial(slide, color=col, size=2, ax=ax, return_fig=False, show=False, frameon=False)
        if i != 0:
            ax.set_title('')
    axes[0].set_title('H&E staining')
    axes[1].set_title('Number of UMIs (log10)')
    axes[2].set_title('Number of genes (log10)')
    os.makedirs(fig_path, exist_ok=True)
    fig.savefig(os.path.join(fig_path, fig_name.format(lesion, sample_id)), bbox_inches='tight')
