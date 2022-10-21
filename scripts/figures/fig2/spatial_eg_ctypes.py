import pandas as pd
import numpy as np
import scanpy as sc
import os
import matplotlib.pyplot as plt


# Define path
fig_path = 'figures/manuscript/fig2/'
fig_name = 'spatial_eg_ctypes.pdf'
plt.rcParams['font.sans-serif'] = 'Arial'

def read_slide(sample_id):

    # Read rna-seq
    slide = sc.read_visium('data/raw/visium/{0}/outs/'.format(sample_id))
    slide.var_names_make_unique()
    sc.pp.filter_genes(slide, min_cells=3)
    sc.pp.filter_cells(slide, min_genes=200)

    # Store raw counts
    slide.layers['counts'] = slide.X.copy()

    # Normalize
    sc.pp.normalize_total(slide, target_sum=1e4)
    sc.pp.log1p(slide)

    # Read props
    props = pd.read_csv('data/prc/visium/{0}/cell_props.csv'.format(sample_id), index_col=0)
    inter = slide.obs.index.intersection(props.index)
    slide.obsm['props'] = props.loc[inter]

    # Annotate by max celltype
    ctype_dict = {
    'Astros': 'Astros',
    'Oligos': 'Oligos',
    'Astros_c': 'Astros',
    'Endothelia': 'Other',
    'Macrophages_f': 'Macro-Microglia',
    'Microglia': 'Macro-Microglia',
    'OPC': 'Oligos',
    'Stroma': 'Other',
    'T-cells': 'Other',
    'B-cells': 'Other',
    'Oligos_d': 'Oligos'
    }

    df = slide.obsm['props'].reset_index().melt(id_vars='index')
    idx = df.groupby(['index'])['value'].transform(max) == df['value']
    df = df[idx].set_index('index')
    df['variable'] = [ctype_dict[i] for i in df['variable']]
    slide.obs['cell_type'] = df['variable']

    return slide

meta = pd.read_csv(os.path.join('data','metadata.csv')).set_index('sample_id')

# Define palette
ctype_palette = {
    'Oligos': '#7DC581',
    'Astros': '#D047B1',
    'Macro-Microglia': '#F9C823',
    'Other': 'gray'
}
samples = ['CO85', 'MS371', 'MS411', 'MS497I']

# Plot
fig, axes = plt.subplots(2, 2, figsize=(6, 6), dpi=150, facecolor='white')
axes = axes.ravel()

for i, sample in enumerate(samples):
    ax = axes[i]
    slide = read_slide(sample)
    sc.pl.spatial(slide, color='cell_type', size=1.6, ax=ax, return_fig=False, show=False, frameon=False, palette=ctype_palette, legend_loc=False)
    lesion = meta.loc[sample, 'lesion_type']
    ax.set_title('{0} {1}'.format(lesion, sample))

# Write
os.makedirs(fig_path, exist_ok=True)
fig.savefig(os.path.join(fig_path, fig_name), bbox_inches='tight')
