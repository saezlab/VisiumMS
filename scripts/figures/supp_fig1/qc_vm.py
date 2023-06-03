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
fig_name = 'qc_vm.pdf'

# Read data
adata = []
meta = pd.read_csv(os.path.join('data','metadata.csv')).set_index('sample_id')
for sample_id in meta.index.values:
    slide = read_slide(sample_id)
    slide.obs['sample_id'] = sample_id
    slide.obs['lesion_type'] = meta.loc[sample_id, 'lesion_type']
    adata.append(slide)
adata = adata[0].concatenate(adata[1:], join='outer')
plt.rcParams['font.sans-serif'] = 'Arial'

# Make palette
palette = pd.read_csv('data/cond_palette.csv')
palette['rgb'] = [(r, g, b) for r, g, b in zip(palette['r'], palette['g'], palette['b'])]
palette = {k: v for k,v in zip(palette['cond'], palette['rgb'])}
df = adata.obs.groupby('sample_id').count()[['n_genes']].rename({'n_genes': 'ncells'}, axis=1).reset_index()
df['ncells'] = np.log10(df['ncells'])
df['lesion'] = [meta.loc[s, 'lesion_type'] for s in df['sample_id']]
df['lesion'] = pd.Categorical(df['lesion'].values, categories=palette.keys(), ordered=True)
df = df.sort_values('lesion')
palette = {s_id: palette[les] for s_id, les in zip(df['sample_id'], df['lesion'])}

# Plot
fig, axes = plt.subplots(3, 1, figsize=(8,8), facecolor='white', sharex=True, tight_layout=True, dpi=150)

ax = axes[0]

sns.barplot(data=df, x="sample_id", y="ncells", ax=ax, palette=palette, order=palette.keys())
ax.set_ylabel('Number of spots (log10)')
ax.tick_params(axis='x', rotation=90)
ax.set_xlabel('')

ax = axes[1]
df = adata.obs[['sample_id', 'total_counts']].copy().reset_index()
df['log_total_counts'] = np.log10(df['total_counts'].values)
sns.boxplot(data=df, x="sample_id", y="log_total_counts", ax=ax, palette=palette, order=palette.keys(), fliersize=0)
ax.set_ylabel('Number of UMIs (log10)')
ax.tick_params(axis='x', rotation=90)
ax.set_xlabel('')

ax = axes[2]
df = adata.obs[['sample_id', 'n_genes_by_counts']].copy().reset_index()
df['log_n_genes_by_counts'] = np.log10(df['n_genes_by_counts'].values)
sns.boxplot(data=df, x="sample_id", y="log_n_genes_by_counts", ax=ax, palette=palette, order=palette.keys(), fliersize=0)
ax.set_ylabel('Number of genes (log10)')
ax.tick_params(axis='x', rotation=90)
ax.set_xlabel('')

os.makedirs(fig_path, exist_ok=True)
fig.savefig(os.path.join(fig_path, fig_name), bbox_inches='tight')
