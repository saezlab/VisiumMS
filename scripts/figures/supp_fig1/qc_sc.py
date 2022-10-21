import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns


# Defina path
fig_path = 'figures/manuscript/supp_fig1/'
fig_name = 'qc_sc.pdf'

# Read data
adata = sc.read_h5ad('data/prc/sc/annotated.h5ad')
plt.rcParams['font.sans-serif'] = 'Arial'
meta = pd.read_csv(os.path.join('data','metadata.csv')).set_index('sample_id')

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
df = df.sort_values('lesion')
sns.barplot(data=df, x="sample_id", y="ncells", ax=ax, palette=palette, order=palette.keys())
ax.set_ylabel('Number of cells (log10)')
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
