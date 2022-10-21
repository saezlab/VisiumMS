import pandas as pd
import numpy as np
import scanpy as sc
import os
import seaborn as sns
import matplotlib.pyplot as plt
import decoupler as dc
from composition_stats import closure
import scipy


# Define path
fig_path = 'figures/manuscript/fig2/'
sfig_path = 'figures/manuscript/supp_fig2/'
fig_name = 'props_niches.pdf'
plt.rcParams['font.sans-serif'] = 'Arial'
alpha = 0.05

# Read
adata = sc.read_h5ad('data/prc/visium/props_niches.h5ad')
adata = adata[adata.obs['leiden'] != '7']

# Open metadata
meta = pd.read_csv(os.path.join('data','metadata.csv')).set_index('sample_id')

# Compute props
df_props = adata.obs.groupby(['leiden', 'sample_id']).count()[['batch']].reset_index()
df_props['leiden'] = df_props['leiden'].astype(str)
df_props['lesion_type'] = [meta.loc[s, 'lesion_type'] for s in df_props['sample_id']]
df_props['prop'] = df_props['batch'] / df_props.groupby('sample_id')['batch'].transform('sum')
sample_ids = np.unique(df_props['sample_id'].values)
for sample_id in sample_ids:
    msk = df_props['sample_id'] == sample_id
    df_props.loc[msk, 'prop'] = closure(df_props['prop'].values[msk])

# Add palette
palette = pd.read_csv('data/cond_palette.csv')
palette['rgb'] = [(r, g, b) for r, g, b in zip(palette['r'], palette['g'], palette['b'])]
palette = {k: v for k,v in zip(palette['cond'], palette['rgb'])}
df_props['lesion_type'] = pd.Categorical(df_props['lesion_type'].values, categories=palette.keys(), ordered=True)

# Test for significant leidens
leidens = np.unique(df_props['leiden'].values)
lesions = list(palette.keys())
kruskal_table = []
for leiden in leidens:
    vals = []
    for lesion in lesions:
        msk = (df_props['leiden'] == leiden) & (df_props['lesion_type'] == lesion)
        prop = df_props.loc[msk, 'prop'].values
        vals.append(prop)
    stat, pval = scipy.stats.kruskal(*vals)
    kruskal_table.append([leiden, stat, pval])
kruskal_table = pd.DataFrame(kruskal_table, columns=['leiden', 'statistic', 'pvalue'])
kruskal_table['adj_pvalue'] = dc.p_adjust_fdr(kruskal_table['pvalue'])
kruskal_table.to_csv(os.path.join(fig_path, 'niches_kruskal_table.csv'), index=False)

# Subset by alpha
leidens = kruskal_table[kruskal_table['adj_pvalue'] < alpha]['leiden'].values

# Plot sign props
g = (sns.catplot(x="lesion_type", col="leiden", y='prop', data=df_props[np.isin(df_props['leiden'].values, leidens)], col_wrap=2, kind='box', sharey=False, height=2, palette=palette)
 .set_titles(template='Niche {col_name}').set_xlabels('').set_ylabels('Proportions').set_xticklabels(rotation=90)
)
os.makedirs(fig_path, exist_ok=True)
g.fig.savefig(os.path.join(fig_path, fig_name), bbox_inches='tight')

# Plot unsign props
g = (sns.catplot(x="lesion_type", col="leiden", y='prop', data=df_props[~np.isin(df_props['leiden'].values, leidens)], col_wrap=4, kind='box', sharey=False, height=2, palette=palette)
 .set_titles(template='Niche {col_name}').set_xlabels('').set_ylabels('Proportions').set_xticklabels(rotation=90)
)
os.makedirs(sfig_path, exist_ok=True)
g.fig.savefig(os.path.join(sfig_path, fig_name), bbox_inches='tight')

# Test inside significant leidens
pairs = [('Control', 'Acute'), ('Control', 'Chronic Active'), ('Control', 'Chronic Inactive')]
ranksum_table = []
for leiden in leidens:
    vals = []
    for lesion in lesions:
        msk = (df_props['leiden'] == leiden) & (df_props['lesion_type'] == lesion)
        prop = df_props.loc[msk, 'prop'].values
        vals.append(prop)
    
    for pair in pairs:
        a, b = pair
        a, b = lesions.index(a), lesions.index(b)
        stat, pval = scipy.stats.ranksums(vals[a], vals[b])
        ranksum_table.append([leiden, '{0}|{1}'.format(pair[0], pair[1]), stat, pval])
ranksum_table = pd.DataFrame(ranksum_table, columns=['leiden', 'pair', 'statistic', 'pvalue'])
ranksum_table['adj_pvalue'] = dc.p_adjust_fdr(ranksum_table['pvalue'])
ranksum_table.to_csv(os.path.join(fig_path, 'niches_ranksum_table.csv'), index=False)
