import pandas as pd
import numpy as np
import scanpy as sc
import os
import seaborn as sns
import matplotlib.pyplot as plt
import decoupler as dc
from composition_stats import closure
import scipy


# Defina path
fig_path = 'figures/manuscript/fig2/'
sfig_path = 'figures/manuscript/supp_fig2/'
fig_name = 'props_vm.pdf'
plt.rcParams['font.sans-serif'] = 'Arial'
alpha = 0.075

# Extract proportions from visium slides
meta = pd.read_csv('data/metadata.csv').set_index('sample_id')
vs = []
for sample_id in meta.index:
    prop = pd.read_csv('data/prc/visium/{0}/cell_props.csv'.format(sample_id), index_col=0)
    cols = prop.columns
    prop = closure(prop.mean(axis=0).values)
    prop = pd.DataFrame(prop, columns=['props'], index=cols).reset_index().rename({'index': 'leiden'}, axis=1)
    prop['sample_id'] = sample_id
    prop['lesion_type'] = meta.loc[sample_id, 'lesion_type']
    vs.append(prop)
vs = pd.concat(vs)

# Add palette
palette = pd.read_csv('data/cond_palette.csv')
palette['rgb'] = [(r, g, b) for r, g, b in zip(palette['r'], palette['g'], palette['b'])]
palette = {k: v for k,v in zip(palette['cond'], palette['rgb'])}
vs['lesion_type'] = pd.Categorical(vs['lesion_type'].values, categories=palette.keys(), ordered=True)

# Test for significant leidens
leidens = np.unique(vs['leiden'].values)
lesions = list(palette.keys())
kruskal_table = []
for leiden in leidens:
    vals = []
    for lesion in lesions:
        msk = (vs['leiden'] == leiden) & (vs['lesion_type'] == lesion)
        prop = vs.loc[msk, 'props'].values
        vals.append(prop)
    stat, pval = scipy.stats.kruskal(*vals)
    kruskal_table.append([leiden, stat, pval])
kruskal_table = pd.DataFrame(kruskal_table, columns=['leiden', 'statistic', 'pvalue'])
kruskal_table['adj_pvalue'] = dc.p_adjust_fdr(kruskal_table['pvalue'])
kruskal_table.to_csv(os.path.join(fig_path, 'vm_kruskal_table.csv'), index=False)

# Subset by alpha
leidens = kruskal_table[kruskal_table['adj_pvalue'] < alpha]['leiden'].values

# Plot sign props
g = (sns.catplot(col='leiden', y='props', x='lesion_type', data=vs[np.isin(vs['leiden'].values, leidens)], col_wrap=3, kind='box', sharey=False, height=2, palette=palette)
 .set_titles(template='{col_name}').set_xlabels('').set_ylabels('Proportions').set_xticklabels(rotation=90)
)
os.makedirs(fig_path, exist_ok=True)
g.fig.savefig(os.path.join(fig_path, fig_name), bbox_inches='tight')

# Plot unsign props
g = (sns.catplot(col='leiden', y='props', x='lesion_type', data=vs[~np.isin(vs['leiden'].values, leidens)], col_wrap=3, kind='box', sharey=False, height=2, palette=palette)
 .set_titles(template='{col_name}').set_xlabels('').set_ylabels('Proportions').set_xticklabels(rotation=90)
)
os.makedirs(sfig_path, exist_ok=True)
g.fig.savefig(os.path.join(sfig_path, fig_name), bbox_inches='tight')

# Test inside significant leidens
pairs = [('Control', 'Acute'), ('Control', 'Chronic Active'), ('Control', 'Chronic Inactive')]
ranksum_table = []
for leiden in leidens:
    vals = []
    for lesion in lesions:
        msk = (vs['leiden'] == leiden) & (vs['lesion_type'] == lesion)
        prop = vs.loc[msk, 'props'].values
        vals.append(prop)
    
    for pair in pairs:
        a, b = pair
        a, b = lesions.index(a), lesions.index(b)
        stat, pval = scipy.stats.ranksums(vals[a], vals[b])
        ranksum_table.append([leiden, '{0}|{1}'.format(pair[0], pair[1]), stat, pval])
ranksum_table = pd.DataFrame(ranksum_table, columns=['leiden', 'pair', 'statistic', 'pvalue'])
ranksum_table['adj_pvalue'] = dc.p_adjust_fdr(ranksum_table['pvalue'])
ranksum_table.to_csv(os.path.join(fig_path, 'vm_ranksum_table.csv'), index=False)
