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
fig_path = 'figures/manuscript/supp_fig2/'
fig_name = 'props_sc.pdf'
plt.rcParams['font.sans-serif'] = 'Arial'

# Read raw sc object
adata = sc.read('data/prc/sc/raw.h5ad')

# Compute props for sc
sc = (adata.obs
 .groupby(['sample_id', 'leiden'])
 .count()
 .reset_index()
 [['sample_id', 'leiden', 'total_counts']]
)
sc['props'] = (sc['total_counts'] /
               sc.groupby(['sample_id'])['total_counts']
               .transform('sum'))

# Add lesion type labels 
meta = pd.read_csv('data/metadata.csv').set_index('sample_id')
sc['lesion_type'] = [meta.loc[s, 'lesion_type'] for s in sc['sample_id']]

# Add palette
palette = pd.read_csv('data/cond_palette.csv')
palette['rgb'] = [(r, g, b) for r, g, b in zip(palette['r'], palette['g'], palette['b'])]
palette = {k: v for k,v in zip(palette['cond'], palette['rgb'])}

# Plot
g = (sns.catplot(col='leiden', y='props', x='lesion_type', data=sc, col_wrap=6, kind='box', sharey=False, height=2, palette=palette, order=palette.keys())
 .set_titles(template='{col_name}').set_xlabels('').set_ylabels('Proportions').set_xticklabels(rotation=90)
)
os.makedirs(fig_path, exist_ok=True)
g.fig.savefig(os.path.join(fig_path, fig_name), bbox_inches='tight')

# Test for significant leidens
leidens = np.unique(sc['leiden'].values)
lesions = list(palette.keys())
kruskal_table = []
for leiden in leidens:
    vals = []
    for lesion in lesions:
        msk = (sc['leiden'] == leiden) & (sc['lesion_type'] == lesion)
        prop = sc.loc[msk, 'props'].values
        vals.append(prop)
    stat, pval = scipy.stats.kruskal(*vals)
    kruskal_table.append([leiden, stat, pval])
kruskal_table = pd.DataFrame(kruskal_table, columns=['leiden', 'statistic', 'pvalue'])
kruskal_table['adj_pvalue'] = dc.p_adjust_fdr(kruskal_table['pvalue'])
kruskal_table.to_csv(os.path.join(fig_path, 'kruskal_table.csv'), index=False)

# Test inside significant leidens
pairs = [('Control', 'Acute'), ('Control', 'Chronic Active'), ('Control', 'Chronic Inactive')]
ranksum_table = []
for leiden in leidens:
    vals = []
    for lesion in lesions:
        msk = (sc['leiden'] == leiden) & (sc['lesion_type'] == lesion)
        prop = sc.loc[msk, 'props'].values
        vals.append(prop)
    
    for pair in pairs:
        a, b = pair
        a, b = lesions.index(a), lesions.index(b)
        stat, pval = scipy.stats.ranksums(vals[a], vals[b])
        ranksum_table.append([leiden, '{0}|{1}'.format(pair[0], pair[1]), stat, pval])
ranksum_table = pd.DataFrame(ranksum_table, columns=['leiden', 'pair', 'statistic', 'pvalue'])
ranksum_table['adj_pvalue'] = dc.p_adjust_fdr(ranksum_table['pvalue'])
ranksum_table.to_csv(os.path.join(fig_path, 'ranksum_table.csv'), index=False)
