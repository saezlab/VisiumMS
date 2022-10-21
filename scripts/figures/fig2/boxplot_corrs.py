import pandas as pd
import numpy as np
import scanpy as sc
import os
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import scipy
import decoupler as dc


# Define path
fig_path = 'figures/manuscript/fig2/'
sfig_path = 'figures/manuscript/supp_fig2/'
fig_name = 'boxplot_corrs.pdf'
plt.rcParams['font.sans-serif'] = 'Arial'

# Read
adata = sc.read_h5ad('data/prc/visium/props_niches.h5ad')
adata = adata[adata.obs['leiden'] != '7']

def diag_melt(df):
    df = df.where(np.triu(np.ones(df.shape)).astype(bool))
    df = df.stack().reset_index()
    df['pair'] = df['level_0'] + '|' + df['level_1']
    del df['level_0']
    del df['level_1']
    df.columns = ['value', 'pair']
    df = df.loc[:, ['pair', 'value']]
    return df


def get_corrs_slide(slide):

    # Extract row and col names
    ctyps_props = slide.obsm['props'].columns.values
    ctyps_signs = slide.obsm['props'].columns.values

    # Init empty mats
    corrs = np.zeros((ctyps_signs.size, ctyps_props.size))
    pvals = np.zeros((ctyps_signs.size, ctyps_props.size))
    for i, c_a in enumerate(ctyps_signs):
        for j, c_b in enumerate(ctyps_props):
            
            if c_a == c_b:
                corrs[i, j], pvals[i, j] = np.nan, np.nan
            else:
                # Compute pearson
                corrs[i, j], pvals[i, j] = stats.pearsonr(slide.obsm['props'][c_a].values, slide.obsm['props'][c_b].values)

    # Transform to dfs
    corrs = pd.DataFrame(corrs, index=ctyps_signs, columns=ctyps_props)
    pvals = pd.DataFrame(pvals, index=ctyps_signs, columns=ctyps_props)

    corrs = diag_melt(corrs)

    return corrs


# Read meta
meta = pd.read_csv(os.path.join('data','metadata.csv')).set_index('sample_id')

def get_corrs_dfs(adata, lesion, meta):
    
    # Subset by lesion
    msk = meta['lesion_type'].values == lesion
    
    # Compute corrs
    ids = meta.index.values[msk]
    corrs = []
    for sample_id in ids:
        slide = adata[adata.obs['sample_id'] == sample_id]
        corr = get_corrs_slide(slide)
        corr['sample_id'] = sample_id
        corrs.append(corr)
    corrs = pd.concat(corrs)

    return corrs


corrs = []
for i, lesion in enumerate(['Control', 'Acute', 'Chronic Active', 'Chronic Inactive']):
    corr = get_corrs_dfs(adata, lesion, meta)
    corr['lesion_type'] = lesion
    corrs.append(corr)
corrs = pd.concat(corrs)

# Open significant interactions
vals = pd.read_csv('figures/manuscript/fig2/props_corrs_table.csv', index_col=0)
sign = pd.read_csv('figures/manuscript/fig2/props_pvals_table.csv', index_col=0)
vals = vals.loc[np.flip(vals.index)]
sign = sign.loc[np.flip(sign.index)]
vals = diag_melt(vals)
sign = diag_melt(sign)
sign['corrs'] = vals['value']
sign = sign[(sign['value'] < 0.10) & (sign['corrs'] > 0.15)]

# Add palette
palette = pd.read_csv('data/cond_palette.csv')
palette['rgb'] = [(r, g, b) for r, g, b in zip(palette['r'], palette['g'], palette['b'])]
palette = {k: v for k,v in zip(palette['cond'], palette['rgb'])}
corrs['lesion_type'] = pd.Categorical(corrs['lesion_type'].values, categories=palette.keys(), ordered=True)

# Test for significant leidens
pairs = np.unique(sign.pair.values)#sign.pair
corrs = corrs[np.isin(corrs['pair'].values, pairs)]
lesions = list(palette.keys())

# Test inside significant leidens
pair_lesions = [('Control', 'Acute'), ('Control', 'Chronic Active'), ('Control', 'Chronic Inactive')]
ranksum_table = []
for pair in pairs:
    vals = []
    for lesion in lesions:
        msk = (corrs['pair'] == pair) & (corrs['lesion_type'] == lesion)
        prop = corrs.loc[msk, 'value'].values
        vals.append(prop)
    
    rows = []
    for l_pair in pair_lesions:
        a, b = l_pair
        a, b = lesions.index(a), lesions.index(b)
        stat, pval = scipy.stats.ranksums(vals[a], vals[b])
        rows.append([pair, '{0}|{1}'.format(l_pair[0], l_pair[1]), stat, pval])
    adj_pvals = dc.p_adjust_fdr([x[-1] for x in rows])
    for i in range(len(pair_lesions)):
        rows[i].append(adj_pvals[i])
    ranksum_table.extend(rows)

ranksum_table = pd.DataFrame(ranksum_table, columns=['pair', 'l_pair', 'statistic', 'pvalue', 'adj_pvalue'])
ranksum_table.to_csv(os.path.join(fig_path, 'inters_ranksum_table.csv'), index=False)

# Plot sign props
sign_pairs = ranksum_table[ranksum_table['adj_pvalue'] < 0.10].sort_values('adj_pvalue')['pair'].values
g = (sns.catplot(col='pair', y='value', x='lesion_type', data=corrs[np.isin(corrs['pair'].values, sign_pairs)], col_wrap=2, kind='box', sharey=False, height=2, palette=palette)
 .set_titles(template='{col_name}').set_xlabels('').set_ylabels('Correlations').set_xticklabels(rotation=90)
)
os.makedirs(fig_path, exist_ok=True)
g.fig.savefig(os.path.join(fig_path, fig_name), bbox_inches='tight')

# Plot sign props
g = (sns.catplot(col='pair', y='value', x='lesion_type', data=corrs[~np.isin(corrs['pair'].values, sign_pairs)], col_wrap=10, kind='box', sharey=False, height=2, palette=palette)
 .set_titles(template='{col_name}').set_xlabels('').set_ylabels('Correlations').set_xticklabels(rotation=90)
)
os.makedirs(sfig_path, exist_ok=True)
g.fig.savefig(os.path.join(sfig_path, fig_name), bbox_inches='tight')
