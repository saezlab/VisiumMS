import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import scipy
import decoupler as dc


def compute_pairwise_test(df, col='region'):
    if col == 'region':
        pairs = [['LC', 'LR'], ['LC', 'NAWM'], ['LR', 'NAWM']]
    else:
        pairs = [['control', 'chronic_active'], ['control', 'chronic_inactive'],
                 ['chronic_active', 'chronic_inactive']]

    pair_test = []
    for pair in pairs:
        a, b = pair
        a = df[df[col] == a]['normalized_spots_area'].values
        b = df[df[col] == b]['normalized_spots_area'].values
        if a.size > 0 and b.size > 0:
            stat, pval = scipy.stats.ranksums(a, b)
            row = ['|'.join(pair), stat, pval]
            pair_test.append(row)

    if len(pair_test) > 0:
        adj_pvals = dc.p_adjust_fdr([r[-1] for r in pair_test])
        for i in range(len(adj_pvals)):
            pair_test[i].append(adj_pvals[i])

        pair_test = pd.DataFrame(pair_test, columns=['pair', 'statistic', 'pvalue', 'adj_pvalue'])
        return pair_test


def plot_boxplots(gene, micro):
    palette = {
    'LC': '#D047B1',
    'LR': '#F9C823',
    'NAWM': '#7DC581',
    }
    pair_tests = []
    
    fig, axes = plt.subplots(1, 3, figsize=(6,3), sharey=True, facecolor='white')
    max_n = np.max(micro[micro['transcript'] == gene]['normalized_spots_area'].values)
    max_n = max_n * 1.3

    ax = axes[0]
    lesion_type = 'control'
    df = micro[(micro['transcript'] == gene) & (micro['subgroup'] == lesion_type)]
    pair_test = compute_pairwise_test(df)
    if pair_test is not None:
        pair_test['lesion_type'] = lesion_type
        pair_tests.append(pair_test)
    sns.boxplot(x='region', y='normalized_spots_area', data=df, ax=ax, palette=palette, order=palette.keys())
    sns.swarmplot(x='region', y='normalized_spots_area', data=df, ax=ax, palette=palette, linewidth=1, order=palette.keys())
    ax.set_xlabel('Control')
    ax.set_ylabel('Normalized spots per area')
    ax.set_axisbelow(True)
    ax.set_xticks([])

    ax = axes[1]
    lesion_type = 'chronic_active'
    df = micro[(micro['transcript'] == gene) & (micro['subgroup'] == lesion_type)]
    pair_test = compute_pairwise_test(df)
    if pair_test is not None:
        pair_test['lesion_type'] = lesion_type
        pair_tests.append(pair_test)
    sns.boxplot(x='region', y='normalized_spots_area', data=df, ax=ax, palette=palette, order=palette.keys())
    sns.swarmplot(x='region', y='normalized_spots_area', data=df, ax=ax, palette=palette, order=palette.keys(), linewidth=1)
    ax.set_ylabel('')
    ax.set_xlabel('Chronic Active')
    ax.set_title(gene)
    ax.yaxis.set_ticks_position('none')
    ax.set_axisbelow(True)
    ax.set_xticks([])

    ax = axes[2]
    lesion_type = 'chronic_inactive'
    df = micro[(micro['transcript'] == gene) & (micro['subgroup'] == lesion_type)]
    pair_test = compute_pairwise_test(df)
    if pair_test is not None:
        pair_test['lesion_type'] = lesion_type
        pair_tests.append(pair_test)
    sns.boxplot(x='region', y='normalized_spots_area', data=df, ax=ax, palette=palette, order=palette.keys())
    sns.swarmplot(x='region', y='normalized_spots_area', data=df, ax=ax, palette=palette, order=palette.keys(), linewidth=1)
    ax.set_ylabel('')
    ax.set_xlabel('Chronic Inactive')
    ax.yaxis.set_ticks_position('none') 
    ax.set_axisbelow(True)
    ax.set_xticks([])
    ax.set_ylim(-0.0005, max_n)

    fig.subplots_adjust(wspace=0.0, hspace=None)
    
    pair_tests = pd.concat(pair_tests)
    pair_tests['transcript'] = gene
    pair_tests = pair_tests[['transcript', 'lesion_type', 'pair', 'statistic', 'pvalue', 'adj_pvalue']]
    
    df = micro[micro['transcript'] == gene]
    lesion_test = compute_pairwise_test(df, col='subgroup')
    lesion_test['transcript'] = gene
    lesion_test = lesion_test[['transcript', 'pair', 'statistic', 'pvalue', 'adj_pvalue']]
    
    return fig, pair_tests, lesion_test


# Defina path
fig_path = 'figures/manuscript/fig3/'
plt.rcParams['font.sans-serif'] = 'Arial'

# Read quantifications
df = pd.read_csv('data/quantification.tsv', sep='\t')
df = df[df['stainings'] == 'FTL']

# Iterate
regions = []
lesions = []
transcripts = np.unique(df['transcript'].values)
for gene in transcripts:

    # Compute plots and tables
    fig, region, lesion = plot_boxplots(gene, df)
    regions.append(region)
    lesions.append(lesion)

    # Write
    os.makedirs(fig_path, exist_ok=True)
    fig.savefig(os.path.join(fig_path, 'boxplot_{0}.pdf'.format(gene)), bbox_inches='tight')

regions = pd.concat(regions)
lesions = pd.concat(lesions)

# Write
regions.to_csv(os.path.join(fig_path, 'quant_regions.csv'))
lesions.to_csv(os.path.join(fig_path, 'quant_lesions.csv'))
