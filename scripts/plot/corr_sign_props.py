import pandas as pd
import numpy as np
from scipy import stats
import os

from composition_stats import closure
from composition_stats import clr

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors


def get_corrs_slide(sample_id, sign_type):

    # Read props
    props = pd.read_csv('data/prc/visium/{0}/cell_props.csv'.format(sample_id), index_col=0)
    props.loc[:, :] = clr(closure(props.values))
    sign = pd.read_csv('data/prc/visium/{0}/cell_{1}_sign.csv'.format(sample_id, sign_type), index_col=0)

    # Intersect
    inter = props.index.intersection(sign.index)
    props = props.loc[inter]
    sign = sign.loc[inter]

    # Extract row and col names
    ctyps_props = props.columns.values
    ctyps_signs = sign.columns.values

    # Init empty mats
    corrs = np.zeros((ctyps_signs.size, ctyps_props.size))
    pvals = np.zeros((ctyps_signs.size, ctyps_props.size))
    for i, c_a in enumerate(ctyps_signs):
        for j, c_b in enumerate(ctyps_props):

            # Compute pearson
            corrs[i, j], pvals[i, j] = stats.pearsonr(sign[c_a].values, props[c_b].values)

    # Transform to dfs
    corrs = pd.DataFrame(corrs, index=ctyps_signs, columns=ctyps_props)
    pvals = pd.DataFrame(pvals, index=ctyps_signs, columns=ctyps_props)

    # Flip to have same order as misty
    corrs = corrs.loc[np.flip(corrs.index)]
    pvals = pvals.loc[np.flip(pvals.index)]

    return corrs, pvals


def aggregate(lst):
    out = np.zeros(lst[0].shape)
    for i in range(lst[0].shape[0]):
        for j in range(lst[0].shape[0]):
            vals = np.array([lst[k].iloc[i, j] for k in range(len(lst))])
            if np.all(~np.isfinite(vals)):
                out[i, j] = np.nan
            else:
                out[i, j] = np.mean(vals[np.isfinite(vals)])
    out = pd.DataFrame(out, index=lst[0].index, columns=lst[0].columns)
    return out


def plot_props_corrs(lesion_type, sign_type):
    # Read meta
    meta = pd.read_csv('data/metadata.csv').set_index('sample_id')
    meta = meta[meta['lesion_type'] == lesion_type]

    corrs = []
    pvals = []
    for sample_id in meta.index.values:
        corr, pval = get_corrs_slide(sample_id, sign_type)
        corrs.append(corr)
        pvals.append(pval)

    corrs = aggregate(corrs)
    pvals = aggregate(pvals)

    # Transform to asterisks
    pvals[np.isfinite(pvals)] = np.where((pvals < 0.05) & (np.abs(corrs) > 0.15), '*', '')

    # Define color map
    cmap = plt.get_cmap('coolwarm').copy()
    cmap.set_bad(color='gray')

    fig, ax = plt.subplots(1, 1, figsize=(4, 4), facecolor='white', dpi=125)
    htm = sns.heatmap(corrs, cmap=cmap, square=True, center=0, vmax=1, vmin=-1, ax=ax, cbar_kws={"shrink": .4, "aspect": 5},
                      annot=pvals.values.astype('U'), fmt='', annot_kws={'fontweight': 'black', 'color': 'black'})
    i = 0
    for _, spine in htm.spines.items():
        if i % 2 == 0:
            spine.set_visible(True)
        i += 1
    ax.set_xlabel('Predictor')
    ax.set_ylabel('Target')
    if sign_type == 'pos':
        sign_name = 'Disease'
    elif sign_type == 'neg':
        sign_name = 'Healthy'
    ax.set_title('Correlation {0} {1}'.format(sign_name, lesion_type), loc='left')
    os.makedirs('figures/corr_sign_prop/', exist_ok=True)
    fig.savefig('figures/corr_sign_prop/{0}_{1}.pdf'.format(sign_type, lesion_type.replace(' ', '')), bbox_inches='tight')

plot_props_corrs('Chronic Active', 'pos')
plot_props_corrs('Chronic Active', 'neg')
plot_props_corrs('Control', 'pos')
plot_props_corrs('Control', 'neg')
