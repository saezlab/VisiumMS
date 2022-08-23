import pandas as pd
import numpy as np
from scipy import stats
import os

import seaborn as sns
import matplotlib.pyplot as plt


def get_corrs_slide(sample_id):

    # Read
    props = pd.read_csv('data/prc/visium/{0}/cell_props.csv'.format(sample_id), index_col=0)
    min_prop = 1 / props.columns.size
    pos = pd.read_csv('data/prc/visium/{0}/cell_pos_sign.csv'.format(sample_id), index_col=0)
    neg = pd.read_csv('data/prc/visium/{0}/cell_neg_sign.csv'.format(sample_id), index_col=0)

    # Intersect
    inter = pos.index.intersection(neg.index).intersection(props.index)
    props = props.loc[inter]
    pos = pos.loc[inter]
    neg = neg.loc[inter]

    # Extract row and col names
    ctyps_signs = pos.columns.values

    # Init empty mats
    corrs = np.zeros((ctyps_signs.size, ))
    for i, c in enumerate(ctyps_signs):
        # Compute pearson
        msk = props[c] > min_prop
        corrs[i], _ = stats.pearsonr(pos[c].values[msk], neg[c].values[msk])
    corrs = pd.DataFrame([corrs], columns=ctyps_signs, index=[sample_id])

    return corrs


#def plot_sign_corrs(lesion_type):
lesion_type = 'Chronic Active'

# Read meta
meta = pd.read_csv('data/metadata.csv').set_index('sample_id')
meta = meta[meta['lesion_type'] == lesion_type]

corrs = []
for sample_id in meta.index.values:
    corr = get_corrs_slide(sample_id)
    corrs.append(corr)

corrs = pd.concat(corrs)
corrs = (
    corrs
    .reset_index()
    .melt(id_vars='index')
    .rename({'index': 'sample_id', 'variable': 'signature', 'value': 'corr'}, axis=1))

fig, ax = plt.subplots(1, 1, figsize=(4, 4), dpi=125, facecolor='white')
sns.boxplot(x='signature', y='corr', data=corrs)
ax.set_title('Correlation Disease-Healthy signs')
ax.set_ylabel('Correlation')
ax.set_xlabel('')
ax.tick_params(axis='x', rotation=45)
ax.axhline(y=0, linestyle='--', color="black")
fig.savefig('figures/corr_pos_neg.pdf', bbox_inches='tight')

