import pandas as pd
import numpy as np
import scanpy as sc
import os
import matplotlib.pyplot as plt
import decoupler as dc
import scipy
import seaborn as sns


# Define path
fig_path = 'figures/manuscript/supp_fig2/'
fig_name = 'violins_pwbysex.pdf'
plt.rcParams['font.sans-serif'] = 'Arial'

# Read
adata = sc.read_h5ad('data/prc/visium/props_niches.h5ad')
adata = adata[adata.obs['leiden'] != '7']
meta = pd.read_csv(os.path.join('data','metadata.csv')).set_index('sample_id')

adata.obs['sex'] = [meta.loc[s, 'sex'] for s in adata.obs['sample_id']]
acts = dc.get_acts(adata, 'progeny')
ol_leidens = ['0', '2', '3', '4']
acts.obs['group'] = ['OL' if x in ol_leidens else 'AS' for x in acts.obs['leiden']]
df = acts.obs[['group', 'sex']].copy()
df = pd.concat([df, acts.obsm['progeny']], axis=1)

def test_diff(acts, pw, leidens):
    isin = acts[np.isin(acts.obs['leiden'], leidens)][:, pw].X.ravel()
    ntin = acts[~np.isin(acts.obs['leiden'], leidens)][:, pw].X.ravel()
    s, p = scipy.stats.ranksums(isin, ntin)
    return s, p

table = []
for sex in ['M', 'F']:
    for pw in ['Androgen', 'Estrogen']:
        s, p = test_diff(acts[acts.obs['sex'] == sex], pw, ol_leidens)
        table.append([pw, sex, s, p])
table = pd.DataFrame(table, columns=['pathway', 'sex', 'statistic', 'pvalue'])
table['adj_pvalue'] = dc.p_adjust_fdr(table['pvalue'].values)
os.makedirs(fig_path, exist_ok=True)
table.to_csv(os.path.join(fig_path, 'pwbysex.csv'), index=False)

fig, ax = plt.subplots(1,2, figsize=(8,3), facecolor='white')
sns.violinplot(x='sex', y='Androgen', hue='group', data=df, ax=ax[0])
ax[0].get_legend().remove()
ax[0].set_xlabel('')
ax[0].set_ylabel('Activity')
ax[0].set_title('Androgen')
sns.violinplot(x='sex', y='Estrogen', hue='group', data=df, ax=ax[1])
ax[1].set_xlabel('')
ax[1].set_ylabel('')
ax[1].set_title('Estrogen')
handles, labels = ax[1].get_legend_handles_labels()
ax[1].legend(handles, labels, bbox_to_anchor=(1,0.5), loc='center left', frameon=False)
fig.savefig(os.path.join(fig_path, 'pwbysex.pdf'), bbox_inches='tight')
