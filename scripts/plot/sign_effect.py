import pandas as pd
import numpy as np
from scipy import stats
import os
import matplotlib.pyplot as plt
import seaborn as sns

meta = pd.read_csv('data/metadata.csv').set_index('sample_id')

# Read areas results
areas = []
for sample_id in meta.index:
    area = pd.read_csv('data/prc/visium/{0}/area.csv'.format(sample_id))
    areas.append(area)
areas = pd.concat(areas)

# Add condition info and format
meta = pd.read_csv('data/metadata.csv').set_index('sample_id')
areas['lesion_type'] = [meta.loc[s, 'lesion_type'] for s in areas['sample_id']]
areas['effect'] = areas['eff_area'] * areas['mean_act']
mean_areas = areas.groupby(['type', 'lesion_type', 'cell_type']).mean().reset_index()
areas = areas[np.isin(areas['lesion_type'], ['Control', 'Chronic Active'])]

# Plots
res_path = 'figures/sign/effect/'
os.makedirs(res_path, exist_ok=True)

# Area
g = (
    sns.catplot(row='type', col='cell_type', y='eff_area', x='lesion_type',
                data=areas, kind='box', sharey=False, height=2,
                palette={'Chronic Active': "#d4242c", 'Control': "#2a9d2a"})
    .set_titles(template='{col_name}')
    .set_xlabels('')
    .set_ylabels('Area')
    .set_xticklabels(rotation=45)
)
g.fig.subplots_adjust(top=0.8)
g.fig.suptitle('Signature area cross cell types')
g.fig.savefig('{0}/boxplot_area.pdf'.format(res_path), bbox_inches='tight')

# Act
g = (
    sns.catplot(row='type', col='cell_type', y='mean_act', x='lesion_type',
                data=areas, kind='box', sharey=False, height=2,
                palette={'Chronic Active': "#d4242c", 'Control': "#2a9d2a"})
    .set_titles(template='{col_name}')
    .set_xlabels('')
    .set_ylabels('Mean activity')
    .set_xticklabels(rotation=45)
)
g.fig.subplots_adjust(top=0.8)
g.fig.suptitle('Signature mean activity cross cell types')
g.fig.savefig('{0}/boxplot_act.pdf'.format(res_path), bbox_inches='tight')

# Efect
g = (
    sns.catplot(row='type', col='cell_type', y='effect', x='lesion_type',
                data=areas, kind='box', sharey=False, height=2,
                palette={'Chronic Active': "#d4242c", 'Control': "#2a9d2a"})
    .set_titles(template='{col_name}')
    .set_xlabels('')
    .set_ylabels('Effect')
    .set_xticklabels(rotation=45)
)
g.fig.subplots_adjust(top=0.8)
g.fig.suptitle('Signature effect cross cell types')
g.fig.savefig('{0}/boxplot_effect.pdf'.format(res_path), bbox_inches='tight')

def compute_t_test(areas, y):
    typs = np.array(['pos', 'neg'])
    cell_types = np.unique(areas['cell_type'].values)
    t_mat = np.zeros((typs.size, cell_types.size))
    p_mat = np.zeros((typs.size, cell_types.size))
    for i in range(typs.size):
        typ = typs[i]
        for j in range(cell_types.size):
            cell_type = cell_types[j]

            # Flip reference depending on typ
            if typ == 'pos':
                grt_cond, lsr_cond = 'Chronic Active', 'Control'
            elif typ == 'neg':
                grt_cond, lsr_cond = 'Control', 'Chronic Active'

            # Subset
            grt = areas[(areas['cell_type'] == cell_type) & (areas['type'] == typ) & (areas['lesion_type'] == grt_cond)][y].values
            lsr = areas[(areas['cell_type'] == cell_type) & (areas['type'] == typ) & (areas['lesion_type'] == lsr_cond)][y].values

            # Run one-sided t-test
            if grt.size >= 2 and lsr.size >= 2:
                t, p = stats.ttest_ind(grt, lsr, alternative='greater')
            else:
                t, p = 0., 1.

            # Store
            t_mat[i, j] = t
            p_mat[i, j] = p

    t_mat = pd.DataFrame(t_mat, index=typs, columns=cell_types)
    p_mat = pd.DataFrame(p_mat, index=typs, columns=cell_types)
    a_mat = p_mat.copy()
    a_mat = np.where((a_mat < 0.05) & (np.abs(t_mat) > 0), '*', '')

    return t_mat, a_mat


# Area
t_mat, a_mat = compute_t_test(areas, 'eff_area')
fig, ax = plt.subplots(1,1)
sns.heatmap(t_mat, annot=a_mat, fmt='', annot_kws={'fontweight': 'black', 'fontsize': 'xx-large'},
            cmap='coolwarm', square=True, center=0, cbar_kws={"shrink": .5, "aspect": 10},
            yticklabels=['Disease', 'Healthy'], ax=ax
           )
ax.tick_params(axis='x', rotation=45)
ax.tick_params(axis='y', rotation=0)
ax.set_title('T-value differences in signature area')
fig.savefig('{0}/heatmap_area.pdf'.format(res_path), bbox_inches='tight')

# Act
t_mat, a_mat = compute_t_test(areas, 'mean_act')
fig, ax = plt.subplots(1,1)
sns.heatmap(t_mat, annot=a_mat, fmt='', annot_kws={'fontweight': 'black', 'fontsize': 'xx-large'},
            cmap='coolwarm', square=True, center=0, cbar_kws={"shrink": .5, "aspect": 10},
            yticklabels=['Disease', 'Healthy'], ax=ax
           )
ax.tick_params(axis='x', rotation=45)
ax.tick_params(axis='y', rotation=0)
ax.set_title('T-value differences in signature mean activity')
fig.savefig('{0}/heatmap_act.pdf'.format(res_path), bbox_inches='tight')

# Effect
t_mat, a_mat = compute_t_test(areas, 'effect')
fig, ax = plt.subplots(1,1)
sns.heatmap(t_mat, annot=a_mat, fmt='', annot_kws={'fontweight': 'black', 'fontsize': 'xx-large'},
            cmap='coolwarm', square=True, center=0, cbar_kws={"shrink": .5, "aspect": 10},
            yticklabels=['Disease', 'Healthy'], ax=ax
           )
ax.tick_params(axis='x', rotation=45)
ax.tick_params(axis='y', rotation=0)
ax.set_title('T-value differences in signature effect')
fig.savefig('{0}/heatmap_effect.pdf'.format(res_path), bbox_inches='tight')
