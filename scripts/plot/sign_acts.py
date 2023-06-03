import pandas as pd
import numpy as np
import scanpy as sc
import os
import matplotlib.pyplot as plt


def read_slide(sample_id):

    # Read rna-seq
    slide = sc.read_visium('data/raw/visium/{0}/outs/'.format(sample_id))
    slide.var_names_make_unique()
    sc.pp.filter_genes(slide, min_cells=3)
    sc.pp.filter_cells(slide, min_genes=200)

    # Store raw counts
    slide.layers['counts'] = slide.X

    # Normalize
    sc.pp.normalize_total(slide, target_sum=1e4)
    sc.pp.log1p(slide)

    # Read props
    props = pd.read_csv('data/prc/visium/{0}/cell_props.csv'.format(sample_id), index_col=0)
    inter = slide.obs.index.intersection(props.index)
    slide.obsm['props'] = props.loc[inter]

    # Read signs
    pos = pd.read_csv('data/prc/visium/{0}/cell_pos_sign.csv'.format(sample_id), index_col=0)
    inter = slide.obs.index.intersection(pos.index)
    slide.obsm['pos'] = pos.loc[inter]
    neg = pd.read_csv('data/prc/visium/{0}/cell_neg_sign.csv'.format(sample_id), index_col=0)
    inter = slide.obs.index.intersection(neg.index)
    slide.obsm['neg'] = neg.loc[inter]

    return slide


def get_area_obs(slide, cell_type, key, min_prop, thr=2):
    msk_prop = slide.obsm['props'][cell_type].values > min_prop
    msk_acts = slide.obsm[key][cell_type].values > thr

    area = np.array([np.nan] * msk_prop.size, dtype=object)
    area[msk_prop] = 'Total area'
    area[msk_prop * msk_acts] = 'Effectve area'
    slide.obs[key + '_area'] = area


def plot_signs_slide(sample_id):

    # Read slide
    slide = read_slide(sample_id)

    # Min prop to consider composition
    min_prop = 1 / slide.obsm['props'].shape[1]

    cell_types = np.unique(slide.obsm['pos'].columns)
    fig, axes = plt.subplots(2, cell_types.size, facecolor='white', dpi=150)
    for i, cell_type in enumerate(cell_types):
        get_area_obs(slide, cell_type, 'pos', min_prop, thr=2)
        get_area_obs(slide, cell_type, 'neg', min_prop, thr=2)

        palette={"Total area": "gainsboro", "Effectve area": "#d4242c"}
        sc.pl.spatial(slide, color='pos_area', palette=palette,
                  ax=axes[0, i], return_fig=False,
                  show=False, frameon=False,
                  na_color=None, size=1.75, na_in_legend=False,
                 )
        axes[0, i].set_title(cell_type)
        palette={"Total area": "gainsboro", "Effectve area": "#2a9d2a"}
        sc.pl.spatial(slide, color='neg_area', palette=palette,
                  ax=axes[1, i], return_fig=False,
                  show=False, frameon=False,
                  na_color=None, size=1.75, na_in_legend=False,
                 )
        axes[1, i].set_title('')
        axes[0, i].legend().remove()
        axes[1, i].legend().remove()

    # Save
    res_path = 'figures/sign/slides'
    os.makedirs(res_path, exist_ok=True)
    fig.subplots_adjust(wspace=0.1, hspace=-0.5)
    fig.savefig('{0}/{1}.pdf'.format(res_path, sample_id))


meta = pd.read_csv('data/metadata.csv').set_index('sample_id')

for sample_id in meta.index:
    print(sample_id)
    plot_signs_slide(sample_id)
