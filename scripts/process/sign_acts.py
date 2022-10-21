import pandas as pd
import numpy as np
import scanpy as sc
import decoupler as dc
import os


# Read and concat
pos = pd.read_csv('data/prc/sign/deg/pos.csv')
neg = pd.read_csv('data/prc/sign/deg/neg.csv')

# Generate stat
pos['stat'] = -np.log10(pos['pvals']) * pos['logFCs']
neg['stat'] = -np.log10(neg['pvals']) * neg['logFCs']

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

    return slide


def compute_signature(slide, net, name):

    # Compute acts
    dc.run_ulm(slide, net, source='contrast', target='name', weight='stat', use_raw=False)
    slide.obsm[name] = slide.obsm['ulm_estimate'].copy()


def compute_area(slide, obsm_key, thr=2):

    # Min prop to consider composition
    min_prop = 1 / slide.obsm['props'].shape[1]

    # Init empty df
    df = []
    for cell_type in slide.obsm[obsm_key]:

        # Extract acts and filter by composition
        acts = slide.obsm[obsm_key][cell_type].values
        prop = slide.obsm['props'][cell_type].values
        prop_msk = prop > min_prop
        acts = acts[prop_msk]

        # Select effective spots by activity
        msk = acts > thr

        # Compute eff area and its mean activity
        area = msk.sum() / msk.size
        mean = np.mean(acts[msk])
        if np.isnan(mean):
            mean = 0.

        # Append
        df.append([cell_type, area, mean])

    # Format as df
    df = pd.DataFrame(df, columns = ['cell_type', 'eff_area', 'mean_act'])
    df['type'] = obsm_key

    return df

meta = pd.read_csv('data/metadata.csv').set_index('sample_id')

for sample_id in meta.index:

    # Load slide
    print(sample_id)
    slide = read_slide(sample_id)

    # Compute signatures
    compute_signature(slide, pos, 'pos')
    compute_signature(slide, neg, 'neg')

    # Compute areas
    p_area = compute_area(slide, 'pos')
    n_area = compute_area(slide, 'neg')
    area = pd.concat([p_area, n_area])
    area['sample_id'] = sample_id

    # Write
    slide.obsm['pos'].to_csv('data/prc/visium/{0}/cell_pos_sign.csv'.format(sample_id))
    slide.obsm['neg'].to_csv('data/prc/visium/{0}/cell_neg_sign.csv'.format(sample_id))
    area.to_csv('data/prc/visium/{0}/area.csv'.format(sample_id), index=False)
