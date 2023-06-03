import pandas as pd
import numpy as np
import scanpy as sc
import decoupler as dc
import os

path_raw = 'data/prc/sc/raw.h5ad'
path_slides = 'data/prc/visium/'

# Load raw object
adata = sc.read(path_raw)

# Get pseudo-bulk profile
padata = dc.get_pseudobulk(adata, sample_col='sample_id', groups_col='cell_type',
                           min_prop=0.2, min_smpls=3, use_raw=False)

# Normalize
sc.pp.normalize_total(padata, target_sum=1e4)
sc.pp.log1p(padata)

# Run contrast
logFCs, pvals = dc.get_contrast(padata,
                                group_col='cell_type',
                                condition_col='lesion_type',
                                condition='Chronic Active',
                                reference='Control',
                                method='t-test'
                               )

# Extract deg
deg = dc.format_contrast_results(logFCs, pvals)

# Filter by basic thrs
deg = deg[(np.abs(deg['logFC']) > 0.5) & (deg['pval'] < 0.05)]

# Filter MT genes
deg = deg[[not g.startswith('MT-') for g in deg['name']]]

# Filter genes that are sign in other cell_types
counts = deg.groupby('name').count()
genes = counts.index[counts['contrast'] == 1].values
deg = deg[np.isin(deg['name'], genes)]

# Split between pos and neg DEG
pos = deg[deg['logFC'] > 0]
neg = deg[deg['logFC'] < 0]

# Switch neg to pos values
neg = neg.assign(logFC=lambda x: np.abs(x.logFC))

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
    dc.run_ulm(slide, net, source='contrast', target='name', weight='logFC', use_raw=False)
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

sample_ids = [name for name in os.listdir('data/raw/visium/') if not name.startswith('.')]

for sample_id in sample_ids:

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

