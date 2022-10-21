import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

# Defina path
fig_path = 'data/cellbrowser/'

# Read data
adata = sc.read_h5ad('data/prc/sc/annotated.h5ad')
adata = adata.raw.to_adata()
raw = sc.read_h5ad('data/prc/sc/raw.h5ad')
inter = adata.var_names.intersection(adata.var_names)
adata = adata[:, inter]
raw = raw[:, inter]
adata.obs = adata.obs[['n_genes', 'n_genes_by_counts', 'total_counts',
                       'total_counts_mt', 'pct_counts_mt', 'doublet_score',
                       'diss_score', 'patient_id', 'sample_id', 'lesion_type',
                       'rin', 'sex', 'age', 'duration_disease', 'batch', 'leiden'
                      ]]
del adata.uns
del adata.obsp
adata.layers['counts'] = raw.X
mgc = sc.read_h5ad('data/prc/sc/microglia.h5ad').obs['leiden']
ast = sc.read_h5ad('data/prc/sc/astros.h5ad').obs['leiden']
adata.obs['MG_states'] = mgc
adata.obs['AS_states'] = ast
os.makedirs(fig_path, exist_ok=True)
adata.write(os.path.join(fig_path, 'snRNA_atlas.h5ad'))
del adata

def read_slide(sample_id):

    # Read rna-seq
    slide = sc.read_visium('data/raw/visium/{0}/outs/'.format(sample_id))
    slide.var_names_make_unique()
    sc.pp.filter_genes(slide, min_cells=3)
    sc.pp.filter_cells(slide, min_genes=200)

    # Store raw counts
    slide.layers['counts'] = slide.X.copy()

    # Normalize
    sc.pp.normalize_total(slide, target_sum=1e4)
    sc.pp.log1p(slide)

    # Read props
    props = pd.read_csv('data/prc/visium/{0}/cell_props.csv'.format(sample_id), index_col=0)
    inter = slide.obs.index.intersection(props.index)
    slide.obs = pd.concat([slide.obs, props.loc[inter]], axis=1)

    return slide


def get_ann_slide(sample_id, adata):
    msk = adata.obs['sample_id'].values == sample_id
    obs = adata[msk].obs.copy()
    progeny = adata.obsm['progeny'].loc[msk].copy()
    obs.index = [''.join([idx.split('-')[0], '-1']) for idx in obs.index]
    progeny.index = [''.join([idx.split('-')[0], '-1']) for idx in progeny.index]
    obs = obs[['sample_id', 'lesion_type', 'leiden']]
    slide = read_slide(sample_id)
    inter = obs.index.intersection(slide.obs.index)
    obs = obs.loc[inter]
    progeny = progeny.loc[inter]
    slide = slide[inter]
    slide.obs = pd.concat([slide.obs, obs, progeny], axis=1)
    return slide


# Visium slides
adata = sc.read_h5ad('data/prc/visium/props_niches.h5ad')
meta = pd.read_csv('data/metadata.csv').set_index('sample_id')

for sample_id in meta.index:
    slide = get_ann_slide(sample_id, adata)
    slide.write(os.path.join(fig_path, 'visium_{0}.h5ad'.format(sample_id)))
