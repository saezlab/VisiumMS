import scanpy as sc
import pandas as pd
import numpy as np
import os
from composition_stats import closure
from composition_stats import ilr
import scanpy.external as sce
import decoupler as dc


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
    slide.obsm['props'] = props.loc[inter]

    return slide

# Open metadata
meta = pd.read_csv(os.path.join('data','metadata.csv')).set_index('sample_id')

# Read all slides
slides = []
for sample_id in meta.index.values:
    slide = read_slide(sample_id)
    slide.obs['sample_id'] = sample_id
    slide.obs['lesion_type'] = meta.loc[sample_id, 'lesion_type']
    slides.append(slide)
slides = slides[0].concatenate(slides[1:], join='outer')

# Generate ILR space
slides.obsm['props'].loc[:, :] = closure(slides.obsm['props'].values)
slides.obsm['ilr'] = ilr(slides.obsm['props'].values)

# Integrate
sce.pp.harmony_integrate(slides, 'batch', basis='ilr', max_iter_harmony=30)

# Generate UMAP features
sc.pp.neighbors(slides, use_rep='X_pca_harmony')
sc.tl.umap(slides)

# Run leiden clustering algorithm
sc.tl.leiden(slides, resolution=0.65)

# progeny
progeny = dc.get_progeny(organism='human', top=300)
dc.run_mlm(slides, progeny, use_raw=False)
slides.obsm['progeny'] = slides.obsm['mlm_estimate'].copy()

# Save
slides.write('data/prc/visium/props_niches.h5ad')
