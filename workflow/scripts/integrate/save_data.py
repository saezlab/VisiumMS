import pandas as pd
import numpy as np
import scanpy as sc
import os
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--sn_path', required=True)
parser.add_argument('-b','--out_path', required=True)
args = vars(parser.parse_args())

sn_path = args['sn_path']
out_path = args['out_path']

# Read sn atlas
adata = sc.read_h5ad(sn_path)

# Modify obs
cols = [
    'Patient id', 'Sample id', 'Condition', 'Lesion type', 'Age', 'Sex', 'RIN', 'PMI (hrs)',
    'Duration (years)', 'MS class', 'Cause of death', 'Batch sn', 'leiden'
]
obs = adata.obs.loc[:, cols].copy()
obs = obs.rename(columns={
    'Patient id': 'patient_id',
    'Sample id': 'sample_id',
    'Condition': 'condition',
    'Lesion type': 'lesion_type',
    'Age': 'age',
    'Sex': 'sex',
    'RIN': 'rin',
    'PMI (hrs)': 'pmi_hrs',
    'Duration (years)': 'duration_y',
    'MS class': 'ms_class',
    'Cause of death': 'cause_death',
    'Batch sn': 'batch_sn',
    'leiden': 'celltype'
})

# Add subtypes metadata
def read_ctype(ctype):
    df = pd.read_csv('data/prc/ctypes/{ctype}_ann.csv'.format(ctype=ctype), index_col=0).reset_index()
    df = df.rename(columns={'cell_states': 'leiden'}).set_index('index')['leiden']
    df.index.name = None
    return df
cstates = []
for ctype in np.unique(obs['celltype'].values.astype('U')):
    cstates.append(read_ctype(ctype))
cstates = pd.concat(cstates)
obs['subtype'] = cstates

# Update obs and trim adata
adata.obs = obs
del adata.uns
del adata.obsp

# Move counts to X
adata.X = adata.layers['counts'].copy()
del adata.layers['counts']

# Write
adata.write(os.path.join(out_path, 'sn_atlas.h5ad'))

def read_slide(sample_id):
    slide = sc.read_h5ad('data/prc/vs/{0}/adata.h5ad'.format(sample_id))
    slide.obs['niches'] = pd.read_csv('data/prc/vs/{0}/niches.csv'.format(sample_id), index_col=0)
    slide.obsm['ctype_abunds'] = pd.read_csv('data/prc/vs/{0}/abunds.csv'.format(sample_id), index_col=0)
    slide.obsm['ctype_props'] = pd.read_csv('data/prc/vs/{0}/props.csv'.format(sample_id), index_col=0)
    slide.obsm['hallmark_scores'] = pd.read_csv('data/prc/vs/{0}/hallmarks.csv'.format(sample_id), index_col=0)
    slide.obsm['reactome_scores'] = pd.read_csv('data/prc/vs/{0}/reactome.csv'.format(sample_id), index_col=0)
    slide.obsm['progeny_scores'] = pd.read_csv('data/prc/vs/{0}/progeny.csv'.format(sample_id), index_col=0)
    slide.obsm['ccc_scores'] = pd.read_csv('data/prc/vs/{0}/ctlr_scores.csv'.format(sample_id), index_col=0)
    del slide.uns['log1p']
    slide.X = slide.layers['counts'].copy()
    del slide.layers['counts']
    # Modify obs
    cols = [
        'array_row', 'array_col', 'Patient id', 'Sample id', 'Condition', 'Lesion type', 'Age', 'Sex', 'RIN',
        'PMI (hrs)', 'Duration (years)', 'MS class', 'Cause of death', 'Batch vs', 'areas', 'niches'
    ]
    obs = slide.obs.loc[:, cols].copy()
    obs = obs.rename(columns={
        'Patient id': 'patient_id',
        'Sample id': 'sample_id',
        'Condition': 'condition',
        'Lesion type': 'lesion_type',
        'Age': 'age',
        'Sex': 'sex',
        'RIN': 'rin',
        'PMI (hrs)': 'pmi_hrs',
        'Duration (years)': 'duration_y',
        'MS class': 'ms_class',
        'Cause of death': 'cause_death',
        'Batch vs': 'batch_vs',
    })
    slide.obs = obs
    # Remove emtpy areas
    msk = ~slide.obs['areas'].isna()
    print(np.any(~msk))
    slide = slide[msk, :].copy()
    return slide

# Read meta
meta = pd.read_csv('config/meta.csv')
st_samples = meta[~meta['Batch vs'].isnull()]['Sample id'].values.astype('U')

# Write each visium sample
for sample_id in st_samples:
    slide = read_slide(sample_id)
    slide.write(os.path.join(out_path, 'visium_{0}.h5ad'.format(sample_id)))
