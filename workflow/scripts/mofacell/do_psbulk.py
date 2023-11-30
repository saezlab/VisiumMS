import sys
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import os
import gc
import decoupler as dc
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--input_sc', required=True)
parser.add_argument('-m','--meta_path', required=True)
parser.add_argument('-a','--sn_output_pb', required=True)
parser.add_argument('-b','--sn_output_coldat', required=True)
parser.add_argument('-c','--vs_output_pb', required=True)
parser.add_argument('-d','--vs_output_coldat', required=True)
args = vars(parser.parse_args())

input_sc = args['input_sc']
meta_path = args['meta_path']
sn_output_pb = args['sn_output_pb']
sn_output_coldat = args['sn_output_coldat']
vs_output_pb = args['vs_output_pb']
vs_output_coldat = args['vs_output_coldat']

# Single nuc
sc_dat = sc.read_h5ad(filename = input_sc)

cell_type_numbers = sc_dat.obs.groupby(["Sample id", "leiden"])["leiden"].count()
cell_type_numbers = cell_type_numbers.to_frame().rename({"leiden":"counts"}, axis=1).reset_index()

padata = dc.get_pseudobulk(sc_dat,
                           sample_col='Sample id', 
                           groups_col='leiden', 
                           layer='counts', 
                           min_prop=0, 
                           min_smpls=0,
                           mode='sum')

pb_dat = pd.DataFrame(padata.X)
pb_dat.columns = padata.var.index.values
pb_dat.index = padata.obs.index.values
pb_dat.to_csv(sn_output_pb)

pb_coldata = padata.obs.copy()
pb_coldata["colname"] = pb_coldata.index.values
pb_coldata = pb_coldata.merge(cell_type_numbers, on = ["leiden","Sample id"], how = "left")
pb_coldata.to_csv(sn_output_coldat)

# Visium
meta = pd.read_csv(meta_path)
vs_samples = meta[~meta['Batch vs'].isnull()]['Sample id'].values.astype('U')
adata = []
for sample_id in vs_samples:
    slide = sc.read_h5ad('data/prc/vs/{0}/adata.h5ad'.format(sample_id))
    slide.obs['niches'] = pd.read_csv('data/prc/vs/{0}/niches.csv'.format(sample_id), index_col=0)
    slide.obs_names = [sample_id + '_' + i for i in slide.obs_names]
    slide = slide[~slide.obs['areas'].isnull()].copy()
    adata.append(slide)
adata = ad.concat(adata, join='outer')

padata = dc.get_pseudobulk(adata,
                           sample_col='Sample id',
                           groups_col='niches',
                           layer='counts',
                           min_prop=0,
                           min_smpls=0,
                           mode='sum')

pb_dat = pd.DataFrame(padata.X)
pb_dat.columns = padata.var.index.values
pb_dat.index = padata.obs.index.values
pb_dat.to_csv(vs_output_pb)

pb_coldata = padata.obs.copy()
pb_coldata["colname"] = pb_coldata.index.values
pb_coldata.to_csv(vs_output_coldat)
