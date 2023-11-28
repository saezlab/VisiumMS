#!/usr/bin/env python
# coding: utf-8
import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import os
import gc
import decoupler as dc
import regex as re

#Inputs
#Maybe unnecesary for snakemake
os.chdir("/Users/ricardoramirez/Dropbox/PostDoc/Research/MS")
input_sc = "./data/sn_annotated.h5ad"
output_meta = "./data/MOFAcell_meta.csv"
output_pb = "./data/MOFAcell_pb_data.csv"
output_coldat = "./data/MOFAcell_pb_coldata.csv"

#Reading object
sc_dat = sc.read_h5ad(filename = input_sc)

#Meta data for MOFAcell
meta_data_cols = ["Patient id", "Sample id", "Condition",
                  "Lesion type", "Age", "Sex", "RIN", "PMI (hrs)", 'Duration (years)',
                 'Cause of death','Batch sn']
meta_data = sc_dat.obs[meta_data_cols].drop_duplicates()
meta_data.to_csv(output_meta)

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
pb_dat.to_csv(output_pb)

pb_coldata = padata.obs.copy()
pb_coldata["colname"] = pb_coldata.index.values
pb_coldata = pb_coldata.merge(cell_type_numbers, on = ["leiden","Sample id"], how = "left")
pb_coldata.to_csv(output_coldat)