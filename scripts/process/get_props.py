import scanpy as sc
import numpy as np

import destvi_utils

from scvi.model import CondSCVI, DestVI
import destvi_utils

import os
import argparse


st_model = DestVI.load('data/prc/visium/MS377T/')
st_model.adata.obsm["proportions"] = st_model.get_proportions()
thrs = destvi_utils.automatic_proportion_threshold(st_model.adata, kind_threshold='secondary', output_file='data/prc/visium/MS377T/thr.html')

