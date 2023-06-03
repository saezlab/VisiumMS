import pandas as pd
import numpy as np
import decoupler as dc
import os

import matplotlib.pyplot as plt

# Read and concat
pos = pd.read_csv('data/prc/sign/deg/pos.csv')
neg = pd.read_csv('data/prc/sign/deg/neg.csv')
neg['logFCs'] = -neg['logFCs'].values 
deg = pd.concat([pos, neg])

# Expand
logFCs = deg.pivot(index='contrast', columns='name', values='logFCs').fillna(0.)
pvals = deg.pivot(index='contrast', columns='name', values='pvals').fillna(1.)

# Find unique cells
cell_types = np.unique([np.unique(pos['contrast'].values), np.unique(neg['contrast'])])
n_cells = cell_types.size

# Plot
res_path = 'figures/sign/deg'
os.makedirs(res_path, exist_ok=True)
for i in range(n_cells):
    ctype = cell_types[i]
    dc.plot_volcano(logFCs, pvals, ctype, lFCs_limit=10, top=10,
                    save='{0}/{1}.pdf'.format(res_path, ctype))
