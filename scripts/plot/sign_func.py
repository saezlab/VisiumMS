import pandas as pd
import numpy as np
import decoupler as dc
import os


# Pathway
ttl = pd.read_csv('data/prc/sign/pws/ttl.csv')
pos = pd.read_csv('data/prc/sign/pws/pos.csv')
neg = pd.read_csv('data/prc/sign/pws/neg.csv')

ttl = ttl.pivot(index='sample', columns='source', values='score')
pos = pos.pivot(index='sample', columns='source', values='score')
neg = neg.pivot(index='sample', columns='source', values='score')

pws_path = 'figures/sign/pws'
os.makedirs(pws_path, exist_ok=True)
ctypes = ttl.index.values
for ctype in ctypes:
    dc.plot_barplot(ttl, ctype, vertical=False,
                    save='{0}/{1}_{2}.pdf'.format(pws_path, ctype, 'ttl'))
for ctype in ctypes:
    dc.plot_barplot(pos, ctype, vertical=False,
                    save='{0}/{1}_{2}.pdf'.format(pws_path, ctype, 'pos'))
for ctype in ctypes:
    dc.plot_barplot(neg, ctype, vertical=False,
                    save='{0}/{1}_{2}.pdf'.format(pws_path, ctype, 'neg'))

# TFs
ttl = pd.read_csv('data/prc/sign/tfs/ttl.csv')
pos = pd.read_csv('data/prc/sign/tfs/pos.csv')
neg = pd.read_csv('data/prc/sign/tfs/neg.csv')

ttl = ttl.pivot(index='sample', columns='source', values='score')
pos = pos.pivot(index='sample', columns='source', values='score')
neg = neg.pivot(index='sample', columns='source', values='score')

tfs_path = 'figures/sign/tfs'
os.makedirs(tfs_path, exist_ok=True)
ctypes = ttl.index.values
for ctype in ctypes:
    dc.plot_barplot(ttl, ctype, vertical=True,
                    save='{0}/{1}_{2}.pdf'.format(tfs_path, ctype, 'ttl'))
for ctype in ctypes:
    dc.plot_barplot(pos, ctype, vertical=True,
                    save='{0}/{1}_{2}.pdf'.format(tfs_path, ctype, 'pos'))
for ctype in ctypes:
    dc.plot_barplot(neg, ctype, vertical=True,
                    save='{0}/{1}_{2}.pdf'.format(tfs_path, ctype, 'neg'))

# Ora
ttl = pd.read_csv('data/prc/sign/ora/ttl.csv')
pos = pd.read_csv('data/prc/sign/ora/pos.csv')
neg = pd.read_csv('data/prc/sign/ora/neg.csv')

ttl = ttl.pivot(index='sample', columns='source', values='pvals')
pos = pos.pivot(index='sample', columns='source', values='pvals')
neg = neg.pivot(index='sample', columns='source', values='pvals')

ttl.values[ttl.values == 0] = np.min(ttl.values[ttl.values != 0])
pos.values[pos.values == 0] = np.min(pos.values[pos.values != 0])
neg.values[neg.values == 0] = np.min(neg.values[neg.values != 0])
ttl = -np.log10(ttl)
pos = -np.log10(pos)
neg = -np.log10(neg)

ora_path = 'figures/sign/ora'
os.makedirs(ora_path, exist_ok=True)
ctypes = ttl.index.values
for ctype in ctypes:
    dc.plot_barplot(ttl, ctype, vertical=True,
                    save='{0}/{1}_{2}.pdf'.format(ora_path, ctype, 'ttl'))
for ctype in ctypes:
    dc.plot_barplot(pos, ctype, vertical=True,
                    save='{0}/{1}_{2}.pdf'.format(ora_path, ctype, 'pos'))
for ctype in ctypes:
    dc.plot_barplot(neg, ctype, vertical=True,
                    save='{0}/{1}_{2}.pdf'.format(ora_path, ctype, 'neg'))
