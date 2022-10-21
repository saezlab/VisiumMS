import pandas as pd
import numpy as np
import decoupler as dc
import os

# Read and concat
pos = pd.read_csv('data/prc/sign/deg/pos.csv')
neg = pd.read_csv('data/prc/sign/deg/neg.csv')
neg['logFCs'] = -neg['logFCs'].values 
deg = pd.concat([pos, neg])

# Expand
logFCs = deg.pivot(index='contrast', columns='name', values='logFCs').fillna(0.)
pvals = deg.pivot(index='contrast', columns='name', values='pvals').fillna(1.)
p_logFCs = pos.pivot(index='contrast', columns='name', values='logFCs').fillna(0.)
p_pvals = pos.pivot(index='contrast', columns='name', values='pvals').fillna(1.)
n_logFCs = -neg.pivot(index='contrast', columns='name', values='logFCs').fillna(0.)
n_pvals = neg.pivot(index='contrast', columns='name', values='pvals').fillna(1.)

# Retrieve functonal networks
progeny = dc.get_progeny(top=300)
dorothea = dc.get_dorothea()
msigdb = dc.get_resource('MSigDB')
msigdb = msigdb[~msigdb.duplicated(['geneset', 'genesymbol'])]
reactome = msigdb[msigdb['collection'] == 'reactome_pathways']

def infer_act(logFCs, pvals, net):

    # Estimate stat
    stat = -np.log10(pvals) * logFCs

    # Infer act with consensus
    acts = dc.run_consensus(mat=stat, net=net)

    # Format
    acts = (
        dc.melt(acts)
        .sort_values(['sample', 'pvals'])
        .drop('method', axis=1)
    )

    return acts

def enrich(deg, net):

    # Infer enrichment with ora
    res = dc.get_ora_df(deg, net, groupby='contrast', features='name',
                        source='geneset', target='genesymbol')

    # Format
    res = (
        dc.melt(res)
        .sort_values(['sample', 'pvals'])
    )

    return res


# Pathways
t_pws = infer_act(logFCs, pvals, progeny)
p_pws = infer_act(p_logFCs, p_pvals, progeny)
n_pws = infer_act(n_logFCs, n_pvals, progeny)

pws_path = 'data/prc/sign/pws'
os.makedirs(pws_path, exist_ok=True)
t_pws.to_csv('{0}/{1}.csv'.format(pws_path, 'ttl'), index=False)
p_pws.to_csv('{0}/{1}.csv'.format(pws_path, 'pos'), index=False)
n_pws.to_csv('{0}/{1}.csv'.format(pws_path, 'neg'), index=False)

# TFs
t_tfs = infer_act(logFCs, pvals, dorothea)
p_tfs = infer_act(p_logFCs, p_pvals, dorothea)
n_tfs = infer_act(n_logFCs, n_pvals, dorothea)

tfs_path = 'data/prc/sign/tfs'
os.makedirs(tfs_path, exist_ok=True)
t_tfs.to_csv('{0}/{1}.csv'.format(tfs_path, 'ttl'), index=False)
p_tfs.to_csv('{0}/{1}.csv'.format(tfs_path, 'pos'), index=False)
n_tfs.to_csv('{0}/{1}.csv'.format(tfs_path, 'neg'), index=False)

# Enrichment
t_ora = enrich(deg, reactome)
p_ora = enrich(pos, reactome)
n_ora = enrich(neg, reactome)

ora_path = 'data/prc/sign/ora'
os.makedirs(ora_path, exist_ok=True)
t_ora.to_csv('{0}/{1}.csv'.format(ora_path, 'ttl'), index=False)
p_ora.to_csv('{0}/{1}.csv'.format(ora_path, 'pos'), index=False)
n_ora.to_csv('{0}/{1}.csv'.format(ora_path, 'neg'), index=False)
