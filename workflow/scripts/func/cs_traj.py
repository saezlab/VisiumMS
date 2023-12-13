import pandas as pd
import numpy as np
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--ann_path', required=True)
parser.add_argument('-d','--deg_path', required=True)
parser.add_argument('-k','--ks_path', required=True)
parser.add_argument('-w','--wt_path', required=True)
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

ann_path = args['ann_path']
deg_path = args['deg_path']
ks_path = args['ks_path']
wt_path = args['wt_path']
out_path = args['out_path']


# Read deg
deg = pd.read_csv(deg_path, index_col=0)
ks_df = pd.read_csv(ks_path)
wt_df = pd.read_csv(wt_path)

# Flip direction contrast
msk = deg['contrast'] == 'CAvsCI'
deg.loc[msk, 'stat'] = - deg.loc[msk, 'stat'].values
deg.loc[msk, 'contrast'] = 'CIvsCA'

# Subset
order = ['CAvsCtrl', 'CIvsCA']
deg = deg[np.isin(deg['contrast'], order)]

def get_cdeg_cats(deg, ctype):
    cde = deg[deg['cell_type'] == ctype]
    cstat = cde.pivot(columns='contrast', values='stat')
    cpval = cde.pivot(columns='contrast', values='padj')
    msk = np.sum(cpval.values < 0.05, axis=1).astype(bool).ravel()
    cstat, cpval = cstat.loc[msk], cpval.loc[msk]
    ccats = pd.DataFrame(index=cstat.index)
    ccats['cat'] = ''
    for i in range(ccats.shape[0]):
        s_a, s_b = cstat.iloc[i].values
        p_a, p_b = cpval.iloc[i].values
        a, b = 0, 0
        if p_a < 0.05:
            a = int(np.sign(s_a))
        if p_b < 0.05:
            b = int(np.sign(s_b))
        ccats.iloc[i, 0] = '{a}_{b}'.format(a=a, b=b)
    return ccats

def get_cdeg(deg, wt_df, ctype):
    cdeg = get_cdeg_cats(deg, ctype=ctype).reset_index(names='names')
    cdeg['ctype'] = ctype
    cs_names = wt_df[(wt_df['type'] == 'cs') & (wt_df['padj'] < 0.05) & (wt_df['ctype'].str.startswith('{ctype}_'.format(ctype=ctype)))]['ctype'].unique()
    mrk = pd.read_csv('data/prc/ctypes/{ctype}_deg.csv'.format(ctype=ctype))
    mrk = mrk[np.isin(mrk['group'].astype('U'), cs_names)]
    mrk = mrk.groupby('names')['group'].apply(lambda x: list(x)).reset_index()
    cdeg = pd.merge(cdeg, mrk).rename(columns={'group': 'cell_states'})
    return cdeg

df = []
for ctype in ['OL', 'AS', 'MG']:
    tmp = get_cdeg(deg, wt_df, ctype=ctype)
    df.append(tmp)
df = pd.concat(df)

# Write
df.to_csv(out_path, index=False)
