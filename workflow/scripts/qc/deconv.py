import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
from scipy import stats
import os
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-m','--meta_path', required=True)
parser.add_argument('-a','--ann_path', required=True)
parser.add_argument('-s','--plot_sample_path', required=True)
parser.add_argument('-l','--plot_leiden_path', required=True)
args = vars(parser.parse_args())

meta_path = args['meta_path']
ann_path = args['ann_path']
plot_sample_path = args['plot_sample_path']
plot_leiden_path = args['plot_leiden_path']

meta = pd.read_csv(meta_path)
sample_ids = meta[['Sample id', 'Batch sn', 'Batch vs']].dropna()['Sample id'].values.astype('U')
obs = sc.read_h5ad(ann_path).obs
obs = obs[np.isin(obs['Sample id'], sample_ids)]

counts = obs.groupby(['Sample id', 'leiden']).count()[['total_counts']]
total = obs.groupby(['Sample id']).count()[['total_counts']]
df = (counts / total).dropna().rename(columns={'total_counts': 'sn'}).reset_index()

vs = []
for sample_id in sample_ids:
    props_path = 'data/prc/vs/{0}/props.csv'.format(sample_id)
    props = (
        pd.read_csv(props_path, index_col=0)
        .mean(0)
        .to_frame(name='vs')
        .reset_index(names='leiden')
    )
    props['Sample id'] = sample_id
    vs.append(props)
vs = pd.concat(vs)
df = pd.merge(df, vs)
df = df[(df['sn'] != 0.) & (df['vs'] != 0.)]
df['sn'] = np.log10(df['sn'].values)
df['vs'] = np.log10(df['vs'].values)
df = pd.merge(df, meta[['Sample id', 'Lesion type']], on='Sample id')

sn_min, sn_max = np.min(df['sn']), np.max(df['sn'])
vs_min, vs_max = np.min(df['vs']), np.max(df['vs'])
tt_min = np.min([vs_min, sn_min])

pdf = matplotlib.backends.backend_pdf.PdfPages(plot_sample_path)
for sample_id in sample_ids:
    tmp = df[df['Sample id'] == sample_id]
    sn = tmp['sn'].values
    vs = tmp['vs'].values
    lb = tmp['leiden'].values.astype('U')

    fig, ax = plt.subplots(1, 1, figsize=(2.5, 2.5), dpi=150)
    for i in range(sn.size):
        x, y, t = sn[i], vs[i], lb[i]
        ax.scatter(x=x, y=y, s=0)
        ax.text(x=x, y=y, s=t, weight='bold')
    ax.set_xlim(tt_min, 0.5)
    ax.set_ylim(tt_min, 0.5)
    ax.set_title(sample_id)
    ax.set_xlabel('SN (log10)')
    ax.set_ylabel('VS (log10)')
    ax.plot(np.unique(sn), np.poly1d(np.polyfit(sn, vs, 1))(np.unique(sn)), linestyle='--', color='gray')
    r, p = stats.pearsonr(sn, vs)
    ax.text(0.05, 0.9, 'R={:.2f}'.format(r), transform=ax.transAxes)
    ax.text(0.05, 0.8, 'p={:.2f}'.format(p), transform=ax.transAxes)
    ax.grid()
    ax.set_axisbelow(True)
    ax.margins(0.25, 0.25)
    pdf.savefig(fig)
pdf.close()

leidens = np.unique(df['leiden'].values.astype('U'))
pdf = matplotlib.backends.backend_pdf.PdfPages(plot_leiden_path)
for leiden in leidens:
    tmp = df[df['leiden'] == leiden]
    sn = tmp['sn'].values
    vs = tmp['vs'].values
    lb = tmp['Sample id'].values.astype('U')

    fig, ax = plt.subplots(1, 1, figsize=(2.5, 2.5), dpi=150)
    for i in range(sn.size):
        x, y, t = sn[i], vs[i], lb[i]
        ax.scatter(x=x, y=y, c='black')
        #ax.text(x=x, y=y, s=t, weight='bold')
    ax.set_xlim(tt_min, 0.5)
    ax.set_ylim(tt_min, 0.5)
    ax.set_title(leiden)
    ax.set_xlabel('SN (log10)')
    ax.set_ylabel('VS (log10)')
    ax.plot(np.unique(sn), np.poly1d(np.polyfit(sn, vs, 1))(np.unique(sn)), linestyle='--', color='gray')
    r, p = stats.pearsonr(sn, vs)
    ax.text(0.05, 0.9, 'R={:.2f}'.format(r), transform=ax.transAxes)
    ax.text(0.05, 0.8, 'p={:.2f}'.format(p), transform=ax.transAxes)
    ax.grid()
    ax.set_axisbelow(True)
    ax.margins(0.25, 0.25)
    pdf.savefig(fig)
pdf.close()

