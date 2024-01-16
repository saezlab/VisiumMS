import pandas as pd
import numpy as np
import decoupler as dc
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import seaborn as sns
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--sn_deg_path', required=True)
parser.add_argument('-b','--ns_deg_path', required=True)
parser.add_argument('-c','--colors_dict', required=True)
parser.add_argument('-d','--plot_path', required=True)
args = vars(parser.parse_args())

sn_deg_path = args['sn_deg_path']
ns_deg_path = args['ns_deg_path']
colors_dict = args['colors_dict']
plot_path = args['plot_path']

# Get palette
palette = dict(item.split(':') for item in colors_dict.strip("'").split(';'))

# Single nuc
deg = pd.read_csv(sn_deg_path, index_col=0)
msk = deg['contrast'] == 'CAvsCI'
deg.loc[msk, 'log2FoldChange'] = - deg.loc[msk, 'log2FoldChange']
deg.loc[msk, 'stat'] = - deg.loc[msk, 'stat']
deg.loc[msk, 'contrast'] = 'CIvsCA'

ctypes = ['OL', 'AS', 'MG']
contrasts = ['CAvsCtrl', 'CIvsCtrl', 'CIvsCA']

def to_int(x):
    return int(np.round(x * 255))

def to_hex(row):
    r, g, b, a = row
    hex = '#%x%x%x' % (to_int(r), to_int(g), to_int(b))
    return hex

sn_figs = []
for contrast in contrasts:
    # Set color
    cond, ref = contrast.split('vs')
    down_color = palette[ref]
    up_color = palette[cond]
    rev_pal = {
        '#808080': '#808080',
        '#1f77b4': down_color,
        '#d62728': up_color,
    }
    for ctype in ctypes:
        fig, ax = plt.subplots(1, 1, figsize=(3, 3), dpi=150)
        dc.plot_volcano_df(
            data=deg[(deg['cell_type'] == ctype) & (deg['contrast'] == contrast)],
            x='log2FoldChange',
            y='pvalue',
            top=10,
            sign_thr=0.05,
            lFCs_thr=1,
            ax=ax,
            lFCs_limit=15,
        )
        ax.set_title('{t} {c}'.format(t=ctype, c=contrast))
        ax.get_children()[0].set_facecolor([rev_pal[to_hex(x)] for x in ax.get_children()[0].get_facecolor()])
        sn_figs.append(fig)

deg = deg[(deg['padj'] < 0.05) & (abs(deg['log2FoldChange']) > 1.)]
deg['condition'] = [c.split('vs')[0] if s > 0 else c.split('vs')[1] for s,
                    c in zip(deg['stat'], deg['contrast'])]
deg['sign'] = ['pos' if s > 0 else 'neg' for s in deg['stat']]

sn_fig, ax = plt.subplots(1, 1, figsize=(2, 2), dpi=150)
(
    deg
    .reset_index()
    .drop_duplicates(['index', 'cell_type', 'condition'], keep=False)
    .groupby(['cell_type', 'condition'])
    .count()[['stat']]
    .reset_index()
    .pivot(index='cell_type', columns='condition', values='stat')
    .fillna(0)
    .loc[['TC', 'SC', 'BC', 'NEU', 'EC', 'OPC', 'AS', 'MG', 'OL']]
    .plot.barh(stacked=True, ax=ax, width=0.85, color=palette)
)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
ax.set_xlabel('Number DEGs')
ax.set_ylabel('')
ax.grid(axis='x')
ax.set_axisbelow(True)

# Niches
# Read
deg = pd.read_csv(ns_deg_path, index_col=0)

ctypes = ['PPWM', 'LR', 'LC', 'VI']
contrasts = ['CIvsCA']

ns_figs = []
for contrast in contrasts:
    # Set color
    cond, ref = contrast.split('vs')
    down_color = palette[ref]
    up_color = palette[cond]
    rev_pal = {
        '#808080': '#808080',
        '#1f77b4': down_color,
        '#d62728': up_color,
    }
    for ctype in ctypes:
        fig, ax = plt.subplots(1, 1, figsize=(3, 3), dpi=150)
        dc.plot_volcano_df(
            data=deg[(deg['niche'] == ctype) & (deg['contrast'] == contrast)],
            x='log2FoldChange',
            y='pvalue',
            top=10,
            sign_thr=0.05,
            lFCs_thr=1,
            ax=ax,
            lFCs_limit=15,
        )
        ax.set_title('{t} {c}'.format(t=ctype, c=contrast))
        ax.get_children()[0].set_facecolor([rev_pal[to_hex(x)] for x in ax.get_children()[0].get_facecolor()])
        ns_figs.append(fig)
deg = deg[(deg['padj'] < 0.05) & (abs(deg['log2FoldChange']) > 1.)]
deg['condition'] = [c.split('vs')[0] if s > 0 else c.split('vs')[1] for s,
                    c in zip(deg['stat'], deg['contrast'])]
deg['sign'] = ['pos' if s > 0 else 'neg' for s in deg['stat']]

ns_fig, ax = plt.subplots(1, 1, figsize=(2, 2), dpi=150)
(
    deg
    .reset_index()
    .drop_duplicates(['index', 'niche', 'condition'], keep=False)
    .groupby(['niche', 'condition'])
    .count()[['stat']]
    .reset_index()
    .pivot(index='niche', columns='condition', values='stat')
    .fillna(0)
    .loc[['VI', 'LC', 'LR', 'PPWM']]
    .plot.barh(stacked=True, ax=ax, width=0.85, color=palette)
)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
ax.set_xlabel('Number DEGs')
ax.set_ylabel('')
ax.grid(axis='x')
ax.set_axisbelow(True)

# Save to pdf
pdf = matplotlib.backends.backend_pdf.PdfPages(plot_path)
for fg in [sn_fig] + sn_figs + [ns_fig] + ns_figs:
    pdf.savefig(fg, bbox_inches='tight')
pdf.close()
