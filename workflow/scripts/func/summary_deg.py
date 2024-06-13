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

gdict = {
    'AS': np.array([
    'PDK1','ACSF2','TANGO2','AQP1','ARL17B','SMAD6','HMGB1','HLA-F-AS1', 'LRP1', 'PRKG2','OSMR',
    'SERPINA3','NFE2L2', 'APP', 'C3', 'DST','TNFRSF1A', 'EEA1','DOCK11','RNF7','CCND2',
    'CNDP2','HSPBP1', 'FOXJ1','ANXA1','SLITRK2','ITGA2','SCRG1','CLU','ITGB1','CFAP97',
    'EBLN3P', 'COL6A1','CTSD','IFITM3','TAGLN','AKT3','CPNE3','KCNK2','PIK3IP1',
    'IFITM2', 'ARL5B', 'NNT', 'CUL4B', 'ALDH1L1']),
    'MG': np.array([
    'PRAM1','ARHGAP12','TMEM119','HAMP','P2RY12',
    'HMOX1','SIRPA','SPP1','TRAF3','PARP9','FLT1','PARVG','MS4A6A','CD14','PPARG','COPB2','C1QB',
    'APOC1','FTL','C1QC','GPNMB','APP','AHI1','FPR3','CD163','CTSS','IQGAP1',
    'PLXNC1','APOE','CD68','ITGA4','ANXA2','LGALS1','TGFBI','MS4A7','ITGB1',
    'PIK3CG','PDCD10','AGPS','LRRK2','CLEC12A','CLEC7A','FRMD4A','CX3CR1', 'FOXP2']),
    'OL': np.array([
    'CAMK2A','ERBB2','CD226','PGM1','SOX13','ADAMTS4','PSEN1','DOCK3','NDE1','EPHB1',
    'CDK18','AXIN1','SRA1','HSPB1','IRF1','EIF5','NFKB2','CALM1','THAP1','HSP90B1',
    'CAMK2D','USP1','NGFR','ATF4','CD274','SLC22A17','EIF2S1','DCC','BRCA2',
    'ITGB1','SOX4','LGALS3','OSMR','ANXA2','MYRIP','MPZ','GMFB','TGFBR2','VWA8','CD2AP',
    'ARNTL','ELOVL6','SVEP1','RORA','SCD']),
}

def to_int(x):
    return int(np.round(x * 255))

def to_hex(row):
    r, g, b, a = row
    hex = '#%x%x%x' % (to_int(r), to_int(g), to_int(b))
    return hex

def plot_volcano_df(data, x, y, top=5, sign_thr=0.05, lFCs_thr=0.5, sign_limit=None, lFCs_limit=None, color_pos='#D62728',
                    color_neg='#1F77B4', color_null='gray', figsize=(7, 5), dpi=100, ax=None, return_fig=False, save=None):

    # Load plotting packages
    import adjustText as at

    # Transform sign_thr
    sign_thr = -np.log10(sign_thr)

    # Extract df
    df = data.copy()
    df['logFCs'] = df[x]
    df['pvals'] = -np.log10(df[y])

    # Define color by up or down regulation and significance
    df['weight'] = color_null
    up_msk = (df['logFCs'] >= lFCs_thr) & (df['pvals'] >= sign_thr)
    dw_msk = (df['logFCs'] <= -lFCs_thr) & (df['pvals'] >= sign_thr)
    df.loc[up_msk, 'weight'] = color_pos
    df.loc[dw_msk, 'weight'] = color_neg

    # Plot
    fig = None
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=dpi)
    df.plot.scatter(x='logFCs', y='pvals', c='weight', sharex=False, ax=ax)
    ax.set_axisbelow(True)

    # Draw sign lines
    ax.axhline(y=sign_thr, linestyle='--', color="black")
    ax.axvline(x=lFCs_thr, linestyle='--', color="black")
    ax.axvline(x=-lFCs_thr, linestyle='--', color="black")

    # Plot top sign features
    signs = df[up_msk | dw_msk].sort_values('pvals', ascending=False)
    signs = signs.loc[top[np.isin(top, signs.index.values.astype('U'))]]

    # Add labels
    ax.set_ylabel('-log10(pvals)')
    texts = []
    for x, y, s in zip(signs['logFCs'], signs['pvals'], signs.index):
        texts.append(ax.text(x, y, s))
    if len(texts) > 0:
        at.adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black'), ax=ax)

    if return_fig:
        return fig

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
        top = gdict[ctype]
        fig, ax = plt.subplots(1, 1, figsize=(3, 3), dpi=150)
        plot_volcano_df(
            data=deg[(deg['cell_type'] == ctype) & (deg['contrast'] == contrast)],
            x='log2FoldChange',
            y='pvalue',
            top=top,
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
top = []
for k in gdict:
    top.extend(list(gdict[k]))
top = np.unique(top)
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
        plot_volcano_df(
            data=deg[(deg['niche'] == ctype) & (deg['contrast'] == contrast)],
            x='log2FoldChange',
            y='pvalue',
            top=top,
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
