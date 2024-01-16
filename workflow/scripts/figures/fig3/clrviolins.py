import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import seaborn as sns
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--dfs_path', required=True)
parser.add_argument('-b','--wt_path', required=True)
parser.add_argument('-c','--plot_path', required=True)
args = vars(parser.parse_args())

dfs_path = args['dfs_path']
wt_path = args['wt_path']
plot_path = args['plot_path']

# Read
dfs = pd.read_csv(dfs_path)
wt = pd.read_csv(wt_path)

def extract_df(df, mode, type, leiden):
    return df[(df['mode'] == mode) & (df['type'] == type) & (df['leiden'] == leiden)]

def plot_violin(df, mode, type, leiden, ax, palette):
    data = extract_df(df, mode=mode, type=type, leiden=leiden)
    sns.violinplot(data=data, x='Lesion type', y='value', ax=ax, palette=palette)
    ax.set_xlabel('')
    ax.set_ylabel(leiden)
    ax.set_title(type)

palette = {
    'Ctrl': '#3aeba2',
    'CA': '#8be9f1',
    'CI': '#ffa4e8',
}

# Make unique
pairs = (
    wt[wt['padj'] < 0.10]
    .sort_values('type', ascending=False)
    .drop_duplicates(['ctype', 'type'])
    [['ctype', 'type']]
)

# Plot
types = pairs['type'].unique()
figs = []
for type in types:
    ctypes = pairs[pairs['type'] == type]['ctype'].values.ravel()
    n_ct = len(ctypes)
    fig, axes = plt.subplots(1, n_ct, figsize=(2 * n_ct, 2), tight_layout=True, dpi=150)
    for i in range(n_ct):
        ctype = ctypes[i]
        ax = axes[i]
        plot_violin(dfs, mode='clr', type=type, leiden=ctype, ax=ax, palette=palette)
    figs.append(fig)

# Save to pdf
pdf = matplotlib.backends.backend_pdf.PdfPages(plot_path)
for fig in figs:
    pdf.savefig(fig, bbox_inches='tight')
pdf.close()
