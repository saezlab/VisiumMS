import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--sn_deg_path', required=True)
parser.add_argument('-b','--colors_dict', required=True)
parser.add_argument('-c','--plot_path', required=True)
args = vars(parser.parse_args())

sn_deg_path = args['sn_deg_path']
colors_dict = args['colors_dict']
plot_path = args['plot_path']

# Get palette
palette = dict(item.split(':') for item in colors_dict.strip("'").split(';'))

# Read
deg = pd.read_csv('data/prc/sn_deg.csv', index_col=0)
deg = deg[(deg['padj'] < 0.05) & (abs(deg['log2FoldChange']) > 1.)]
deg['condition'] = [c.split('vs')[0] if s > 0 else c.split('vs')[1] for s,
                    c in zip(deg['stat'], deg['contrast'])]
deg['sign'] = ['pos' if s > 0 else 'neg' for s in deg['stat']]
tab = deg.reset_index().groupby('condition')['index'].apply(lambda x: set(x))
labels = tab.index
sets = tab.values

# Plot
fig, ax = plt.subplots(1, 1, figsize=(4, 4), dpi=150)
v = venn3(sets, labels, set_colors=[palette[l] for l in labels], alpha=0.75, ax=ax)
c = venn3_circles(sets, linewidth=1, ax=ax)

# Save
fig.savefig(plot_path, bbox_inches='tight')
