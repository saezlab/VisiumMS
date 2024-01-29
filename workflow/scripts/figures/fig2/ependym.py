import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import seaborn as sns
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--slide_a', required=True)
parser.add_argument('-b','--slide_b', required=True)
parser.add_argument('-c','--niches_a', required=True)
parser.add_argument('-d','--niches_b', required=True)
parser.add_argument('-e','--plot_path', required=True)
args = vars(parser.parse_args())

slide_a = args['slide_a']
slide_b = args['slide_b']
niches_a = args['niches_a']
niches_b = args['niches_b']
plot_path = args['plot_path']


def read_slide(slide_path, niches_path):
    slide = sc.read_h5ad(slide_path)
    slide.obs['niches'] = pd.read_csv(niches_path, index_col=0)
    return slide

ms549h = read_slide(slide_a, niches_a)
ms549t = read_slide(slide_b, niches_b)

ms549h.obs['ependym'] = np.where(ms549h.obs['niches'] == 'Ependym', 'Ependym', 'Rest')
ms549t.obs['ependym'] = np.where(ms549t.obs['niches'] == 'Ependym', 'Ependym', 'Rest')

def plot_slide(slide, palette={'Ependym': '#e64040', 'Rest': '#73b360'}):
    fig, axes = plt.subplots(1, 3, figsize=(9, 3), dpi=150, tight_layout=True)
    axes = axes.ravel()
    sc.pl.spatial(slide, color=['ependym'], size=1.5, frameon=False, ax=axes[0], show=False, palette=palette, legend_loc=None)
    sc.pl.spatial(slide, color=['SPAG17'], size=1.5, frameon=False, ax=axes[1], show=False)
    sc.pl.violin(slide, keys='SPAG17', groupby='ependym', stripplot=False, rotation=90, inner='box', ax=axes[2], show=False, palette=palette)
    return fig

fig1 = plot_slide(ms549t)
fig2 = plot_slide(ms549h)

# Test for significance
sc.tl.rank_genes_groups(ms549t, groupby='ependym')
sc.tl.rank_genes_groups(ms549h, groupby='ependym')
print(sc.get.rank_genes_groups_df(ms549t, group='Ependym').set_index('names').loc['SPAG17'])
print(sc.get.rank_genes_groups_df(ms549h, group='Ependym').set_index('names').loc['SPAG17'])

# Save to pdf
pdf = matplotlib.backends.backend_pdf.PdfPages(plot_path)
for fig in [fig1, fig2]:
    pdf.savefig(fig, bbox_inches='tight')
pdf.close()
