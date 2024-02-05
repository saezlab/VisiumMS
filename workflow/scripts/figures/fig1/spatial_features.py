import pandas as pd
import numpy as np
import scanpy as sc
import decoupler as dc
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-d','--colors_dict', required=True)
parser.add_argument('-p','--plot_path', required=True)
args = vars(parser.parse_args())

colors_dict = args['colors_dict']
plot_path = args['plot_path']


def read_slide(sample_id):
    slide = sc.read_h5ad('data/prc/vs/{0}/adata.h5ad'.format(sample_id))
    slide.obs['niches'] = pd.read_csv('data/prc/vs/{0}/niches.csv'.format(sample_id), index_col=0)
    slide.obsm['reactome'] = pd.read_csv('data/prc/vs/{0}/reactome.csv'.format(sample_id), index_col=0)
    slide.obsm['props'] = pd.read_csv('data/prc/vs/{0}/props.csv'.format(sample_id), index_col=0)
    return slide

# Read slides
slide_ct = read_slide('CO40')
slide_ca = read_slide('MS377T')
slide_ci = read_slide('MS497I')

# Get palette
colors_dict = dict(item.split('_') for item in colors_dict.strip("'").split(';'))

# Plot areas
fig1, axes = plt.subplots(3, 3, dpi=300, figsize=(5, 6))

ax = axes[0, 0]
sc.pl.spatial(slide_ct, color=None, size=1.5, frameon=False, legend_loc=None,
              title='', ax=ax, show=False, palette=colors_dict)

ax = axes[1, 0]
sc.pl.spatial(slide_ca, color=None, size=1.5, frameon=False, legend_loc=None,
              title='', ax=ax, show=False, palette=colors_dict)

ax = axes[2, 0]
sc.pl.spatial(slide_ci, color=None, size=1.5, frameon=False, legend_loc=None,
              title='', ax=ax, show=False, palette=colors_dict)

ax = axes[0, 1]
sc.pl.spatial(slide_ct, color='areas', size=1.5, frameon=False, legend_loc=None,
              title='', ax=ax, show=False, palette=colors_dict)

ax = axes[1, 1]
sc.pl.spatial(slide_ca, color='areas', size=1.5, frameon=False, legend_loc=None,
              title='', ax=ax, show=False, palette=colors_dict)

ax = axes[2, 1]
sc.pl.spatial(slide_ci, color='areas', size=1.5, frameon=False, legend_loc=None,
              title='', ax=ax, show=False, palette=colors_dict)

ax = axes[0, 2]
sc.pl.spatial(slide_ct, color='niches', size=1.5, frameon=False, legend_loc=None,
              title='', ax=ax, show=False, palette=colors_dict)

ax = axes[1, 2]
sc.pl.spatial(slide_ca, color='niches', size=1.5, frameon=False, legend_loc=None,
              title='', ax=ax, show=False, palette=colors_dict)

ax = axes[2, 2]
sc.pl.spatial(slide_ci, color='niches', size=1.5, frameon=False, legend_loc=None,
              title='', ax=ax, show=False, palette=colors_dict)

fig1.subplots_adjust(hspace=0, wspace=0.05)

# Plot props
props_ct = dc.get_acts(slide_ct, 'props')
props_ca = dc.get_acts(slide_ca, 'props')
props_ci = dc.get_acts(slide_ci, 'props')

fig2, axes = plt.subplots(3, 3, dpi=300, figsize=(5, 6))

ax = axes[0, 0]
sc.pl.spatial(props_ct, color='OL', size=1.5, frameon=False, legend_loc=None, cmap='magma',
              title='', ax=ax, show=False, vmin=0, vmax=1, colorbar_loc=None)

ax = axes[1, 0]
sc.pl.spatial(props_ca, color='OL', size=1.5, frameon=False, legend_loc=None, cmap='magma',
              title='', ax=ax, show=False, vmin=0, vmax=1, colorbar_loc=None)

ax = axes[2, 0]
sc.pl.spatial(props_ci, color='OL', size=1.5, frameon=False, legend_loc=None, cmap='magma',
              title='', ax=ax, show=False, vmin=0, vmax=1, colorbar_loc=None)

ax = axes[0, 1]
sc.pl.spatial(props_ct, color='MG', size=1.5, frameon=False, legend_loc=None, cmap='magma',
              title='', ax=ax, show=False, vmin=0, vmax=1, colorbar_loc=None)

ax = axes[1, 1]
sc.pl.spatial(props_ca, color='MG', size=1.5, frameon=False, legend_loc=None, cmap='magma',
              title='', ax=ax, show=False, vmin=0, vmax=1, colorbar_loc=None)

ax = axes[2, 1]
sc.pl.spatial(props_ci, color='MG', size=1.5, frameon=False, legend_loc=None, cmap='magma',
              title='', ax=ax, show=False, vmin=0, vmax=1, colorbar_loc=None)

ax = axes[0, 2]
sc.pl.spatial(props_ct, color='AS', size=1.5, frameon=False, legend_loc=None, cmap='magma',
              title='', ax=ax, show=False, vmin=0, vmax=1, colorbar_loc=None)

ax = axes[1, 2]
sc.pl.spatial(props_ca, color='AS', size=1.5, frameon=False, legend_loc=None, cmap='magma',
              title='', ax=ax, show=False, vmin=0, vmax=1, colorbar_loc=None)

ax = axes[2, 2]
sc.pl.spatial(props_ci, color='AS', size=1.5, frameon=False, legend_loc=None, cmap='magma',
              title='', ax=ax, show=False, vmin=0, vmax=1, colorbar_loc=None)

fig2.subplots_adjust(hspace=0, wspace=0.05)

# Plot pathways
rc_ct = dc.get_acts(slide_ct, 'reactome').copy()
rc_ca = dc.get_acts(slide_ca, 'reactome').copy()
rc_ci = dc.get_acts(slide_ci, 'reactome').copy()

inter = np.intersect1d(rc_ct.var_names, rc_ca.var_names)
inter = np.intersect1d(inter, rc_ci.var_names)
rc_ct, rc_ca, rc_ci = rc_ct[:, inter].copy(), rc_ca[:, inter].copy(), rc_ci[:, inter].copy()

# Scale
X = np.concatenate([rc_ct.X, rc_ca.X, rc_ci.X], axis=0)
means = np.mean(X, axis=0)
stds = np.std(X, axis=0, ddof=1)
rc_ct.X = (rc_ct.X - means) / stds
rc_ca.X = (rc_ca.X - means) / stds
rc_ci.X = (rc_ci.X - means) / stds

def plot_element_rc(adata, gset, ax):
    sc.pl.spatial(adata, color=gset, size=1.5, frameon=False, legend_loc=None, cmap='RdBu_r',
              title='', ax=ax, show=False, vmin=-2, vmax=2, colorbar_loc=None)

def plot_column_rc(adata1, adata2, adata3, axes, col_i, gset):
    plot_element_rc(adata1, gset, ax=axes[0, col_i])
    plot_element_rc(adata2, gset, ax=axes[1, col_i])
    plot_element_rc(adata3, gset, ax=axes[2, col_i])


fig3, axes = plt.subplots(3, 3, dpi=150, figsize=(5, 6))
plot_column_rc(rc_ct, rc_ca, rc_ci, axes, 0, 'EGR2 AND SOX10 MEDIATED INITIATION OF SCHWANN CELL MYELINATION')
plot_column_rc(rc_ct, rc_ca, rc_ci, axes, 1, 'INTERFERON GAMMA SIGNALING')
plot_column_rc(rc_ct, rc_ca, rc_ci, axes, 2, 'ASSEMBLY OF COLLAGEN FIBRILS AND OTHER MULTIMERIC STRUCTURES') #'MUSCLE CONTRACTION')
fig3.subplots_adjust(hspace=0, wspace=0.05)

# Save to pdf
pdf = matplotlib.backends.backend_pdf.PdfPages(plot_path)
for fig in [fig1, fig2, fig3]:
    pdf.savefig(fig)
pdf.close()
