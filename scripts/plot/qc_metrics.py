
# python scripts/plot/qc_metrics.py --output cellbender
# python scripts/plot/qc_metrics.py --output cellranger

import numpy as np
import pandas as pd
import pickle

import os
import argparse
from pathlib import Path
from plotting import *

"""
Script to plot different QC metrics after filtering the data.
"""

# add command line flag arguments to specify either "cellbender" or "cellranger" output
parser = argparse.ArgumentParser()
parser.add_argument("--output", type=str, required=True)
args = parser.parse_args()

# set up relative paths within the project
current_folder = Path(__file__).parent
if args.output == "cellbender":
    input_dir = current_folder / ".." / ".." / "data" / "prc" / "sc" / "cellbender_qc"
    output_dir = current_folder / ".." / ".." / "out" / "cellbender_qc"
elif args.output == "cellranger":
    input_dir = current_folder / ".." / ".." / "data" / "prc" / "sc" / "cellranger_qc"
    output_dir = current_folder / ".." / ".." / "out" / "cellranger_qc"
else:
    raise ValueError("output must be either 'cellbender' or 'cellranger'")
output_dir.mkdir(parents=True, exist_ok=True)

samples = [s.split(".")[0] for s in os.listdir(input_dir) if s.endswith('.pkl')]

# Initialize sumary variables
total_df = pd.DataFrame(columns=['n_genes_by_counts','total_counts',
                                   'pct_counts_mt','doublet_score'])
summary_df = []
total_n_rem = np.zeros(5)

# Run QC plots for each sample and store summary
for sample in samples:
    print(sample)
    plot_data = pickle.load(open(input_dir / (sample + ".pkl"), "rb"))

    # Filter params
    mt_thr = plot_data['mt_thr']
    gene_qnt = plot_data['gene_qnt']
    doublet_thr = plot_data['doublet_thr']
    diss_qnt = plot_data['diss_qnt']
    df = plot_data['df']
    gene_thr = np.quantile(df.n_genes_by_counts, gene_qnt)
    diss_thr = np.quantile(df.diss_score, diss_qnt)

    # Create QC plots
    fig, axes = plt.subplots(2,3, figsize=(11,8.5), dpi=150)
    fig.suptitle('QC metrics {0}'.format(sample), fontsize=11)
    axes = axes.flatten()

    # Plot MT
    plot_mt_vs_counts(df, axes[0], mt_thr=mt_thr)

    # Plot ngenes
    plot_ngenes_vs_counts(df, axes[1], gene_thr=gene_thr)

    # Plot lost cells
    labels = ['MT', 'Gene', 'Doublet', 'Diss', 'Total']
    msks = np.array([
        df.pct_counts_mt > mt_thr,
        df.n_genes_by_counts > gene_thr,
        df.doublet_score > doublet_thr,
        df.diss_score > diss_thr
    ])
    msks = np.vstack((msks, [np.sum(msks, axis=0) > 0]))
    n_rem = np.sum(msks, axis=1)
    total_n_rem += n_rem
    plot_ncell_diff(df, axes[2], labels=labels, n_rem=n_rem)
    
    # Plot doublet score
    plot_doublet_scores(df, axes[3], doublet_thr=doublet_thr)

    # Plot diss score
    plot_diss_scores(df, axes[4], diss_thr=diss_thr)
    
    # Filter
    msk = (df.n_genes_by_counts < gene_thr) & \
      (df.pct_counts_mt < mt_thr) &  \
      (df.doublet_score < doublet_thr) & \
      (df.diss_score < diss_thr)
    df = df[msk]

    # Plot diff cells
    pd.DataFrame([[sample, len(msk), sum(msk)]], 
                 columns=["Sample ID", "Before", "After"]).plot.bar(x=0, width=0.5, ax=axes[5])
    axes[5].set_title('Lost cells', fontsize=11)
    
    # Adjust plots
    fig.tight_layout()
    fig.subplots_adjust(top=0.88)
    fig.set_facecolor('white')

    # Write to png
    fig.savefig(output_dir / ("qc_"+sample+".png"))

    # Append
    total_df = pd.concat([total_df, df], ignore_index=True)
    summary_df.append([sample, len(msk), sum(msk)])
        
# Summary plots
fig, axes = plt.subplots(2,3, figsize=(11,8.5), dpi=150)
fig.suptitle('Merged QC metrics', fontsize=11)
axes = axes.flatten()

plot_mt_vs_counts(total_df, axes[0], mt_thr=np.nan)

plot_ngenes_vs_counts(total_df, axes[1], gene_thr=np.nan)

plot_ncell_diff(df, axes[2], labels=labels, n_rem=total_n_rem)

plot_doublet_scores(total_df, axes[3], doublet_thr=np.nan)

plot_diss_scores(total_df, axes[4], diss_thr=np.nan)

pd.DataFrame(summary_df, columns=["Sample ID", "Before", "After"]).plot.bar(x=0, ax=axes[5])
axes[5].set_title('Lost cells', fontsize=11)

fig.tight_layout()
fig.subplots_adjust(top=0.88)
fig.set_facecolor('white')

# Save
fig.savefig(output_dir / 'qc_summary.png')
