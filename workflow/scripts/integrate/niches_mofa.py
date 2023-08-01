import os
import pandas as pd
import numpy as np
from mofapy2.run.entry_point import entry_point
import mofax as mfx
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-m','--meta_path', required=True)
parser.add_argument('-n','--n_factors', required=True)
parser.add_argument('-p','--plot_path', required=True)
parser.add_argument('-o', '--out_path', required=True)
args = vars(parser.parse_args())

meta_path = args['meta_path']
n_factors = int(args['n_factors'])
plot_path = args['plot_path']
out_path = args['out_path']

# Read meta
meta = pd.read_csv(meta_path)
vs_samples = meta[~meta['Batch vs'].isnull()]['Sample id'].values.astype('U')

# Process hallmarks and abunds
df = []
for sample_id in vs_samples:
    # Read data
    hallmarks = pd.read_csv('data/prc/vs/{0}/hallmarks.csv'.format(sample_id), index_col=0)
    abunds = pd.read_csv('data/prc/vs/{0}/abunds.csv'.format(sample_id), index_col=0)
    hallmarks.index = [sample_id + '_' + i for i in hallmarks.index]
    abunds.index = [sample_id + '_' + i for i in abunds.index]

    # Transform to melted df
    abunds = (
        abunds
        .reset_index(names='sample')
        .melt(id_vars='sample', var_name="feature", value_name="value")
        .assign(view=lambda x: 'abunds', group=lambda x: sample_id)
        [["sample", "group", "feature", "value", "view"]]
    )
    hallmarks = (
        hallmarks
        .reset_index(names='sample')
        .melt(id_vars='sample', var_name="feature", value_name="value")
        .assign(view=lambda x: 'hallmarks', group=lambda x: sample_id)
        [["sample", "group", "feature", "value", "view"]]
    )

    # Append
    df.append(abunds)
    df.append(hallmarks)
df = pd.concat(df)

# MOFA+
ent = entry_point()
ent.set_data_options(
    scale_views = True
)
ent.set_data_df(df, likelihoods = ['gaussian', 'gaussian'])
ent.set_model_options(
    factors = n_factors,
    spikeslab_weights = True,
    ard_weights = True,
    ard_factors = True
)
ent.set_train_options(
    convergence_mode = "fast",
    dropR2 = 0.001,
    gpu_mode = False,
    seed = 1
)
ent.build()
ent.run()

# Write and read model
ent.save(outfile=out_path)
model = mfx.mofa_model(out_path)

# Save plot
fig = mfx.plot.plot_r2(model, y="Group", x="Factor")
fig.tight_layout()
fig.savefig(plot_path, bbox_inches='tight')

