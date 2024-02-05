import pandas as pd
import matplotlib.pyplot as plt
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--cstates_path', required=True)
parser.add_argument('-b','--colors_dict', required=True)
parser.add_argument('-c','--plot_path', required=True)
args = vars(parser.parse_args())

cstates_path = args['cstates_path']
colors_dict = args['colors_dict']
plot_path = args['plot_path']

# Define palette
palette = dict(item.split(':') for item in colors_dict.strip("'").split(';'))

# Read and remove NA states
df = pd.read_csv(cstates_path)
df = df[~df['name'].str.contains('_NA')]

# Count n states per ctype
df[['ctype', 'cstate']] = df['name'].str.split('_', n=1, expand=True)
df = df.groupby('ctype').count()[['name']]

# Plot
order = ['OL', 'BC', 'MG', 'TC', 'NEU', 'OPC', 'AS', 'EC', 'SC']
df = df.loc[order]
fig, ax = plt.subplots(1, 1, figsize=(2, 4), dpi=150)
bars = ax.barh(df.index, df['name'], color=[palette[ctype] for ctype in df.index])
for bar, value in zip(bars, df['name']):
    ax.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height() / 2, str(value), ha='left', va='center')
ax.set_ylabel('')
ax.set_xlabel('n of cell states')
ax.set_xlim(0, 17)

fig.savefig(plot_path, bbox_inches='tight')
