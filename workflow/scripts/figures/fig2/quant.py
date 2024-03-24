import pandas as pd
import numpy as np
import scipy
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
from scipy.stats import ranksums
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--input_path', required=True)
parser.add_argument('-o','--output_path', required=True)
args = vars(parser.parse_args())

input_path = args['input_path']
output_path = args['output_path']

df = pd.read_csv(input_path)

# Plot and test
fig1, ax = plt.subplots(1, 1, figsize=(2, 2), dpi=150)
sns.boxplot(
    data=df,
    x='group',
    y='asc_mm2',
    ax=ax
)

s, p = ranksums(
    x=df.loc[df['group'] == 'LC', 'asc_mm2'],
    y=df.loc[df['group'] == 'rest', 'asc_mm2']
)
print('asc_mm2: ', p)

# Plot and test
fig2, ax = plt.subplots(1, 1, figsize=(2, 2), dpi=150)
sns.boxplot(
    data=df,
    x='group',
    y='asc_prop',
    ax=ax
)
s, p = ranksums(
    x=df.loc[df['group'] == 'LC', 'asc_prop'],
    y=df.loc[df['group'] == 'rest', 'asc_prop']
)
print('asc_prop: ', p)

# Save to pdf
pdf = matplotlib.backends.backend_pdf.PdfPages(output_path)
for fig in [fig1, fig2]:
    pdf.savefig(fig, bbox_inches='tight')
pdf.close()
