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

# CA
tmp = df[df['type'] == 'CA'].copy()
fig1, ax = plt.subplots(1, 1, figsize=(2, 2), dpi=150)
sns.boxplot(
    data=tmp,
    x='group',
    y='asc_mm2',
    ax=ax,
    order=['Border', 'LC']
)
ax.set_title('CA')

s, p = ranksums(
    x=tmp.loc[tmp['group'] == 'LC', 'asc_mm2'],
    y=tmp.loc[tmp['group'] == 'Border', 'asc_mm2']
)
print('asc_mm2: ', p)


fig2, ax = plt.subplots(1, 1, figsize=(2, 2), dpi=150)
sns.boxplot(
    data=tmp,
    x='group',
    y='asc_prop',
    ax=ax,
    order=['Border', 'LC']
)
ax.set_title('CA')
s, p = ranksums(
    x=tmp.loc[tmp['group'] == 'LC', 'asc_prop'],
    y=tmp.loc[tmp['group'] == 'Border', 'asc_prop']
)
print('asc_prop: ', p)

# CI
tmp = df[df['type'] == 'CI'].copy()
fig3, ax = plt.subplots(1, 1, figsize=(2, 2), dpi=150)
sns.boxplot(
    data=tmp,
    x='group',
    y='asc_mm2',
    ax=ax,
    order=['Border', 'LC']
)
ax.set_title('CI')


s, p = ranksums(
    x=tmp.loc[tmp['group'] == 'LC', 'asc_mm2'],
    y=tmp.loc[tmp['group'] == 'Border', 'asc_mm2']
)
print('asc_mm2: ', p)


fig4, ax = plt.subplots(1, 1, figsize=(2, 2), dpi=150)
sns.boxplot(
    data=tmp,
    x='group',
    y='asc_prop',
    ax=ax,
    order=['Border', 'LC']
)
ax.set_title('CI')
s, p = ranksums(
    x=tmp.loc[tmp['group'] == 'LC', 'asc_prop'],
    y=tmp.loc[tmp['group'] == 'Border', 'asc_prop']
)
print('asc_prop: ', p)

# Save to pdf
pdf = matplotlib.backends.backend_pdf.PdfPages(output_path)
for fig in [fig1, fig2, fig3, fig4]:
    pdf.savefig(fig, bbox_inches='tight')
pdf.close()
