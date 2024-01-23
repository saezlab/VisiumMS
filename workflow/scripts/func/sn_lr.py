import numpy as np
import pandas as pd
import scanpy as sc
import decoupler as dc
import plotnine as p9
import liana as li
import matplotlib.backends.backend_pdf
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-d','--deg_path', required=True)
parser.add_argument('-a','--ann_path', required=True)
parser.add_argument('-p','--plot_path', required=True)
parser.add_argument('-o', '--out_path', required=True)
parser.add_argument('-s','--thr_adjpval', required=True)
parser.add_argument('-g','--thr_gexprop', required=True)
args = vars(parser.parse_args())

deg_path = args['deg_path']
ann_path = args['ann_path']
plot_path = args['plot_path']
out_path = args['out_path']
thr_adjpval = float(args['thr_adjpval'])
thr_gexprop = float(args['thr_gexprop'])

# Get DEG results
deg = pd.read_csv(deg_path, index_col=0)
deg = deg.rename(columns={'cell_type': 'leiden'})  # Need same name as in adata
contrasts = deg['contrast'].unique().astype(str)

# Read atlas
adata = sc.read_h5ad(ann_path)

# For each contrast
df_lr = []
figs = []
for contrast in contrasts:

    print(contrast)

    # Subset deg
    cdeg = deg[(deg['contrast'] == contrast)].copy()

    # Compute LR scores based on DEG results
    c_lr = li.multi.df_to_lr(
        adata,
        dea_df=cdeg,
        resource_name='consensus',
        expr_prop=thr_gexprop, # calculated for adata as passed - used to filter interactions
        groupby='leiden',
        stat_keys=['stat', 'pvalue', 'padj'],
        use_raw=False,
        complex_col='stat', # NOTE: we use the Wald Stat to deal with complexes
        verbose=False,
        return_all_lrs=False,
    )

    # Format df
    c_lr = c_lr[['source', 'target', 'ligand_complex', 'receptor_complex', 'interaction', 'ligand_stat', 'ligand_padj',
                     'receptor_stat', 'receptor_padj', 'interaction_stat', 'interaction_padj']]
    
    # Remove interactions with missmatched sign and not significant
    c_lr = (
        c_lr[c_lr['interaction_padj'] < thr_adjpval]
        .sort_values('interaction_padj')
        .assign(both=lambda x: ((x['ligand_padj'] < thr_adjpval) & (x['receptor_padj'] < thr_adjpval)) & \
                (np.sign(x['ligand_stat']) == np.sign(x['receptor_stat'])))
    )
    c_lr = (
        c_lr[c_lr['both']]
        .reset_index(drop=True)
        .drop(columns='both')
        .assign(sign=lambda x: np.sign(x['ligand_stat']))
    )
    c_lr.insert(0, 'contrast', contrast)

    # Compute LR scores based on expr for cell types with no contrast
    cond, ref = contrast.split('vs')
    if ref == 'Ctrl':
        # Subset adata
        sub_adata = adata[adata.obs['Lesion type'] == cond].copy()

        # Run rank_aggregate
        li.mt.rank_aggregate(
            sub_adata,
            groupby='leiden',
            expr_prop=thr_gexprop,
            verbose=True,
            use_raw=False,
            n_perms=None
        )

        # Extract results and format
        df = sub_adata.uns['liana_res'].copy()
        df['interaction_padj'] = dc.p_adjust_fdr(df['magnitude_rank'])
        df = df[((df['source'].str.contains('BC|TC|SC')) | (df['target'].str.contains('BC|TC|SC'))) & (df['interaction_padj'] < thr_adjpval)]
        df['contrast'] = contrast
        df['interaction'] = [x + '^' + y for x, y in zip(df['ligand_complex'], df['receptor_complex'])]
        df['ligand_stat'] = df['expr_prod']
        df['receptor_stat'] = df['expr_prod']
        df['interaction_stat'] = df['expr_prod']
        df['ligand_padj'] = df['interaction_padj']
        df['receptor_padj'] = df['interaction_padj']
        df['sign'] = +1
        df = df[['contrast', 'source', 'target', 'ligand_complex', 'receptor_complex', 'interaction', 'ligand_stat',
                 'ligand_padj', 'receptor_stat', 'receptor_padj', 'interaction_stat', 'interaction_padj', 'sign']]
        c_lr = pd.concat([c_lr, df])
        # Remove duplicates
        c_lr = c_lr.drop_duplicates(['contrast', 'source', 'target', 'interaction'])

    # Plot
    fig = li.pl.tileplot(
        liana_res=c_lr,
        fill = 'interaction_stat',
        label='padj',
        label_fun = lambda x: '*' if x < thr_adjpval else np.nan,
        top_n=50,
        orderby = 'interaction_padj',
        orderby_ascending = True,
        orderby_absolute = False,
        cmap='coolwarm',
        figure_size=(7, 14),
        return_fig=True
    ) + p9.ggtitle(contrast)

    # Store
    df_lr.append(c_lr)
    figs.append(fig)
df_lr = pd.concat(df_lr)

# Save to pdf
p9.save_as_pdf_pages(plots=figs, filename=plot_path)

# Write csv
df_lr.to_csv(out_path, index=False)
