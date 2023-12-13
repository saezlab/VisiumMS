import pandas as pd
import numpy as np
import decoupler as dc
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import matplotlib.colors
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--deg_path', required=True)
parser.add_argument('-b','--seg_path', required=True)
parser.add_argument('-c','--gmt_path', required=True)
parser.add_argument('-d','--plot_path', required=True)
args = vars(parser.parse_args())

deg_path = args['deg_path']
seg_path = args['seg_path']
gmt_path = args['gmt_path']
plot_path = args['plot_path']

# Read DEG at sn
deg = pd.read_csv(deg_path, index_col=0)
deg = deg[(deg['padj'] < 0.05) & (abs(deg['log2FoldChange']) > 1.)]
deg['condition'] = [c.split('vs')[0] if np.sign(s) > 0 else c.split('vs')[1] for c, s in zip(deg['contrast'], deg['stat'])]

# Read REACTOME
gmt = dc.read_gmt(gmt_path)
gmt['source'] = [s.split('REACTOME')[1].replace('_', ' ').lstrip() for s in gmt['source']]
msk = ~gmt['source'].str.contains('FETAL|INFECTION|SARS', case=False)
gmt = gmt[msk]
gmt = gmt.drop_duplicates(['source', 'target'])

# Compute ORA for ctype and condition
er = []
for ctype in ['OL', 'AS', 'MG']:
    for ltype in ['Ctrl', 'CA', 'CI']:
        print(ctype, ltype)
        ldeg = deg[(deg['cell_type'] == ctype) & (deg['condition'] == ltype)]
        ldeg = ldeg.reset_index().drop_duplicates(['index']).set_index('index')
        
        tmp = dc.get_ora_df(
            df=ldeg,
            net=gmt,
        )
        tmp['ctype'], tmp['ltype'] = ctype, ltype
        er.append(tmp)
er = pd.concat(er)

def get_act(er, ctype):
    padj = (
        er
        [er['ctype'] == ctype]
        .rename(columns={'ltype': 'index'})
        .pivot(columns='Term', index='index', values='FDR p-value')
        .fillna(1)
    )
    
    odds = (
        er
        [er['ctype'] == ctype]
        .rename(columns={'ltype': 'index'})
        .pivot(columns='Term', index='index', values='Odds ratio') # Overlap ratio, Odds ratio
        .fillna(0)
    ) / 100
    
    act = ad.AnnData(padj, dtype=float, layers={'odds': odds})
    act.obs['ltype'] = act.obs_names
    act.obs['ctype'] = ctype
    act.layers['score'] = - np.log10(act.X) * act.layers['odds']
    act.obs_names = [i + '_' + ctype for i in act.obs_names]
    return act

act = ad.concat([
    get_act(er, ctype='OL'),
    get_act(er, ctype='AS'),
    get_act(er, ctype='MG'),
], join='inner')
act.obs['join'] = act.obs['ltype'] + '_' + act.obs['ctype']

def plot_dotplot(act, var_names, mode='sn'):
    if mode == 'sn':
        categories_order = ['Ctrl_OL', 'CA_OL', 'CI_OL', 'Ctrl_AS', 'CA_AS', 'CI_AS',  'Ctrl_MG', 'CA_MG', 'CI_MG']
    else:
        categories_order = ['CA_PPWM', 'CA_LR', 'CA_LC', 'CA_VI', 'CI_PPWM', 'CI_LR', 'CI_LC', 'CI_VI']
    dot_color_df=act.to_df().loc[:, var_names]
    dot_size_df=dc.swap_layer(act, 'odds').to_df().loc[:, var_names]
    g = sc.pl.DotPlot(
        act,
        var_names,
        groupby='join',
        categories_order=categories_order,
        dot_size_df=dot_size_df,
    )
    g.swap_axes()
    g.style(cmap='viridis_r', dot_max=0.12)
    g.legend(size_title='Odds ratio', colorbar_title='FDR p-value')
    g.show()
    return g.fig

var_names=[
    'ATF6 ATF6 ALPHA ACTIVATES CHAPERONE GENES', 'FOLDING OF ACTIN BY CCT TRIC', 'HSF1 DEPENDENT TRANSACTIVATION',  # CA_OL    
    'ANTIGEN PRESENTATION FOLDING ASSEMBLY AND PEPTIDE LOADING OF CLASS I MHC', 'ENDOSOMAL VACUOLAR PATHWAY',  # CA_OL
    'GP1B IX V ACTIVATION SIGNALLING', 'PLATELET ADHESION TO EXPOSED COLLAGEN', 'NUCLEAR SIGNALING BY ERBB4', 'SIGNALING BY ERBB2',  # CI_AS
    'EPHRIN SIGNALING', 'INTERLEUKIN RECEPTOR SHC SIGNALING', 'ACTIVATION OF C3 AND C5', 'SIGNALING BY VEGF',  # Ctrl_MG
    'SELENOAMINO ACID METABOLISM', 'SIGNALING BY ROBO RECEPTORS', 'CELLULAR RESPONSE TO STARVATION',  # CA_MG and CI_MG
]

fig1 = plot_dotplot(act, var_names=var_names)

# Read niches DEG at ST
seg = pd.read_csv(seg_path, index_col=0)
seg = seg.reset_index().pivot(index='niche', columns='index', values='stat').fillna(0)

# Run ulm enrichment
act, pvl = dc.dense_run(
    func=dc.run_ulm,
    mat=seg,
    net=gmt,
    weight=None,
)
act = ad.AnnData(act.fillna(0))
act.obs['join'] = act.obs_names
pvl = pvl.fillna(0)

# Plot
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#8be9f1","white","#ffa4e8"])
var_names=[
    'ADAPTIVE IMMUNE SYSTEM',
    'CYTOKINE SIGNALING IN IMMUNE SYSTEM',
    'HSF1 DEPENDENT TRANSACTIVATION',
    'SIGNALING BY TGF BETA RECEPTOR COMPLEX IN CANCER',
    'LOSS OF FUNCTION OF SMAD2 3 IN CANCER',
    'CARGO TRAFFICKING TO THE PERICILIARY MEMBRANE',
    'SMOOTH MUSCLE CONTRACTION',
    'TYPE I HEMIDESMOSOME ASSEMBLY',
    'LAMININ INTERACTIONS',
    'BIOTIN TRANSPORT AND METABOLISM',
    'EGR2 AND SOX10 MEDIATED INITIATION OF SCHWANN CELL MYELINATION',
    'CITRIC ACID CYCLE TCA CYCLE',
    'SELENOAMINO ACID METABOLISM', 'CELLULAR RESPONSE TO STARVATION',
    'RESPIRATORY ELECTRON TRANSPORT',
    'COMPLEX I BIOGENESIS',
]

fig2 = sc.pl.matrixplot(
    adata=act,
    var_names=var_names,
    groupby='join',
    cmap=cmap,
    vmax=3,
    vmin=-3,
    swap_axes=True,
    categories_order=['PPWM', 'LR', 'LC', 'VI'],
    colorbar_title='Activity (CA-, CI+)',
    return_fig=True,
    show=False
)
fig2.show()
fig2 = fig2.fig

# Save to pdf
pdf = matplotlib.backends.backend_pdf.PdfPages(plot_path)
for fig in [fig1, fig2]:
    pdf.savefig(fig, bbox_inches='tight')
pdf.close()
