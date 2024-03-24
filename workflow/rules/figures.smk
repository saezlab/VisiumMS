rule fig1_spatial_features:
    input:
        expand("data/prc/vs/{vs_sample}/niches.csv", vs_sample=vs_samples)
    output:
        'results/figures/fig1/spatial_features.pdf'
    params:
        colors_dict=config['colors_areas']
    shell:
        """
        python workflow/scripts/figures/fig1/spatial_features.py -d '{params.colors_dict}' -p {output}
        """

rule fig1_ncstates:
    input:
        x=expand("data/prc/ctypes/{ctype}_ann.csv", ctype=ctypes),
        a='config/cellstates.csv'
    output:
        'results/figures/fig1/ncstates.pdf'
    params:
        colors_dict=config['colors_ctypes']
    shell:
        """
        python workflow/scripts/figures/fig1/ncstates.py -a {input.a} -b '{params.colors_dict}' -c {output}
        """

rule fig1_niches_markers:
    input:
        m='config/meta.csv',
        n=expand("data/prc/vs/{vs_sample}/niches.csv", vs_sample=vs_samples)
    output:
        'results/figures/fig1/niches_markers.pdf'
    shell:
        """
        python workflow/scripts/figures/fig1/niches_markers.py -m '{input.m}' -p {output}
        """

rule fig2_asc:
    input:
        o='data/prc/ctypes/AS.h5ad',
        a='data/prc/ctypes/AS_ann.csv',
        m='config/meta.csv',
        r='config/c2.cp.reactome.v2023.1.Hs.symbols.gmt',
        c='config/collectri.csv',
    output:
        plot='results/figures/fig2/asc.pdf',
        csv='results/figures/fig2/asc.csv',
    shell:
        """
        python workflow/scripts/figures/fig2/asc.py \
        -o '{input.o}' \
        -a '{input.a}' \
        -m '{input.m}' \
        -r '{input.r}' \
        -c '{input.c}' \
        -m '{input.m}' \
        -d '{output.csv}' \
        -p {output.plot}
        """

rule fig2_ascquant:
    input:
        'config/asc_counts.csv',
    output:
        'results/figures/fig2/ascquant.pdf'
    shell:
        """
        python workflow/scripts/figures/fig2/quant.py -i {input} -o {output}
        """


rule fig2_ependym:
    input:
        a='data/prc/vs/MS549H/adata.h5ad',
        b='data/prc/vs/MS549T/adata.h5ad',
        c='data/prc/vs/MS549H/niches.csv',
        d='data/prc/vs/MS549T/niches.csv',
    output:
        'results/figures/fig2/ependym.pdf'
    shell:
        """
        python workflow/scripts/figures/fig2/ependym.py -a {input.a} -b {input.b} -c {input.c} -d {input.d} -e {output}
        """

rule fig3_mds:
    input:
        a=expand("data/prc/vs/{vs_sample}/niches.csv", vs_sample=vs_samples),
        b='data/prc/sn_annotated.h5ad',
        c='config/meta.csv',
        d='data/prc/mofacell/sn_factors.csv',
        e='data/prc/mofacell/vs_factors.csv',
    output:
        'results/figures/fig3/mds.pdf'
    params:
        cdict=config['colors_conds']
    shell:
        """
        python workflow/scripts/figures/fig3/mds.py \
        -a '{params.cdict}' -b {input.b} -c {input.c} -d {input.d} -e {input.e} -f {output}
        """

rule fig3_pwheatmaps:
    input:
        a='data/prc/sn_deg.csv',
        b='data/prc/ns_deg.csv',
        c='config/c2.cp.reactome.v2023.1.Hs.symbols.gmt',
    output:
        'results/figures/fig3/pwheatmaps.pdf'
    shell:
        """
        python workflow/scripts/figures/fig3/pwheatmaps.py \
        -a {input.a} -b {input.b} -c {input.c} -d {output}
        """

rule fig3_venn:
    input:
        sn='data/prc/sn_deg.csv',
    output:
        'results/figures/fig3/venn.pdf'
    params:
        cdict=config['colors_conds'],
    shell:
        """
        python workflow/scripts/figures/fig3/venn.py -a {input.sn} -b '{params.cdict}' -c {output}
        """

rule fig3_gmarkers:
    input:
        a='data/prc/sn_annotated.h5ad',
        b='data/prc/cs_traj.csv',
    output:
        'results/figures/fig3/gmarkers.pdf'
    shell:
        """
        python workflow/scripts/figures/fig3/gmarkers.py \
        -a {input.a} -b {input.b} -c {output}
        """

rule fig3_clrviolins:
    input:
        a='data/prc/comps/dfs.csv',
        b='data/prc/comps/wilcoxon_table.csv',
    output:
        'results/figures/fig3/clrviolins.pdf'
    shell:
        """
        python workflow/scripts/figures/fig3/clrviolins.py \
        -a {input.a} -b {input.b} -c {output}
        """

rule fig4_summary:
    input:
        a='data/prc/cs_ctlr.csv',
    output:
        'results/figures/fig4/summary.pdf'
    params:
        cdict=config['colors_conds'],
    shell:
        """
        python workflow/scripts/figures/fig4/summary.py -a {input} -b '{params.cdict}' -c {output}
        """

rule fig4_examples:
    input:
        a='data/prc/cs_ctlr.csv',
        b='data/prc/sn_annotated.h5ad',
        c='data/prc/corr_pw_scores.csv',
        d='config/meta.csv',
    output:
        'results/figures/fig4/examples.pdf'
    params:
        cdict=config['colors_conds'],
    shell:
        """
        python workflow/scripts/figures/fig4/examples.py -a {input.a} -b {input.b} -c {input.c} -d {input.d} -e '{params.cdict}' -f {output}
        """




