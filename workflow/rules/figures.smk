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
        'results/figures/fig2/asc.pdf'
    shell:
        """
        python workflow/scripts/figures/fig2/asc.py \
        -o '{input.o}' \
        -a '{input.a}' \
        -m '{input.m}' \
        -r '{input.r}' \
        -c '{input.c}' \
        -m '{input.m}' \
        -p {output}
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
