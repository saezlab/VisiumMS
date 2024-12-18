rule sn_deg:
    input:
        'data/prc/sn_annotated.h5ad'
    output:
        'data/prc/sn_deg.csv'
    resources:
        partition='cpu-multi',
        slurm='ntasks-per-node=64'        
    shell:
        """
        python workflow/scripts/func/sn_deg.py -i {input} -o {output}
        """

rule vs_deg:
    input:
        m='config/meta.csv',
        n=expand("data/prc/vs/{vs_sample}/niches.csv", vs_sample=vs_samples)
    output:
        'data/prc/vs_deg.csv'
    resources:
        partition='cpu-multi',
        slurm='ntasks-per-node=64'        
    shell:
        """
        python workflow/scripts/func/vs_deg.py -m {input.m} -o {output}
        """

rule ns_deg:
    input:
        m='config/meta.csv',
        n=expand("data/prc/vs/{vs_sample}/niches.csv", vs_sample=vs_samples)
    output:
        'data/prc/ns_deg.csv'
    resources:
        partition='cpu-multi',
        slurm='ntasks-per-node=64'        
    shell:
        """
        python workflow/scripts/func/ns_deg.py -m {input.m} -o {output}
        """

rule summary_deg:
    input:
        sn='data/prc/sn_deg.csv',
        ns='data/prc/ns_deg.csv',
    output:
        'results/deg/summary.pdf'
    params:
        cdict=config['colors_conds'],
    shell:
        """
        python workflow/scripts/func/summary_deg.py -a {input.sn} -b {input.ns} -c '{params.cdict}' -d {output}
        """

rule sn_pathway:
    input:
        inp='data/prc/sn_deg.csv',
        gmt='config/c2.cp.reactome.v2023.1.Hs.symbols.gmt'
    output:
        'data/prc/sn_pathway.csv'
    shell:
        """
        python workflow/scripts/func/sn_pathway.py -i {input.inp} -g {input.gmt} -o {output}
        """

rule vs_traj:
    input:
        m='config/meta.csv',
        n=expand("data/prc/vs/{vs_sample}/niches.csv", vs_sample=vs_samples),
    output:
        o='data/prc/vs_traj.csv',
        p='results/traj/vs_traj.pdf'
    shell:
        """
        python workflow/scripts/func/vs_traj.py -m {input.m} -o {output.o} -p {output.p}
        """

rule sn_traj:
    input:
        a='data/prc/sn_annotated.h5ad',
        n='config/c2.cp.reactome.v2023.1.Hs.symbols.gmt',
        r='data/prc/sn_pathway.csv'
    output:
        o='data/prc/sn_traj.csv',
        p='results/traj/sn_traj.pdf'
    shell:
        """
        python workflow/scripts/func/sn_traj.py -a {input.a} -n {input.n} -r {input.r} -o {output.o} -p {output.p}
        """

rule vs_pathway_hallmarks:
    input:
        inp='data/prc/vs/{vs_sample}/adata.h5ad',
        gmt='config/h.all.v2023.1.Hs.symbols.gmt'
    output:
        'data/prc/vs/{vs_sample}/hallmarks.csv'
    shell:
        """
        python workflow/scripts/func/vs_pathway.py -s {input.inp} -g {input.gmt} -n 'HALLMARK' -o {output}
        """

rule vs_pathway_reactome:
    input:
        inp='data/prc/vs/{vs_sample}/adata.h5ad',
        gmt='config/c2.cp.reactome.v2023.1.Hs.symbols.gmt'
    output:
        'data/prc/vs/{vs_sample}/reactome.csv'
    shell:
        """
        python workflow/scripts/func/vs_pathway.py -s {input.inp} -g {input.gmt} -n 'REACTOME' -o {output}
        """

rule vs_pathway_progeny:
    input:
        inp='data/prc/vs/{vs_sample}/adata.h5ad',
        gmt='config/progeny.csv'
    output:
        'data/prc/vs/{vs_sample}/progeny.csv'
    shell:
        """
        python workflow/scripts/func/vs_pathway.py -s {input.inp} -g {input.gmt} -n 'progeny' -o {output}
        """

rule coda_compositions:
    input:
        meta='config/meta.csv',
        ann='data/prc/sn_annotated.h5ad',
        cs='config/cellstates.csv',
        deg=expand("data/prc/ctypes/{ctype}_deg.csv", ctype=ctypes)
    output:
        plot='results/composition/props.pdf',
        table='results/composition/table.csv'
    resources:
        partition='cpu-multi',
        slurm='ntasks-per-node=64',
        runtime=120
    shell:
        """
        python workflow/scripts/func/compositions.py -m {input.meta} -a {input.ann} -s {input.cs} -p {output.plot} -t {output.table}
        """

rule clr_comps:
    input:
        ann='data/prc/sn_annotated.h5ad',
        meta='config/meta.csv',
        cs=expand("data/prc/ctypes/{ctype}_ann.csv", ctype=ctypes)
    output:
        k='data/prc/comps/krustal_table.csv',
        w='data/prc/comps/wilcoxon_table.csv',
        d='data/prc/comps/dfs.csv',
    shell:
        """
        python workflow/scripts/func/clr_comps.py \
        -m {input.meta} -a {input.ann} \
        -k {output.k} -w {output.w} -d {output.d}
        """

rule cs_traj:
    input:
        ann='data/prc/sn_annotated.h5ad',
        d='data/prc/sn_deg.csv',
        k='data/prc/comps/krustal_table.csv',
        w='data/prc/comps/wilcoxon_table.csv',
        cs=expand("data/prc/ctypes/{ctype}_deg.csv", ctype=ctypes)
    output:
        'data/prc/cs_traj.csv'
    shell:
        """
        python workflow/scripts/func/cs_traj.py -a {input.ann} -d {input.d} -k {input.k}  -w {input.w} -o {output}
        """

rule sn_lr:
    input:
        deg='data/prc/sn_deg.csv',
        ann='data/prc/sn_annotated.h5ad'
    output:
        plot='results/ccc/sn_lr.pdf',
        df='data/prc/sn_lr.csv'
    resources:
        mem_mb=32000
    shell:
        """
        python workflow/scripts/func/sn_lr.py -d {input.deg} -a {input.ann} -p {output.plot} -o {output.df} -s 0.15 -g 0.05
        """

rule vs_ctlr:
    input:
        slide='data/prc/vs/{vs_sample}/adata.h5ad',
        props='data/prc/vs/{vs_sample}/props.csv',
        sn_lr='data/prc/sn_lr.csv'
    output:
        lr='data/prc/vs/{vs_sample}/ctlr_scores.csv'
    shell:
        """
        python workflow/scripts/func/vs_ctlr.py -s {input.slide} -p {input.props} -n {input.sn_lr} -b 150 -o {output.lr}
        """

rule test_ctlr:
    input:
        slides=expand('data/prc/vs/{vs_sample}/ctlr_scores.csv', vs_sample=vs_samples),
        sn_lr='data/prc/sn_lr.csv',
        meta='config/meta.csv'
    output:
        'data/prc/vs_diff_ctlr.csv'
    shell:
        """
        python workflow/scripts/func/test_ctlr.py -n {input.sn_lr} -m {input.meta} -p 0.15 -o {output}
        """

rule corr_ctlr_pw:
    input:
        c='data/prc/vs_diff_ctlr.csv',
        p='data/prc/sn_pathway.csv',
        m='config/meta.csv'
    output:
        'data/prc/corr_ctlr_pw.csv'
    shell:
        """
        python workflow/scripts/func/corr_ctlr_pw.py -c {input.c} -p {input.p} -m {input.m} -t 0.10 -a 0.15 -o {output}
        """

rule corr_pw_scores:
    input:
        a='data/prc/cs_ctlr.csv',
        b='config/meta.csv'
    output:
        'data/prc/corr_pw_scores.csv'
    shell:
        """
        python workflow/scripts/func/corr_pw_scores.py -a {input.a} -b {input.b} -c {output}
        """

rule cs_ctlr:
    input:
        g='config/markers.csv',
        t='data/prc/comps/wilcoxon_table.csv',
        s='data/prc/sn_lr.csv',
        c='data/prc/vs_diff_ctlr.csv',
        m='config/meta.csv',
        states=expand('data/prc/ctypes/{ctype}_deg.csv', ctype=ctypes),
    output:
        o='data/prc/cs_ctlr.csv',
        p='results/ccc/cs_ctlr.pdf'
    params:
        cdict=config['colors_conds'],
    shell:
        """
        python workflow/scripts/func/cs_ctlr.py \
        -g {input.g} \
        -t {input.t} \
        -s {input.s} \
        -c {input.c} \
        -m {input.m} \
        -o {output.o} \
        -q '{params.cdict}' \
        -p {output.p}
        """

rule run_misty:
    input:
        sn_lr='data/prc/sn_lr.csv',
        sn_pw='data/prc/sn_pathway.csv',
        slide='data/prc/vs/{vs_sample}/adata.h5ad',
        vs_lr='data/prc/vs/{vs_sample}/lr_scores.csv',
        props='data/prc/vs/{vs_sample}/props.csv',
        react='data/prc/vs/{vs_sample}/reactome.csv',
    output:
        inters='data/prc/vs/{vs_sample}/misty_inters.csv',
        metrics='data/prc/vs/{vs_sample}/misty_metrics.csv'
    shell:
        """
        python workflow/scripts/func/run_misty.py \
        -d {input.sn_lr} \
        -w {input.sn_pw} \
        -s {input.slide} \
        -l {input.vs_lr} \
        -p {input.props} \
        -r {input.react} \
        -a 0.15 \
        -b 100 \
        -i {output.inters} \
        -t {output.metrics}
        """

