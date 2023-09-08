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

rule compositions:
    input:
        meta='config/meta.csv',
        ann='data/prc/sn_annotated.h5ad'
    output:
        plot='results/composition/props.pdf',
        table='results/composition/tests.csv'
    shell:
        """
        python workflow/scripts/func/compositions.py -m {input.meta} -a {input.ann} -p {output.plot} -t {output.table}
        """

rule sn_lr:
    input:
        deg='data/prc/sn_deg.csv',
        ann='data/prc/sn_annotated.h5ad'
    output:
        plot='results/ccc/sn_lr.csv',
        df='data/prc/sn_lr.csv'
    shell:
        """
        python workflow/scripts/func/sn_lr.py -d {input.deg} -a {input.ann} -p {output.plot} -o {output.df} -s 0.15 -g 0.05
        """

rule vs_lr:
    input:
        slide='data/prc/vs/{vs_sample}/adata.h5ad'
    output:
        lr='data/prc/vs/{vs_sample}/lr_scores.csv'
    shell:
        """
        python workflow/scripts/func/vs_lr.py -s {input.slide} -t 0.05 -b 100 -o {output.lr}
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
