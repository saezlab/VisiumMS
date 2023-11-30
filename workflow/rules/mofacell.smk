rule create_psbulks:
    input:
        i='data/prc/sn_annotated.h5ad',
        m='config/meta.csv',
        n=expand('data/prc/vs/{vs_sample}/niches.csv', vs_sample=vs_samples)
    output:
        a='data/prc/mofacell/sn_pb_data.csv',
        b='data/prc/mofacell/sn_pb_meta.csv',
        c='data/prc/mofacell/vs_pb_data.csv',
        d='data/prc/mofacell/vs_pb_meta.csv',       
    shell:
        """
        python workflow/scripts/mofacell/do_psbulk.py \
        -i {input.i} \
        -m {input.m} \
        -a {output.a} \
        -b {output.b} \
        -c {output.c} \
        -d {output.d}
        """

rule edger_markers:
    input:
        a='data/prc/mofacell/sn_pb_data.csv',
        b='data/prc/mofacell/sn_pb_meta.csv',
    output:
        'data/prc/mofacell/edgermarkers.csv'
    shell:
        """
        Rscript workflow/scripts/mofacell/edgermarkers.R {input.a} {input.b} {output}
        """

rule mofacell_sn_run:
    input:
        a='data/prc/mofacell/sn_pb_data.csv',
        b='data/prc/mofacell/sn_pb_meta.csv',
        c='config/meta.csv',
        d='data/prc/mofacell/edgermarkers.csv',
    output:
        e='data/prc/mofacell/sn_factors.csv',
        f='data/prc/mofacell/sn_loadings.csv',
        g='data/prc/mofacell/sn_model.hdf5',
        i='results/mofacell/sn_hmap.pdf',
        j='results/mofacell/sn_fspace.pdf',
        k='results/mofacell/sn_bplots.pdf',
        l='results/mofacell/sn_datacompleteness.pdf',
    shell:
        """
        Rscript workflow/scripts/mofacell/mofacell_sn_run.R \
        {input.a} {input.b} {input.c} {input.d} \
        {output.e} {output.f} {output.g} {output.i} {output.j} {output.k} {output.l}
        """

rule mofacell_vs_run:
    input:
        a='data/prc/mofacell/vs_pb_data.csv',
        b='data/prc/mofacell/vs_pb_meta.csv',
    output:
        c='data/prc/mofacell/vs_factors.csv',
        d='data/prc/mofacell/vs_loadings.csv',
        e='data/prc/mofacell/vs_model.hdf5',
        f='results/mofacell/vs_hmap.pdf',
        g='results/mofacell/vs_fspace.pdf',
        i='results/mofacell/vs_bplots.pdf',
        j='results/mofacell/vs_datacompleteness.pdf',
    shell:
        """
        Rscript workflow/scripts/mofacell/mofacell_vs_run.R \
        {input.a} {input.b} \
        {output.c} {output.d} {output.e} {output.f} {output.g} {output.i} {output.j}
        """

