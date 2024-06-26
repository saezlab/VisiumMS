rule sn_merge:
    input:
        pdf='results/qc/sn_all.pdf',
        input_path='data/prc/sn/',
        meta='config/meta.csv'
    output:
        plot_path='results/integrate/sn_merge.pdf',
        out_path='data/prc/sn_merged.h5ad'
    params:
        n_hvg=config['sn_integrate']['n_hvg']
    resources:
        partition='cpu-multi',
        mem_mb=64000,
        slurm='ntasks-per-node=64'
    shell:
        """
        python workflow/scripts/integrate/sn_merge.py -i {input.input_path} -m {input.meta} -n {params.n_hvg} -p {output.plot_path} -o {output.out_path}
        """

rule sn_integrate:
    input:
        "data/prc/sn_merged.h5ad"
    output:
        plot_path='results/integrate/sn_integrate.pdf',
        out_path='data/prc/sn_integrated.h5ad'
    resources:
        partition='cpu-multi',
        slurm='ntasks-per-node=64'
    shell:
        """
        python workflow/scripts/integrate/sn_integrate.py -i {input} -p {output.plot_path} -o {output.out_path}
        """

rule sn_annotate:
    input:
        "data/prc/sn_integrated.h5ad"
    output:
        plot_path='results/integrate/sn_annotate.pdf',
        out_path='data/prc/sn_annotated.h5ad'
    params:
        markers=config['annotate']['markers'],
        resolution=config['annotate']['resolution'],
        annotation=config['annotate']['annotation']
    shell:
        """
        python workflow/scripts/integrate/sn_annotate.py -i {input} -m {params.markers} -r {params.resolution} -a '{params.annotation}' -p {output.plot_path} -o {output.out_path}
        """

rule niches_mofa:
    input:
        slide="data/prc/vs/{vs_sample}/adata.h5ad",
        pathway="data/prc/vs/{vs_sample}/hallmarks.csv",
        props="data/prc/vs/{vs_sample}/props.csv"
    output:
        plot='results/integrate/niches_mofa/{vs_sample}.pdf',
        model='data/prc/vs/{vs_sample}/niches_mofa.hdf5',
        out='data/prc/vs/{vs_sample}/niches.csv'
    params:
        n_hvg=config['sn_integrate']['n_hvg'],
        colors_dict=config['colors_areas'],
        annotation=lambda w: config['niches_ann'][w.vs_sample]
    resources:
        partition='cpu-multi',
        mem_mb=16000,
        slurm='ntasks-per-node=32'
    shell:
        """
        python workflow/scripts/integrate/niches_mofa.py -s {input.slide} -n {params.n_hvg} -m {output.model} -r 1.0 -c '{params.colors_dict}' -a '{params.annotation}' -p {output.plot} -o {output.out}
        """

rule cell_states:
    input:
        ann='data/prc/sn_annotated.h5ad'
    output:
        plot='results/integrate/ctypes/{ctype}.pdf',
        out='data/prc/ctypes/{ctype}.h5ad'
    resources:
        partition='cpu-multi',
        mem_mb=32000,
        slurm='ntasks-per-node=32'
    shell:
        """
        python workflow/scripts/integrate/cell_states.py -i {input} -p {output.plot} -o {output.out}
        """

rule cs_annotate:
    input:
        adata='data/prc/ctypes/{ctype}.h5ad',
        cstates='config/cellstates.csv'
    output:
        plot='results/integrate/cell_states/{ctype}.pdf',
        ann='data/prc/ctypes/{ctype}_ann.csv',
        deg='data/prc/ctypes/{ctype}_deg.csv'
    shell:
        """
        python workflow/scripts/integrate/cs_annotate.py \
        -d {input.adata} \
        -c {input.cstates} \
        -p {output.plot} \
        -a {output.ann} \
        -g {output.deg}
        """

rule save_data:
    input:
        sn='data/prc/sn_annotated.h5ad',
        cc=expand('data/prc/vs/{vs_sample}/ctlr_scores.csv', vs_sample=vs_samples),
        pr=expand('data/prc/vs/{vs_sample}/progeny.csv', vs_sample=vs_samples),
        hl=expand('data/prc/vs/{vs_sample}/hallmarks.csv', vs_sample=vs_samples),
        rc=expand('data/prc/vs/{vs_sample}/reactome.csv', vs_sample=vs_samples),
        co=expand('data/prc/ctypes/{ctype}.h5ad', ctype=ctypes),
        ca=expand('data/prc/ctypes/{ctype}_ann.csv', ctype=ctypes),
    output:
        d=directory('data/final/'),
        sn='data/final/sn_atlas.h5ad',
        vs=expand('data/final/visium_{vs_sample}.h5ad', vs_sample=vs_samples),
        ct=expand('data/final/ctype_{ctype}.h5ad', ctype=ctypes),
    shell:
        """
        python workflow/scripts/integrate/save_data.py \
        -a {input.sn} \
        -b {output.d}
        """



