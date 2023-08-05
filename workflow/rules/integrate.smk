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
        pathway="data/prc/vs/{vs_sample}/pathway.csv",
        props="data/prc/vs/{vs_sample}/props.csv"
    output:
        plot='results/integrate/vs_{vs_sample}_mofa.pdf',
        model='data/prc/vs/{vs_sample}/niches_mofa.hdf5'
    params:
        n_hvg=config['sn_integrate']['n_hvg'],
        colors_dict=config['colors_areas']
    resources:
        partition='cpu-multi',
        mem_mb=16000,
        slurm='ntasks-per-node=32'
    shell:
        """
        python workflow/scripts/integrate/niches_mofa.py -s {input.slide} -n {params.n_hvg} -m {output.model} -r 1.0 -c '{params.colors_dict}' -p {output.plot}
        """
