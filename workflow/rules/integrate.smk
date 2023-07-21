rule sn_merge:
    input:
        pdf='results/qc/sn_all.pdf',
        input_path='data/prc/sn/',
        meta='config/meta.csv'
    output:
        plot_path='results/integrate/sn_merge.pdf',
        out_path='data/prc/sn/merged.h5ad'
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
        "data/prc/sn/merged.h5ad"
    output:
        plot_path='results/integrate/sn_integrate.pdf',
        out_path='data/prc/sn/integrated.h5ad'
    resources:
        partition='cpu-multi',
        slurm='ntasks-per-node=64'
    shell:
        """
        python workflow/scripts/integrate/sn_integrate.py -i {input} -p {output.plot_path} -o {output.out_path}
        """
