rule vs_hallmarks:
    input:
        'data/prc/vs/{vs_sample}/adata.h5ad'
    output:
        'data/prc/vs/{vs_sample}/hallmarks.csv'
    params:
        hallm_path=config['hallmarks']
    shell:
        """
        python workflow/scripts/func/vs_hallmarks.py -s {input} -g {params.hallm_path} -o {output}
        """
