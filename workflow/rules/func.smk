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

rule vs_pathway:
    input:
        inp='data/prc/vs/{vs_sample}/adata.h5ad',
        gmt='config/c2.cp.reactome.v2023.1.Hs.symbols.gmt'
    output:
        'data/prc/vs/{vs_sample}/pathway.csv'
    shell:
        """
        python workflow/scripts/func/vs_pathway.py -s {input.inp} -g {input.gmt} -o {output}
        """
