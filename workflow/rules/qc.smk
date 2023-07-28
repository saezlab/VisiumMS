rule sn_qc:
    input:
         "data/prc/sn/{sn_sample}/cellbender_matrix.h5"
    output:
         "data/prc/sn/{sn_sample}/adata.h5ad",
         "results/qc/sn_{sn_sample}.pdf"
    params:
        min_genes=config['qc']['min_genes'],
        min_cells=config['qc']['min_cells'],
        mt_thr=config['qc']['mt_thr'],
        ct_thr=config['qc']['ct_thr'],
        db_thr=config['qc']['db_thr']
    shell:
        """
        python workflow/scripts/qc/sn.py -i {input} -g {params.min_genes} -c {params.min_cells} -m {params.mt_thr} -n {params.ct_thr} -d {params.db_thr}
        """

rule merge_sn_qc:
    input:
        expand('results/qc/sn_{sn_sample}.pdf', sn_sample=sn_samples)
    output:
        'results/qc/sn_all.pdf'
    shell:
        """
        python workflow/scripts/general/merge_pdfs.py -i {input} -o {output}
        """

rule vm_qc:
    input:
         "data/raw/vs/{vs_sample}/"
    output:
         "data/raw/vs/{vs_sample}/adata.h5ad",
         "results/qc/vs_{vs_sample}.pdf"
    params:
        min_genes=config['qc']['min_genes'],
        min_cells=config['qc']['min_cells'],
        colors_dict=config['colors_areas']
    shell:
        """
        python workflow/scripts/qc/vs.py -i {input} -g {params.min_genes} -c {params.min_cells} -d '{params.colors_dict}'
        """

rule merge_vm_qc:
    input:
        expand('results/qc/vs_{vs_sample}.pdf', vs_sample=vs_samples)
    output:
        'results/qc/vs_all.pdf'
    shell:
        """
        python workflow/scripts/general/merge_pdfs.py -i {input} -o {output}
        """
