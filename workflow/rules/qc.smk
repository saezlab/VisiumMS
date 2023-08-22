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

rule vs_qc:
    input:
         sample="data/raw/vs/{vs_sample}/",
         meta='config/meta.csv'
    output:
         "data/prc/vs/{vs_sample}/adata.h5ad",
         "results/qc/vs_{vs_sample}.pdf"
    params:
        min_genes=config['qc']['min_genes'],
        min_cells=config['qc']['min_cells'],
        colors_dict=config['colors_areas']
    shell:
        """
        python workflow/scripts/qc/vs.py -i {input.sample} -m {input.meta} -g {params.min_genes} -c {params.min_cells} -d '{params.colors_dict}'
        """

rule merge_vs_qc:
    input:
        expand('results/qc/vs_{vs_sample}.pdf', vs_sample=vs_samples)
    output:
        'results/qc/vs_all.pdf'
    shell:
        """
        python workflow/scripts/general/merge_pdfs.py -i {input} -o {output}
        """

rule summary:
    input:
        sn_inp_path='data/prc/sn_annotated.h5ad',
        meta='config/meta.csv',
        sn_qc='results/qc/sn_all.pdf',
        vs_qc='results/qc/vs_all.pdf'
    output:
        'results/qc/summary.pdf'
    shell:
        """
        python workflow/scripts/qc/summary.py -s {input.sn_inp_path} -v 'data/prc/vs/' -m {input.meta} -p {output}
        """

rule deconv_qc:
    input:
        meta_path='config/meta.csv',
        ann_path='data/prc/sn_annotated.h5ad',
        samples=expand('data/prc/vs/{vs_sample}/props.csv', vs_sample=vs_samples)
    output:
        leiden='results/qc/deconv_leiden.pdf',
        sample='results/qc/deconv_sample.pdf'
    shell:
        """
        python workflow/scripts/qc/deconv.py -m {input.meta_path} -a {input.ann_path} -s {output.sample} -l {output.leiden}
        """

rule qc_niches_mofa:
    input:
        meta='config/meta.csv'
    params:
        colors_dict=config['colors_areas']
    output:
        'results/qc/niches_mofa.pdf'
    shell:
        """
        python workflow/scripts/qc/niches_mofa.py -m {input.meta} -d '{params.colors_dict}' -p {output}
        """
