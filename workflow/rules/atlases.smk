rule download_absinta2021:
    output:
        annot="data/raw/annot_absinta2021.txt",
        counts="data/raw/count_absinta2021.csv"
    params:
        url_annot=config['atlases']['absinta_2021']['annot'],
        url_count=config['atlases']['absinta_2021']['count'],
    shell:
        """
        bash workflow/scripts/general/download_data.sh '{params.url_annot}' {output.annot} -u
        bash workflow/scripts/general/download_data.sh '{params.url_count}' {output.counts} -u
        """

rule process_absinta2021:
    resources:
        mem_mb=128000,
        partition='cpu-multi'
    input:
        ann_path='data/raw/annot_absinta2021.txt',
        cnt_path='data/raw/count_absinta2021.csv'
    output:
        "data/prc/absinta2021.h5ad"
    shell:
        """
        python workflow/scripts/atlases/process_absinta2021.py -a {input.ann_path} -c {input.cnt_path} -o {output}
        """

rule comp_absinta2021:
    resources:
        mem_mb=64000,
        partition='cpu-multi'
    input:
        lerma='data/prc/sn_annotated.h5ad',
        absinta='data/prc/absinta2021.h5ad'
    output:
        'results/qc/comp_absinta2021.pdf'
    shell:
        """
        python workflow/scripts/atlases/comp_absinta2021.py -c {input.lerma} -a {input.absinta} -p {output}
        """
