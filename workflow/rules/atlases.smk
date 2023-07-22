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
