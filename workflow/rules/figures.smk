rule fig1_spatial_features:
    input:
        expand("data/prc/vs/{vs_sample}/niches.csv", vs_sample=vs_samples)
    output:
        'results/figures/fig1/spatial_features.pdf'
    params:
        colors_dict=config['colors_areas']
    shell:
        """
        python workflow/scripts/figures/fig1/spatial_features.py -d '{params.colors_dict}' -p {output}
        """

rule fig1_niches_markers:
    input:
        m='config/meta.csv',
        n=expand("data/prc/vs/{vs_sample}/niches.csv", vs_sample=vs_samples)
    output:
        'results/figures/fig1/niches_markers.pdf'
    shell:
        """
        python workflow/scripts/figures/fig1/niches_markers.py -m '{input.m}' -p {output}
        """
