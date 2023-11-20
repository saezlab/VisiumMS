rule fig1_spatial_features:
    input:
        expand("data/prc/vs/{vs_sample}/niches.csv", vs_sample=vs_samples)
    output:
        'results/figures/fig1/spatial_features.pdf'
    params:
        colors_dict=config['colors_areas']
    shell:
        """
        python workflow/scripts/figures/spatial_features.py -d '{params.colors_dict}' -p {output}
        """
