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

rule fig2_asc:
    input:
        o='data/prc/ctypes/AS.h5ad',
        a='data/prc/ctypes/AS_ann.csv',
        m='config/meta.csv',
        r='config/c2.cp.reactome.v2023.1.Hs.symbols.gmt',
        c='config/collectri.csv',
    output:
        'results/figures/fig2/asc.pdf'
    shell:
        """
        python workflow/scripts/figures/fig2/asc.py \
        -o '{input.o}' \
        -a '{input.a}' \
        -m '{input.m}' \
        -r '{input.r}' \
        -c '{input.c}' \
        -m '{input.m}' \
        -p {output}
        """