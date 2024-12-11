# VisiumMS
Scripts to reproduce results from "Spatial cell type mapping of multiple sclerosis lesions".

This pipeline uses `Snakemake` to ensure reproducibility

## Data
Processed data in `h5ad` format can be downloaded from the Human Cell Atlas [portal](https://explore.data.humancellatlas.org/projects/4c8e9d75-d85a-47de-9598-06549cf44b91), while raw data is available on [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279183).

## Installation

Clone repo:
```
git clone git@github.com:saezlab/VisiumMS.git
cd VisiumMS
```

Install `mamba` (this might take a while) to install packages faster:
```
conda install -n base -c conda-forge mamba
```

Then create a new enviroment specific for `Snakemake`:
```
mamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
```

To run specific outputs you can run:
```
snakemake --profile config/slurm/  path/to/output/file 
```

To run the complete processing of the data (this can take a while), run:
```
snakemake --profile config/slurm/
```

## Reference
Lerma-Martin, C., Badia-i-Mompel, P., Ramirez Flores, R.O. et al. Cell type mapping reveals tissue niches and interactions in subcortical multiple sclerosis lesions. Nat Neurosci 27, 2354–2365 (2024). https://doi.org/10.1038/s41593-024-01796-z
