# VisiumMS
Scripts to reproduce results from "Spatial cell type mapping of multiple sclerosis lesions".

This pipeline uses `Snakemake` to ensure reproducibility

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
Lerma-Martin, C., Badia-i-Mompel, P. et al. **Spatial cell type mapping of multiple sclerosis lesions**. bioRxiv 2022.11.03.514906 (2022) doi:10.1101/2022.11.03.514906.
