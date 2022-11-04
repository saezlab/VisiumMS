# VisiumMS
Scripts to reproduce results from "Spatial cell type mapping of multiple sclerosis lesions".

## Structure

Inside `scripts` there are different subdirectories:

- `figures/`: Python and R scripts to reproduce the mansucript's figures
- `launchers/`: Scripts to send jobs to the cluster using SLURM
- `pipeline/`: Scripts to run processing scripts in `process` sequentaly
- `plot/`: General plotting scripts
- `process/`: Scripts to process the raw snRNA-seq and Visium data
- `smFISH/`: Scripts to reproduce the smFISH imaging quantification

## Reference
Lerma-Martin, C., Badia-i-Mompel, P. et al. **Spatial cell type mapping of multiple sclerosis lesions**. bioRxiv 2022.11.03.514906 (2022) doi:10.1101/2022.11.03.514906.
