#!/bin/bash
#SBATCH -p single
#SBATCH --ntasks 1 --cpus-per-task 24
#SBATCH --time=03:00:00
#SBATCH --mem=60000
#SBATCH --job-name="proc_sc_atlas"
#SBATCH --chdir=/net/data.isilon/ag-saez/bq_pbadia/VisiumMS/
#SBATCH --output="logs/proc_sc_atlas.txt"

echo '# Processing sc atlas #'
bash scripts/pipeline/qc.sh
bash scripts/pipeline/merge.sh
bash scripts/pipeline/integrate.sh
#bash scripts/pipeline/annotate.sh
