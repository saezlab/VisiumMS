#!/bin/bash
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH --time=03:00:00
#SBATCH --chdir=/net/data.isilon/ag-saez/bq_pbadia/VisiumMS/
#SBATCH --exclude="cln146g"
#SBATCH -J regr_m
#SBATCH -o "logs/regression_model.out"

echo '# Generating raw object from samples #'
/net/data.isilon/ag-saez/bq_pbadia/programs/conda/envs/ms/bin/python scripts/process/get_raw.py -s data/raw/sc/ -a data/prc/sc/annotated.h5ad -o data/prc/sc/raw.h5ad
echo '# Training regression model #'
/net/data.isilon/ag-saez/bq_pbadia/programs/conda/envs/ms/bin/python scripts/process/regression_model.py -r data/prc/sc/raw.h5ad -n leiden -s sample_id -p 1 -o data/prc/visium/

