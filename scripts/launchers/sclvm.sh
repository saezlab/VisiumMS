#!/bin/bash
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH --time=03:00:00
#SBATCH --mem=60000
#SBATCH --job-name="sclvm_train"
#SBATCH --chdir=/net/data.isilon/ag-saez/bq_pbadia/VisiumMS/
#SBATCH --output="logs/sclvm_train.txt"

bash scripts/pipeline/sclvm.sh

