#!/bin/bash
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH --time=06:00:00
#SBATCH --chdir=/net/data.isilon/ag-saez/bq_pbadia/VisiumMS/
#SBATCH --exclude=cln146g

# Deconvolute a single slide
sample_id=$(basename $path_to_slide)
/net/data.isilon/ag-saez/bq_pbadia/programs/conda/envs/ms/bin/python scripts/process/deconv.py -r $path_to_regres -s $path_to_slide -o $path_to_output/$sample_id
