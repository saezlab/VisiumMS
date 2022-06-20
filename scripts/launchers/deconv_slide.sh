#!/bin/bash
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH --time=03:00:00
#SBATCH --chdir=/net/data.isilon/ag-saez/bq_pbadia/VisiumMS/


# Deconvolute a single slide
/net/data.isilon/ag-saez/bq_pbadia/programs/conda/envs/ms/bin/python scripts/process/deconv.py -r $path_to_regres -s $path_to_slide -o $path_to_output
