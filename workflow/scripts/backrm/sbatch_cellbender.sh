#!/bin/bash
#SBATCH --job-name=example_gpu_job    # Job name
#SBATCH --gres=gpu:1                  # Request one GPU
#SBATCH --mem=16G                     # Request memory (adjust as needed)
#SBATCH --time=240                    # Expected job duration (hh:mm:ss)
#SBATCH --partition=single

# Load any necessary modules
module load system/singularity/3.9.2
module load devel/cuda/11.6

sample_id=$(basename $(dirname $out_path))

# Run your GPU-accelerated task here
singularity exec --nv --bind $(pwd) workflow/envs/cellbender.sif bash -c "mkdir -p data/prc/sn/$sample_id && cd data/prc/sn/$sample_id && cellbender remove-background --input ../../../../$inp_path --output cellbender_matrix.h5 --model $model --cuda --expected-cells $expected_cells  --total-droplets-included $total_droplets_included --fpr $fpr --epochs $epochs --posterior-batch-size $posterior_batch_size --cells-posterior-reg-calc $cells_posterior_reg_calc && rm cellbender_matrix.h5 && rm cellbender_matrix.log && rm cellbender_matrix_cell_barcodes.csv && mv cellbender_matrix_filtered.h5 cellbender_matrix.h5 && mv cellbender_matrix.pdf ../../../../results/backrm/$sample_id.pdf"
