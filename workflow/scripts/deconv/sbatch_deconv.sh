#!/bin/bash
#SBATCH --job-name=deconv             # Job name
#SBATCH --gres=gpu:1                  # Request one GPU
#SBATCH --mem=32G                     # Request memory (adjust as needed)
#SBATCH --time=240                    # Expected job duration (hh:mm:ss)
#SBATCH --partition=single

# Load any necessary modules
module load system/singularity/3.9.2
module load devel/cuda/11.6

# Run your GPU-accelerated task here
singularity exec --nv --bind $(pwd) workflow/envs/cell2loc.sif python workflow/scripts/deconv/deconv.py -s $slide_path -r $reg_path -n $n_cells_spot -a $d_alpha -e $max_epochs -p $plot_path -o $out_path

