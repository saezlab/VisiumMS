#!/bin/bash
#SBATCH --job-name=example_gpu_job    # Job name
#SBATCH --gres=gpu:1                  # Request one GPU
#SBATCH --mem=64G                     # Request memory (adjust as needed)
#SBATCH --time=120                    # Expected job duration (hh:mm:ss)
#SBATCH --partition=single

#######################
# Important variables #
inp_path='data/prc/sn_annotated.h5ad'
out_path='data/prc/vm_regmodel.csv'
plot_path='results/deconv/loss_reg_model.pdf'
max_epochs=250
batch_size=2500
lr=0.002
num_samples=1000
#                     #
#######################

# Load any necessary modules
module load system/singularity/3.9.2
module load devel/cuda/11.6

mkdir -p $(dirname $plot_path)
mkdir -p $(dirname $out_path)
singularity exec --nv --bind $(pwd) workflow/envs/cell2loc.sif python workflow/scripts/deconv/reg_model.py -i $inp_path -e $max_epochs -b $batch_size -l $lr -n $num_samples -p $plot_path -o $out_path
