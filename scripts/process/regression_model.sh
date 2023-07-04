#!/bin/sh

# usage
# sbatch --output=/net/data.isilon/ag-saez/bq_pschaefer/VisiumMS/logs/regression_model_cb_gpu.log --partition=gpu --nodes=1 --ntasks-per-node=4 --gres=gpu:1 --time=24:00:00 --mem=16G --export=OUTPUT="cellbender" /net/data.isilon/ag-saez/bq_pschaefer/VisiumMS/scripts/process/regression_model.sh
# sbatch --output=/net/data.isilon/ag-saez/bq_pschaefer/VisiumMS/logs/regression_model_cr_gpu.log --partition=gpu --nodes=1 --ntasks-per-node=4 --gres=gpu:1 --time=24:00:00 --mem=16G --export=OUTPUT="cellranger" /net/data.isilon/ag-saez/bq_pschaefer/VisiumMS/scripts/process/regression_model.sh
# sbatch --output=/net/data.isilon/ag-saez/bq_pschaefer/VisiumMS/logs/regression_model_cb_gpusaez.log --partition=gpusaez --nodes=1 --ntasks-per-node=4 --gres=gpu:1 --time=24:00:00 --mem=16G --export=OUTPUT="cellbender" /net/data.isilon/ag-saez/bq_pschaefer/VisiumMS/scripts/process/regression_model.sh
# sbatch --output=/net/data.isilon/ag-saez/bq_pschaefer/VisiumMS/logs/regression_model_cr_gpusaez.log --partition=gpusaez --nodes=1 --ntasks-per-node=4 --gres=gpu:1 --time=24:00:00 --mem=16G --export=OUTPUT="cellranger" /net/data.isilon/ag-saez/bq_pschaefer/VisiumMS/scripts/process/regression_model.sh

source /home/bq_pschaefer/.bashrc

conda activate ms_env_1

python /net/data.isilon/ag-saez/bq_pschaefer/VisiumMS/scripts/process/deconv.py --output $OUTPUT