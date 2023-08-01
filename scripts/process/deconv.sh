#!/bin/bash

source /home/bq_pschaefer/.bashrc

conda activate ms_env_1

python /net/data.isilon/ag-saez/bq_pschaefer/VisiumMS/scripts/process/deconv.py --output $OUTPUT --model $REG_MODEL --sample $SAMPLE --n_cells_spot $N_CELLS_SPOT --d_alpha $D_ALPHA --recompute $RECOMPUTE
