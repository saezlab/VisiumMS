#!/bin/bash
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH --time=03:00:00
##SBATCH --mem=60000
#SBATCH --job-name="get_props"
#SBATCH --chdir=/net/data.isilon/ag-saez/bq_pbadia/VisiumMS/
#SBATCH --output="logs/get_props.txt"


path_to_slides="/net/data.isilon/ag-saez/bq_pbadia/VisiumMS/data/prc/visium/"

for path_to_slide in ${path_to_slides}*; do
    sample_id=$(basename $path_to_slide);
    echo $sample_id;
    if [[ $sample_id != sclvm ]]
    then
        /net/data.isilon/ag-saez/bq_pbadia/programs/conda/envs/ms/bin/python scripts/process/get_props.py -s $path_to_slide
    fi
done
