echo '## Merging of samples ##'
python scripts/process/merge.py --output cellbender
python scripts/process/merge.py --output cellranger
python scripts/plot/merge_metrics.py --output cellbender
python scripts/plot/merge_metrics.py --output cellranger
