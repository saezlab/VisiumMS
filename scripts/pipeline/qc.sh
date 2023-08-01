echo '## QC ##'
python scripts/process/qc.py --output cellbender
python scripts/process/qc.py --output cellranger
python scripts/plot/qc_metrics.py --output cellbender
python scripts/plot/qc_metrics.py --output cellranger
