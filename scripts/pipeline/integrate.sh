echo '## Integration of samples ##'
python scripts/process/integrate.py --output cellbender
python scripts/process/integrate.py --output cellranger
python scripts/plot/integrate_metrics.py --output cellbender
python scripts/plot/integrate_metrics.py --output cellranger
