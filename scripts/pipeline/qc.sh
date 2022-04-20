echo '## QC ##'
python scripts/process/qc.py -i data/raw/sc/ -o data/prc/sc/
python scripts/plot/qc_metrics.py -i data/prc/sc/ -o figures/
