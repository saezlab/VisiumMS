echo '## Merging of samples ##'
python scripts/process/merge.py -i data/prc/sc/ -m data/metadata.csv  -o data/prc/sc/
python scripts/plot/merge_metrics.py -i data/prc/sc/merged.h5ad -o figures/
