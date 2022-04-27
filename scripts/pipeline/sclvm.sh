echo '# Training scLVM model #'
python scripts/process/get_raw.py -s data/raw/sc/ -a data/prc/sc/annotated.h5ad -o data/prc/sc/raw.h5ad
python scripts/process/sclvm.py -r data/prc/sc/raw.h5ad -s data/raw/visium/ -n leiden -d sample_id -g 2000 -o data/prc/visium/
