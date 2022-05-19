echo '# Training scLVM model #'
/net/data.isilon/ag-saez/bq_pbadia/programs/conda/envs/ms/bin/python scripts/process/get_raw.py -s data/raw/sc/ -a data/prc/sc/annotated.h5ad -o data/prc/sc/raw.h5ad
/net/data.isilon/ag-saez/bq_pbadia/programs/conda/envs/ms/bin/python scripts/process/sclvm.py -r data/prc/sc/raw.h5ad -s data/raw/visium/ -n leiden -d sample_id -g 4096 -o data/prc/visium/
