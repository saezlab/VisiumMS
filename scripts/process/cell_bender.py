
from pathlib import Path
import os
import subprocess
import scanpy as sc

#current_folder = globals()['_dh'][0] # for jupyter notebook
current_folder = Path(__file__).parent
raw_location = current_folder / ".." / ".." / "data" / "raw_old_sn"
filtered_location = current_folder / ".." / ".." / "data" / "raw" / "visium"
samples = [sample for sample in os.listdir(raw_location) if not sample.startswith(".")]

for sample in samples:
    print(sample)

    # TODO: Here using as expected cell number the filtered barcode matrix
    filtered_adata = sc.read_10x_h5(filtered_location / sample / "outs" / "filtered_feature_bc_matrix.h5")
    expected_cells = filtered_adata.n_obs
    del filtered_adata

    input_h5 = raw_location / sample / "raw_feature_bc_matrix.h5"
    print(input_h5)
    subprocess.run(["cellbender", "remove-background", 
                    "--input", input_h5, # input file
                    "--output", "cell_bender_matrix.h5",
                    "--model", "full",
                    "--cuda", # cuda enables
                    "--expected-cells", str(expected_cells), # number of expected cells, no default
                    "--total-droplets-included", str(25000), # number of droplets from the rank-ordered UMI plot that will be analyzed, default 25_000
                    "--fpr", str(0.01), # target false positive rate, default 0.01
                    "--epochs", str(150)] # how many epochs to train, default 150
                    )

# call: python scripts/process/cell_bender.py