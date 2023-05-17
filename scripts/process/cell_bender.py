
from pathlib import Path
import os
from subprocess import Popen, STDOUT, PIPE
import scanpy as sc
import argparse

# TODO: temporary command line argument
parser = argparse.ArgumentParser()
parser.add_argument("recompute", help='whether cellbender should be rerun if output is already present', default=True)
args = parser.parse_args()
recompute = args.recompute in [True, "True", "true"]
print(f"{recompute=}")

#current_folder = globals()['_dh'][0] # for jupyter notebook
current_folder = Path(__file__).parent
raw_location = current_folder / ".." / ".." / "data" / "raw_old_sn"
filtered_location = current_folder / ".." / ".." / "data" / "raw" / "visium"
samples = [sample for sample in os.listdir(raw_location) if not sample.startswith(".")]

for sample in samples:

    if ("cell_bender_matrix.h5" in os.listdir(raw_location / sample)) and (not recompute):
        print(f"Cellbender was run before for sample {sample}")
        continue

    print(sample)

    # TODO: Here using as expected cell number the filtered barcode matrix
    filtered_adata = sc.read_10x_h5(filtered_location / sample / "outs" / "filtered_feature_bc_matrix.h5")
    expected_cells = filtered_adata.n_obs
    del filtered_adata

    input_h5 = raw_location / sample / "raw_feature_bc_matrix.h5"
    print(input_h5)

    # create log file
    log_file = raw_location / sample / "cell_bender_log.txt"
    with open(log_file, "w") as f:

        p = Popen(["cellbender", "remove-background", 
            "--input", input_h5, # input file
            "--output", "cell_bender_matrix.h5",
            "--model", "full",
            "--cuda", # cuda enables
            "--expected-cells", str(expected_cells), # number of expected cells, no default
            "--total-droplets-included", str(25000), # number of droplets from the rank-ordered UMI plot that will be analyzed, default 25_000
            "--fpr", str(0.01), # target false positive rate, default 0.01
            "--epochs", str(150), # how many epochs to train, default 150
            "--posterior-batch-size",  str(5), # batch size for creating the posterior, default 20 (important if running of of GPU RAM)
            "--cells-posterior-reg-calc", str(50)], # number of cells used to estimate posterior regularization lambda, default 100
            cwd=raw_location / sample,
            stderr=STDOUT,
            stdout=f
            )
        
        p.wait(timeout=None)    # https://docs.python.org/3/library/subprocess.html#subprocess.Popen.wait
                                # don't spawn another process if p has not finished yet (sequential handling)

# install cellbender: pip install -e .
# call: python scripts/process/cell_bender.py