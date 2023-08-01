
from pathlib import Path
import os
from subprocess import Popen, STDOUT
import scanpy as sc

# harcoded configs
recompute = False

#current_folder = globals()['_dh'][0] # for jupyter notebook
current_folder = Path(__file__).parent
input_dir = current_folder / ".." / ".." / "data" / "raw" / "sc"
output_dir = current_folder / ".." / ".." / "data" / "prc" / "sc" / "cell_bender"

# verbose, helpful for debugging
print("recompute: ", recompute)
print("input_dir: ", input_dir)
print("output_dir: ", output_dir)

samples = [sample for sample in os.listdir(input_dir) if not sample.startswith(".")]

for sample in samples:

    output_dir_sample = output_dir / sample
    output_dir_sample.mkdir(parents=True, exist_ok=True)
    print(output_dir_sample)

    if ("cell_bender_matrix.h5" in os.listdir(output_dir_sample)) and (not recompute):
        print(f"Cellbender was run before for sample {sample}")
        continue

    # TODO: Here using as expected cell number the filtered barcode matrix
    filtered_adata = sc.read_10x_h5(input_dir / sample / "filtered_feature_bc_matrix.h5")
    expected_cells = filtered_adata.n_obs
    del filtered_adata

    input_h5 = input_dir / sample / "raw_feature_bc_matrix.h5"
    print(input_h5)

    # create log file
    log_file = output_dir / sample / "cell_bender_custom_log.txt"
    with open(log_file, "w") as f:

        f.write("input_h5: " + str(input_h5) + "\n")

        f.write("output_dir_sample: " + str(output_dir_sample) + "\n")

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
            cwd=output_dir_sample,
            stderr=STDOUT,
            stdout=f
            )
        
        p.wait(timeout=None)    # https://docs.python.org/3/library/subprocess.html#subprocess.Popen.wait
                                # don't spawn another process if p has not finished yet (sequential handling)

# install cellbender: pip install -e .
# call: python scripts/process/cell_bender.py 