
# usage
# python scripts/process/cell_bender_c2l_submit.py --model cond --recompute False
# python scripts/process/cell_bender_c2l_submit.py --model all --recompute False

import os
from pathlib import Path
import argparse

# get cmd line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--model", type=str, required=True, help="['all', 'cond']")
parser.add_argument("--n_cells_spot", type=int, required=False, default=5, help="prior for number of cells per spot")
parser.add_argument("--d_alpha", type=int, required=False, default=20, help="prior for heterogeneity?")
parser.add_argument("--recompute", type=str, required=False, default="False", help="prior for heterogeneity?")
args = parser.parse_args()

# check the arguments
if args.model not in ["all", "cond"]:
    raise ValueError("Model must be either 'all' or 'cond'")

reg_model = args.model
n_cells_spot = args.n_cells_spot
d_alpha = args.d_alpha
recompute = args.recompute in ["True", "true"]

current_folder = Path(__file__).parent
visium_dir = current_folder / ".." / ".." / "data" / "raw" / "visium"
samples = [f for f in os.listdir(visium_dir) if not f.startswith('.')]
print(samples)
script = current_folder / "cell_bender_c2l.sh"
script = str(script.resolve())
print(script)
log_dir = current_folder / ".." / ".." / "logs"
log_dir.mkdir(parents=True, exist_ok=True)
log_dir = log_dir.resolve()

# sbatch configs
MEM_GB = "16"
TIME = "24:00:00"
N_CPU = "4"

# submit GPU batch jobs for each sample
# how to deal with the right conda env?
for sample in samples:
    print("submitting job for sample: ", sample)
    log_file = log_dir / f"$(date +'%Y-%m-%d-%H-%M')-{sample}-%j.log"
    command = f"sbatch --output={log_file} --partition=gpu --nodes=1 --ntasks-per-node={N_CPU} --gres=gpu:1 --time={TIME} --mem={MEM_GB}G --export=REG_MODEL={reg_model},SAMPLE={sample},N_CELLS_SPOT={n_cells_spot},D_ALPHA={d_alpha},RECOMPUTE={recompute} {script}"
    print(command)
    os.system(command)
