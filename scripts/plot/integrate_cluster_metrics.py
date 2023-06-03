
# python scripts/plot/integrate_metrics.py --output cellbender
# python scripts/plot/integrate_metrics.py --output cellranger

import scanpy as sc
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import argparse
from pathlib import Path
import os

"""
Script to plot different metrics after integration.
"""

# add command line flag arguments to specify either "cellbender" or "cellranger" output
parser = argparse.ArgumentParser()
parser.add_argument("--output", type=str, required=True)
args = parser.parse_args()

# set up relative paths within the project
current_folder = Path(__file__).parent
output_dir = current_folder / ".." / ".." / "data" / "prc" / "sc"
if args.output == "cellbender":
    input_path = current_folder / ".." / ".." / "data" / "prc" / "sc" / "cellbender_annotated.h5ad"
    output_dir = current_folder / ".." / ".." / "out" / "cellbender_integrated"
elif args.output == "cellranger":
    input_path = current_folder / ".." / ".." / "data" / "prc" / "sc" / "cellranger_annotated.h5ad"
    output_dir = current_folder / ".." / ".." / "out" / "cellranger_integrated"
else:
    raise ValueError("output must be either 'cellbender' or 'cellranger'")
output_dir.mkdir(parents=True, exist_ok=True)

# verbose
print("input_path: ", input_path)
print("output_dir: ", output_dir)

version = "0"

###############################

# Read merged object
adata = sc.read_h5ad(input_path)

