
# python scripts/process/annotate.py --output cellbender --resolution 0.75 
# python scripts/process/annotate.py --output cellranger --resolution 0.75

import scanpy as sc

import argparse
from pathlib import Path
from annotate_dicts import annotation_dict

"""
Cluster cells and optionally annotate them.
"""

# Read command line and set args
parser = argparse.ArgumentParser(prog='ra', description='Reanotates atlas')
parser.add_argument("--output", type=str, required=True)
parser.add_argument('--resolution', help='Resoulution for leiden clustering algorithm', required=True)
args = vars(parser.parse_args())

# set up relative paths within the project
current_folder = Path(__file__).parent
output_dir = current_folder / ".." / ".." / "data" / "prc" / "sc"
if args['output'] == "cellbender":
    input_path = current_folder / ".." / ".." / "data" / "prc" / "sc" / "cellbender_integrated.h5ad"
    out_name = "cellbender_annotated.h5ad"
elif args['output'] == "cellranger":
    input_path = current_folder / ".." / ".." / "data" / "prc" / "sc" / "cellranger_integrated.h5ad"
    out_name = "cellranger_annotated.h5ad"
else:
    raise ValueError("output must be either 'cellbender' or 'cellranger'")

resolution = float(args['resolution'])
dictionary = annotation_dict.get(args['output']).get(args['resolution'])
###############################

# Read merged object
adata = sc.read_h5ad(input_path)

# Cluster cells
sc.tl.leiden(adata, resolution=resolution, key_added='leiden')
print(adata.obs.leiden.value_counts())

# Annotate if dict is available
if dictionary is not None:
    
    # Update names
    adata.obs['leiden'] = [dictionary[clust] for clust in adata.obs['leiden']]
    
# Write to file
adata.write(output_dir / out_name)
