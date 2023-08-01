
# custom file to remove clusters

import scanpy as sc
from pathlib import Path

resolution = 0.75
clust_to_remove = ["5", "18", "24"]

current_folder = Path(__file__).parent
input_file = current_folder / ".." / ".." / "data" / "prc" / "sc" / "cellranger_integrated.h5ad"

adata = sc.read_h5ad(input_file)
print(adata)
sc.tl.leiden(adata, resolution=resolution, key_added='leiden')
adata = adata[~adata.obs['leiden'].isin(clust_to_remove), :]
print(adata)

adata.write(current_folder / ".." / ".." / "data" / "prc" / "sc" / "cellranger_integrated_removed_5_18_24.h5ad")
