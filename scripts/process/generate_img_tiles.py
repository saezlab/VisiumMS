
import scanpy as sc
import pandas as pd
pd.set_option('mode.chained_assignment', None) # TODO: currently ignoring chained assignments
import numpy as np
import os
from pathlib import Path
import re
import json
from PIL import Image
# note PIL uses cartesian coordinates with (0, 0) in the upper left corner
# https://pillow.readthedocs.io/en/stable/handbook/concepts.html#coordinate-system

# TODO: harcoded configs -> size should be suitable for DL
crop_size = 40
target_size = 299

current_folder = Path(__file__).parent
#current_folder = globals()['_dh'][0]
visium_dir = current_folder / ".." / ".." / "data" / "uscsc_dump"
image_dir = current_folder / ".." / ".." / "data" / "raw" / "visium"

samples = [f for f in os.listdir(visium_dir) if f.startswith("visium")]

for smp in samples:

    base_name = re.sub(r"\.h5ad$", "", smp)
    base_name = re.sub(r"^visium_", "", base_name)
    vis_adata = sc.read_h5ad(visium_dir / smp)

    image_out_smp =  current_folder / ".." / ".." / "data" / "spot_images" / base_name
    image_out_smp.mkdir(parents=True, exist_ok=True)

    print("Processing sample", base_name, "saving in", image_out_smp)

    # check here: https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/spatial
    smp_img_dir = image_dir / base_name / "outs" / "spatial"
    scale_facs = json.load(open(smp_img_dir / "scalefactors_json.json"))

    tissue_pos = pd.read_csv(smp_img_dir / "tissue_positions_list.csv", header=None, names=["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"])
    tissue_pos.set_index("barcode", inplace=True)
    tissue_pos = tissue_pos.loc[vis_adata.obs_names.to_numpy()]

    # generate highres pixel coordinates from fullres pixel coordinates
    coord = tissue_pos[["pxl_row_in_fullres", "pxl_col_in_fullres"]]
    coord["imagerow"] = coord["pxl_row_in_fullres"] * scale_facs["tissue_hires_scalef"]
    coord["imagecol"] = coord["pxl_col_in_fullres"] * scale_facs["tissue_hires_scalef"]

    # from: https://github.com/BiomedicalMachineLearning/stLearn/blob/09c8d7a79979268fe78273fcb25726c5e94c2c6c/stlearn/image_preprocessing/image_tiling.py#L13
    img_key = list(vis_adata.uns["spatial"].keys())
    if len(img_key) == 1:
        img_key = img_key[0]
        print("Single image key found:", img_key)
    else:
        raise ValueError(f"Multiple image keys found in adata.uns['spatial']: {img_key}")
    image = vis_adata.uns["spatial"][img_key]["images"]["hires"] # or lowres?
    if image.dtype == np.float32 or image.dtype == np.float64:
        image = (image * 255).astype(np.uint8)
    img_pillow = Image.fromarray(image)

    if img_pillow.mode == "RGBA":
        img_pillow = img_pillow.convert("RGB")

    for (barcode, imagerow, imagecol) in zip(coord.index, coord["imagerow"], coord["imagecol"]):
        imagerow_down = imagerow - crop_size / 2
        imagerow_up = imagerow + crop_size / 2
        imagecol_left = imagecol - crop_size / 2
        imagecol_right = imagecol + crop_size / 2
        tile = img_pillow.crop(
            (imagecol_left, imagerow_down, imagecol_right, imagerow_up)
        )
        tile.thumbnail((target_size, target_size), Image.Resampling.LANCZOS)
        tile.resize((target_size, target_size))
        tile.save(image_out_smp / (barcode + ".jpeg"), "JPEG")

