
from pathlib import Path
import numpy as np
from PIL import Image
import pandas as pd
import scanpy
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.image import imread
import json
import os
import scanpy as sc
from sklearn.decomposition import PCA
from umap import UMAP
from tqdm import tqdm

# not use the GPU (because I temporarily lack GPU memory)
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'

from tensorflow.keras.applications.resnet50 import ResNet50
from tensorflow.keras.applications.resnet50 import preprocess_input as keras_preprocess_input


def Read10X(path, quality = "hires"):  # use quality = "lowres" for low resolution images

    adata = scanpy.read_10x_h5(path / "filtered_feature_bc_matrix.h5")

    tissue_positions_file = (
        path / "spatial/tissue_positions.csv"
        if (path / "spatial/tissue_positions.csv").exists()
        else path / "spatial/tissue_positions_list.csv"
    )

    library_id = "id"  # TODO: placeholder for library id
    adata.uns["spatial"] = {}
    adata.uns["spatial"][library_id] = {}

    files = dict(
        tissue_positions_file=tissue_positions_file,
        scalefactors_json_file=path / "spatial/scalefactors_json.json",
        hires_image=path / "spatial/tissue_hires_image.png",
        lowres_image=path / "spatial/tissue_lowres_image.png",
)

    # check if files exists, continue if images are missing
    adata.uns["spatial"][library_id]["images"] = dict()
    for res in ["hires", "lowres"]:
        try:
            adata.uns["spatial"][library_id]["images"][res] = imread(
                str(files[f"{res}_image"])
            )
        except Exception:
            raise OSError(f"Could not find '{res}_image'")

    # read json scalefactors
    adata.uns["spatial"][library_id]["scalefactors"] = json.loads(
        files["scalefactors_json_file"].read_bytes()
    )

    # read coordinates
    positions = pd.read_csv(files["tissue_positions_file"], header=None)
    positions.columns = [
        "barcode",
        "in_tissue",
        "array_row",
        "array_col",
        "pxl_col_in_fullres",
        "pxl_row_in_fullres",
    ]
    positions.index = positions["barcode"]

    adata.obs = adata.obs.join(positions, how="left")

    adata.obsm["spatial"] = (
        adata.obs[["pxl_row_in_fullres", "pxl_col_in_fullres"]]
        .to_numpy()
        .astype(int)
    )
    adata.obs.drop(
        columns=["barcode", "pxl_row_in_fullres", "pxl_col_in_fullres"],
        inplace=True,
    )

    adata.var_names_make_unique()

    if library_id is None:
        library_id = list(adata.uns["spatial"].keys())[0]

    scale = adata.uns["spatial"][library_id]["scalefactors"][
        "tissue_" + quality + "_scalef"
    ]
    image_coor = adata.obsm["spatial"] * scale

    adata.obs["imagecol"] = image_coor[:, 0]
    adata.obs["imagerow"] = image_coor[:, 1]
    adata.uns["spatial"][library_id]["use_quality"] = quality

    adata.obs["array_row"] = adata.obs["array_row"].astype(int)
    adata.obs["array_col"] = adata.obs["array_col"].astype(int)
    adata.obsm["spatial"] = adata.obsm["spatial"].astype("int64")

    return adata


def get_tiles(adata, out_path, crop_size = 40, target_size = 299):

    library_id = list(adata.uns["spatial"].keys())[0]

    out_path.mkdir(parents=True, exist_ok=True)

    image = adata.uns["spatial"][library_id]["images"][
        adata.uns["spatial"][library_id]["use_quality"]
    ]
    if image.dtype == np.float32 or image.dtype == np.float64:
        image = (image * 255).astype(np.uint8)
    img_pillow = Image.fromarray(image)

    if img_pillow.mode == "RGBA":
        img_pillow = img_pillow.convert("RGB")

    tile_names = []

    for imagerow, imagecol in zip(adata.obs["imagerow"], adata.obs["imagecol"]):
        imagerow_down = imagerow - crop_size / 2
        imagerow_up = imagerow + crop_size / 2
        imagecol_left = imagecol - crop_size / 2
        imagecol_right = imagecol + crop_size / 2
        tile = img_pillow.crop(
            (imagecol_left, imagerow_down, imagecol_right, imagerow_up)
        )
        tile.thumbnail((target_size, target_size), Image.Resampling.LANCZOS)
        tile.resize((target_size, target_size))
        tile_name = str(imagecol) + "-" + str(imagerow) + "-" + str(crop_size)
        out_tile = Path(out_path) / (tile_name + ".jpeg")
        tile_names.append(str(Path(out_path) / (tile_name + ".jpeg")))
        tile.save(out_tile, "JPEG")

    adata.obs["tile_path"] = tile_names


def get_features(adata):

    if not "tile_path" in adata.obs.columns:
        raise ValueError("Please run tiling before extracting features.")

    model = ResNet50(include_top=False, weights="imagenet", pooling="avg")

    features = np.zeros((adata.shape[0], 2048))

    for i, tile_path in enumerate(tqdm(adata.obs.tile_path)):
        tile = Image.open(tile_path)
        tile = np.asarray(tile, dtype="int32")
        tile = tile.astype(np.float32)
        tile = np.stack([tile])
        tile = keras_preprocess_input(tile)
        features[i, :] = model.predict(tile, batch_size=1, verbose=False)

    adata.obsm["image_features"] = features
    pca_coord = PCA(n_components=50).fit_transform(features)
    adata.obsm["image_features_pca"] = pca_coord
    umap_coord = UMAP(n_components=2).fit_transform(pca_coord)
    adata.obsm["image_features_umap"] = umap_coord


current_folder = Path(__file__).parent
image_dir = current_folder / ".." / ".." / "data" / "raw" / "visium"
processed_dir = current_folder / ".." / ".." / "data" / "uscsc_dump" 
tile_out = current_folder / ".." / ".." / "data" / "tiles"
tile_out.mkdir(parents=True, exist_ok=True)
image_features_out = current_folder / ".." / ".." / "data" / "image_features"
image_features_out.mkdir(parents=True, exist_ok=True)
samples = [f for f in os.listdir(image_dir) if not f.startswith(".")]

for sample in samples:
    print(sample)
    adata = Read10X(image_dir / sample / "outs")
    get_tiles(adata, out_path = tile_out / samples[0])
    get_features(adata)
    adata.write(image_features_out / (sample + ".h5ad"))

    annotated_file = processed_dir / ("visium_" + sample + ".h5ad")
    if annotated_file.exists():
        adata_annotated = sc.read_h5ad(annotated_file)
        adata.obs = adata.obs.join(adata_annotated.obs.loc[:, ["leiden"]], how="left")
        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        sns.scatterplot(
            x=adata.obsm["image_features_umap"][:, 0],
            y=adata.obsm["image_features_umap"][:, 1],
            hue=adata.obs["leiden"],
            s=20)

        # save the figure
        fig.savefig(image_features_out / (sample + "img_feature_umap.png"), dpi=300)
