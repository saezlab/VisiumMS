
# to run the script from the project directory:
# python scripts/process/deep_image_features.py --quality hires
# python scripts/process/deep_image_features.py --quality fullres

# constants
DEFAULT_CROP_SIZES = {
    "lowres": 40,  # not sure whether this makes sense
    "hires": 40   # this is what they use by default in stlearn
}

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
import argparse

# not use the GPU (because I temporarily lack GPU memory)
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'

from tensorflow.keras import backend as K
# check here: https://www.tensorflow.org/api_docs/python/tf/keras/applications/resnet50/ResNet50
from tensorflow.keras.applications.resnet50 import ResNet50
# check here: https://www.tensorflow.org/api_docs/python/tf/keras/applications/resnet50/preprocess_input
from tensorflow.keras.applications.resnet50 import preprocess_input as keras_preprocess_input


# command line arguments for quality to use
parser = argparse.ArgumentParser()
parser.add_argument("--quality", type=str, default="fullres")
parser.add_argument("--crop_size", type=int, default=None)
args = parser.parse_args()

if args.quality not in ["lowres", "hires", "fullres"]:
    raise ValueError("quality must be one of 'lowres', 'hires' or 'fullres'")
if type(args.crop_size) is not int and args.crop_size is not None:
    raise ValueError("crop_size must be an integer or not specified (None)")

QUALITY = args.quality
CROP_SIZE = args.crop_size


def Read10X(spaceranger_path, img_path, library_id="id"):  # use quality = "lowres" for low resolution images

    adata = scanpy.read_10x_h5(spaceranger_path / "filtered_feature_bc_matrix.h5")

    # compute QC metrics
    adata.obs["detected_genes"] = np.sum(adata.X > 0, axis=1)
    adata.obs["total_umis"] = np.sum(adata.X, axis=1)

    tissue_positions_file = (
        spaceranger_path / "spatial/tissue_positions.csv"
        if (spaceranger_path / "spatial/tissue_positions.csv").exists()
        else spaceranger_path / "spatial/tissue_positions_list.csv"
    )

    adata.uns["spatial"] = {}
    adata.uns["spatial"][library_id] = {}

    files = dict(
        tissue_positions_file = tissue_positions_file,
        scalefactors_json_file = spaceranger_path / "spatial/scalefactors_json.json",
        hires_image = spaceranger_path / "spatial/tissue_hires_image.png",
        lowres_image = spaceranger_path / "spatial/tissue_lowres_image.png",
)

    # add lowres and hires images
    adata.uns["spatial"][library_id]["images"] = dict()
    for res in ["lowres", "hires"]:
        try:
            adata.uns["spatial"][library_id]["images"][res] = imread(str(files[f"{res}_image"]))
        except Exception:
            raise OSError(f"Could not find '{res}_image'")
    # add fullres image
    try:
        adata.uns["spatial"][library_id]["images"]["fullres"] = imread(img_path)
    except FileNotFoundError:
        raise OSError(f"Could not find fullres image at {img_path}")
    
    # check type of images and convert to uint8 if necessary
    for key in adata.uns["spatial"][library_id]["images"].keys():
        image = adata.uns["spatial"][library_id]["images"][key]
        if image.dtype == np.float32 or image.dtype == np.float64:
            adata.uns["spatial"][library_id]["images"][key] = (image * 255).astype(np.uint8)

    # read json scalefactors
    adata.uns["spatial"][library_id]["scalefactors"] = json.loads(files["scalefactors_json_file"].read_bytes())

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

    adata.obsm["spatial"] = (adata.obs[["pxl_row_in_fullres", "pxl_col_in_fullres"]].to_numpy().astype(int))

    adata.var_names_make_unique()

    # get scaling for each reslution
    lowres_scale = adata.uns["spatial"][library_id]["scalefactors"]["tissue_lowres_scalef"]
    hires_scale = adata.uns["spatial"][library_id]["scalefactors"]["tissue_hires_scalef"]
    fullres_scale = 1

    # scale image coordiantes for each resolution
    lowres_image_coor = adata.obsm["spatial"] * lowres_scale
    hires_image_coor = adata.obsm["spatial"] * hires_scale
    fullres_image_coor = adata.obsm["spatial"] * fullres_scale

    adata.obs["imagecol_lowress"] = lowres_image_coor[:, 0].astype(int)
    adata.obs["imagerow_lowres"] = lowres_image_coor[:, 1].astype(int)
    adata.obs["imagecol_hires"] = hires_image_coor[:, 0].astype(int)
    adata.obs["imagerow_hires"] = hires_image_coor[:, 1].astype(int)
    adata.obs["imagecol_fullres"] = fullres_image_coor[:, 0].astype(int)
    adata.obs["imagerow_fullres"] = fullres_image_coor[:, 1].astype(int)
    adata.obsm["spatial"] = adata.obsm["spatial"].astype("int64")

    return adata


def get_tiles(adata, out_path, quality, crop_size = 40, target_size = 299):

    library_id = list(adata.uns["spatial"].keys())[0]

    out_path.mkdir(parents=True, exist_ok=True)

    image = adata.uns["spatial"][library_id]["images"][quality]
    if image.dtype == np.float32 or image.dtype == np.float64:
        image = (image * 255).astype(np.uint8)
        raise Warning("Image was converted to uint8")
    img_pillow = Image.fromarray(image)

    if img_pillow.mode == "RGBA":
        img_pillow = img_pillow.convert("RGB")
        raise Warning("Image was converted form RGBA (including alpha channel) to RGB")

    tile_names = []

    for imagerow, imagecol in zip(adata.obs["imagerow_" + quality], adata.obs["imagecol_" + quality]):
        imagerow_down = imagerow - crop_size / 2
        imagerow_up = imagerow + crop_size / 2
        imagecol_left = imagecol - crop_size / 2
        imagecol_right = imagecol + crop_size / 2
        tile = img_pillow.crop(
            (imagecol_left, imagerow_down, imagecol_right, imagerow_up)
        )
        # NOTE: if I understand correctly thumbnail is only needed if the tile is larger than the target size
        tile.thumbnail((target_size, target_size), Image.Resampling.LANCZOS)
        # NOTE: here we resize the tile to the target size which usually means upsampling
        # NOTE: there was a bug in the original code see here: https://github.com/BiomedicalMachineLearning/stLearn/issues/238
        tile = tile.resize((target_size, target_size))
        tile_name = str(imagecol) + "-" + str(imagerow) + "-" + str(crop_size)
        out_tile = Path(out_path) / (tile_name + ".jpeg")
        tile_names.append(str(Path(out_path) / (tile_name + ".jpeg")))
        tile.save(out_tile, "JPEG")

    adata.obs["tile_path"] = tile_names


def get_features(adata):

    if not "tile_path" in adata.obs.columns:
        raise ValueError("Please run tiling before extracting features.")

    model = ResNet50(include_top=False,  # removing fully-connected output layer at the end
                     weights="imagenet", # weights from pre-training on ImageNet
                     pooling="avg"       # global average pool is applied to the output of the last convolutional block
                                         # and thus the output of the model will be a 2D tensor.
                     )
    data_format = K.image_data_format()

    features = np.zeros((adata.shape[0], 2048))

    for i, tile_path in enumerate(tqdm(adata.obs.tile_path)):
        tile = Image.open(tile_path)
        tile = np.asarray(tile, dtype="int32")
        tile = tile.astype(np.float32)
        tile = np.stack([tile])
        # features = encode(tile, model) # this is the original line
        if data_format == "channels_first":
            tile = tile.transpose(0, 3, 1, 2)
        tile = keras_preprocess_input(tile)
        tile = keras_preprocess_input(tile.astype(K.floatx()))
        features[i, :] = model.predict(tile, batch_size=1, verbose=False)

    adata.obsm["image_features"] = features
    pca_coord = PCA(n_components=50).fit_transform(features)
    adata.obsm["image_features_pca"] = pca_coord
    umap_coord = UMAP(n_components=2).fit_transform(pca_coord)
    adata.obsm["image_features_umap"] = umap_coord

current_folder = Path(__file__).parent
image_dir = current_folder / ".." / ".." / "data" / "raw" / "images"
visium_dir = current_folder / ".." / ".." / "data" / "raw" / "vis"
tile_out = current_folder / ".." / ".." / "data" / "prc" / "images" / ("tiles_" + QUALITY)
tile_out.mkdir(parents=True, exist_ok=True)
image_features_out = current_folder / ".." / ".." / "data" / "prc" / "images" / ("deep_image_features_" + QUALITY)
image_features_out.mkdir(parents=True, exist_ok=True)
samples = [f for f in os.listdir(image_dir) if not f.startswith(".")]

for sample in samples:
    print(sample)
    adata = Read10X(spaceranger_path=(visium_dir / sample / "outs"), img_path=(image_dir / sample / (sample+"_pic.tif")))
    
    # verbose output
    for key, val in adata.uns["spatial"]["id"]["images"].items():
        print(key, val.shape)
        print(val.dtype)
    # setting the crop size per tile
    spot_diameter_fullres = int(adata.uns["spatial"]["id"]["scalefactors"]["spot_diameter_fullres"])
    if CROP_SIZE is None and QUALITY == "fullres":
        crop_size = spot_diameter_fullres
    elif CROP_SIZE is None and QUALITY != "fullres":
        crop_size = DEFAULT_CROP_SIZES[QUALITY]
    else:
        crop_size = CROP_SIZE

    # extract tiles and compute features
    get_tiles(adata, out_path = (tile_out / sample), quality=QUALITY, crop_size=crop_size)
    get_features(adata)

    adata.obs["detected_genes"] = np.sum(adata.X > 0, axis=1)
    adata.obs["total_umis"] = np.sum(adata.X, axis=1)
    adata.write(image_features_out / (sample + ".h5ad"))

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    sns.scatterplot(
        x=adata.obsm["image_features_umap"][:, 0],
        y=adata.obsm["image_features_umap"][:, 1],
        hue=adata.obs["detected_genes"],
        s=20)
    plt.suptitle(sample)
    fig.savefig(image_features_out / (sample + "_img_feature_umap_detected_genes.png"), dpi=300)
    plt.close()

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    sns.scatterplot(
        x=adata.obsm["image_features_umap"][:, 0],
        y=adata.obsm["image_features_umap"][:, 1],
        hue=adata.obs["total_umis"],
        s=20)
    plt.suptitle(sample)
    fig.savefig(image_features_out / (sample + "_img_feature_umap_total_umis.png"), dpi=300)
    plt.close()
