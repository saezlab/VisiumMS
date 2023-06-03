# Import local libraries
## Many parts copied from base script by Miguel Ibarra (https://twitter.com/miguelib93)

import time
from pathlib import Path
import argparse
from argparse import ArgumentParser as AP

# Import external libraries
from aicsimageio import AICSImage
from aicsimageio.writers import OmeTiffWriter

import numpy as np
import pandas as pd
from skimage import measure

def get_args():
    # Script description
    description = ""

    # Initialize parser
    parser = AP(description=description, formatter_class=argparse.RawDescriptionHelpFormatter)

    # Parser sections
    tool = parser.add_argument_group(title="Input", description="Required tool input")
    tool.add_argument("-i", "--image", dest="image", action="store", required=True, type=str,
                      help="Path to base image in .ome.tif or .tif format.")
    tool.add_argument("-m", "--mask", dest="mask", action="store", required=True, type=str,
                      help="Path to a labeled mask image in .ome.tif or .tif format.")
    tool.add_argument("-s", "--spots", dest="spots", action="store", required=True, type=str,
                      help="Path to a table with spot positions as produced by RS-fish.")
    tool.add_argument("-c", "--cellmasks", dest="cellmasks", action="store", required=False, type=str,
                      help="Path to labeled cell masks from a segmentation tool.")

    out = parser.add_argument_group(title="Output", description="Tool output")
    out.add_argument("-o", "--output", dest="output", action="store", required=True, type=str,
                     help="Path to output table with accumulated spots in the respective mask labels.")
    out.add_argument("-oi", "--output_image", dest="output_image", action="store", required=False, type=str,
                     help="Path to output tif file with spots.")
    # Parse arguments
    arg = parser.parse_args()

    # Resolve filenames
    arg.mask = Path(arg.mask).resolve()
    arg.image = Path(arg.image).resolve()
    #arg.cellmasks = Path(arg.cellmasks).resolve()
    arg.output = Path(arg.output).resolve()
    arg.output_image = Path(arg.output_image).resolve()

    return arg

def mask_spots(image,mask,spots):
    ## Takes a mask 2d image and a list of spots (x-y positions required) and computes the sum of 
    ## spots in each labeled object
    
    ## Overlap spots with mask
    def f(x):
        return int(round(x))
    f2 = np.vectorize(f)
    
    spot_arr = np.zeros_like(image, dtype= 'uint8')
    spot_list = np.array(spots.values.tolist())
    spot_arr[f2(spot_list[:,1]), f2(spot_list[:,0])] = 1
    
    def sum_intensity(regionmask, intensity_image):
        return np.sum(intensity_image[regionmask])
    
    props = measure.regionprops_table(mask, intensity_image = spot_arr,
                                      properties = ['label', 'centroid', 'orientation', 'eccentricity' , 'area'], 
                                      extra_properties=(sum_intensity,))
    
    props_pd = pd.DataFrame(props)
    return(props_pd)

def main(args):
    def f(x):
        return int(round(x))
    f2 = np.vectorize(f)
    
    # Read mask
    mask = AICSImage(args.mask)
    mask = mask.data[0,0,0,:,:]
    
    image = AICSImage(args.image)
    image = image.data[0,0,0,:,:]

    spots = pd.read_csv(args.spots)
    
    ## Filter spots based on intensity
    spots = spots[["x","y"]]
    
    ## If user also specifies that he wants a tif with spot locations (mostly to make sure spot coordinates match images)
    if args.output_image:
        spot_arr = np.zeros_like(image, dtype= 'uint8')
        spot_list = np.array(spots.values.tolist())
        spot_arr[f2(spot_list[:,1]), f2(spot_list[:,0])] = 255
        OmeTiffWriter.save(spot_arr, args.output_image, dim_order="YX")
    
    # Print the number of labeled objects in the mask file
    print(f"Number of objects in mask: {mask.max()}")
    print(f"Number of spots in table: {len(spots.index)}")
    
    spots_counts = mask_spots(image,mask,spots)
    
    ## If user also gives nuclear or cell masks to count number of cells in annotated region
    if args.cellmasks:
        cellmask = AICSImage(args.cellmasks)
        cellmask = cellmask.data[0,0,0,:,:]
        cell_props = measure.regionprops_table(cellmask, 
                                               properties = ['label', 'centroid', 'orientation', 'eccentricity' , 'area'])
        
        ## Only keep x-y centroids to overlay with mask regions
        cell_props = pd.DataFrame(cell_props)
        cell_props = cell_props[["centroid-1","centroid-0"]]
        cell_counts = mask_spots(image,mask,cell_props)
        spots_counts[["cell_count"]] = cell_counts[["sum_intensity"]]

    ## Write table to csv file
    spots_counts.to_csv(args.output, sep = '\t', index=False)


if __name__ == '__main__':
    # Get arguments
    arg = get_args()

    # Run script and time it
    st = time.time()
    main(arg)
    rt = time.time() - st
    print(f"Script finish in {rt//60:.0f}m {rt%60:.0f}s")