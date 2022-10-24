## Run the assign spots script for each sample:
## 1) Assign and count spots to the manual masks
## 3) Assign and count nuclei to manual masks

## 
script_dir=../scripts

img_dir=../data/single_channel_images
spots_dir=../data/rsfish_spots
mask_dir=../data/masks

out_dir=../data/result_tables
out_dir_img=../data/spot_images

## Make a list of all channel names with spots
declare -a channelarray=("1" "2" "3")

for FILE in $mask_dir/*.tif
do
    file_base=$(basename $FILE .ch_0mask.tif)
    for channel in ${channelarray[@]}
    do
        spots_file=$spots_dir"/"$file_base".ch_"$channel".csv"
        echo $file_base $file_base".ch_"$channel".csv"

        python $script_dir"/assign_spots.py" \
        -i $img_dir"/"$file_base".ch_"$channel".tif" \
        -m $FILE \
        -s $spots_file \
        -o $out_dir"/"$file_base".results.ch_"$channel".csv" \
        -oi $out_dir_img"/"$file_base".results.ch_"$channel".tif"
    done
done