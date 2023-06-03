## Run RS-FISH spot counting for all tif images in a given folder 
## threshold and sigma parameters were manually determined for each channel prior to this analysis

## binary for bfconvert and channels to process
declare -a channelarray=("1" "2" "3")
bfconvert_bin="/tools/bftools/bfconvert"

## Input output dirs for RS-FISH
input="../images"
single_channels="../single_channel_images"
output="../rsfish_spots"

## Loop over all ome tiff files in dir
for FILE in $input/*.lif
do
    file_base=$(basename $FILE .lif)
    #file_base=(${FILE//./ })
    # file_name=${file_base[0]}
    # outfile_name_test=$file_name""
    for channel in ${channelarray[@]}
    do
        if [ ! -f $output/$file_base".ch_"$channel".csv" ]
        then

            outfile_name_tif=$file_base".ch_"$channel".tif"
            outfile_name_csv=$file_base".ch_"$channel".csv"

            echo $outfile_name_tif

            # bfconvert seems to work fine on converting .lif file types
            $bfconvert_bin -channel $channel -overwrite $FILE $single_channels/$outfile_name_tif

            ## Set the parameters for spot counting based on channel
            if [ "$channel" == "1" ];
            then 
                set_threshold="0.006"
                set_sigma="1.43"
            elif [ "$channel" == "2" ];
            then
                set_threshold="0.012"
                set_sigma="1.4"
            elif [ "$channel" == "3" ];
            then
                set_threshold="0.006"
                set_sigma="1.46"
            fi

            if [ "$channel" != "0" ];
            then
            ## Run RS-FISH for spot detection
            docker run -v $single_channels:/input \
                -v $output:/output \
                rs-fish:2.3.1 /RS-FISH/rs-fish \
                --threshold $set_threshold \
                --sigma $set_sigma \
                --ransac 1 \
                --image=/input/$outfile_name_tif \
                --output=/output/$outfile_name_csv
            fi
        else
            echo $FILE" has already been processed!"
        fi
    done
done
