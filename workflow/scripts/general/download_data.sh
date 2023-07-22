#!/bin/bash

# Check if the correct number of arguments is provided
if [ $# -lt 2 ] || [ $# -gt 3 ]; then
    echo "Usage: $0 <url> <output_path> [-u|--unzip]"
    exit 1
fi

# Parse the arguments
url="$1"
output_path="$2"
unzip_flag=false

if [ "$3" == "-u" ] || [ "$3" == "--unzip" ]; then
    unzip_flag=true
fi

# Function to download the file from the given URL
download_file() {
    wget -q --show-progress -O "$output_path" "$url"
}

# Function to unzip the downloaded file
unzip_file() {
    gzip -d "$output_path"
}

# Main script logic
download_file

if [ "$unzip_flag" = true ]; then
    unzip_file
fi

echo "Download completed!"

