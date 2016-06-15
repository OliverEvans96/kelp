#!/bin/bash

# movie.sh
# Create mp4 from png files in img subdirectory of specified directory

dir=$1
name=$2

ffmpeg -framerate 10 -i ${dir}/img/*.png -vcodec mpeg4 -b 800k -r 10 ${dir}/${name}.mp4 -y
