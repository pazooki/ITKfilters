#!/bin/bash

# Usage itk_crop input 10 10 10 2 2 2
# will divide input image (20x20x20) into two images of 10x10x10 with no overlap
# There is no checking of valid sizes.
# function itk_crop
# {
input=${1}
base_name="${input%.*}"
suffix="${input##*.}"
sizex=${2}
sizey=${3}
sizez=${4}
geomx=${5}
geomy=${6}
geomz=${7}
indx=0
indy=0
indz=0
printf "input: $input \nsuffix: $suffix \n"
printf "sizes: $sizex , $sizey , $sizez\n"
printf "geoms: $geomx , $geomy , $geomz\n"
linear_index=0
z=0
while [ $z -lt $geomz ]; do
    y=0
    while [ $y -lt $geomy ]; do
        x=0
        while [ $x -lt $geomx ]; do
            linear_index=$(($linear_index+1))
            # echo $linear_index
            echo python extract_image_filter.py $input ${base_name}_${geomx}x${geomy}x${geomz}_tile_${linear_index}.${suffix} \
                $((${indx}+${sizex}*${x})) \
                $((${indy}+${sizey}*${y})) \
                $((${indz}+${sizez}*${z})) \
                $sizex \
                $sizey \
                $sizez
            python extract_image_filter.py $input ${base_name}_${geomx}x${geomy}x${geomz}_tile_${linear_index}.${suffix} \
                $((${indx}+${sizex}*${x})) \
                $((${indy}+${sizey}*${y})) \
                $((${indz}+${sizez}*${z})) \
                $sizex \
                $sizey \
                $sizez
            x=$(($x+1))
        done
        y=$(($y+1))
    done
    z=$(($z+1))
done


