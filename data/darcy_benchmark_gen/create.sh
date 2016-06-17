#!/bin/bash

for layers in 2 ; do
    for depth in {0..12..2} ; do
        size=$(( 2 ** ($depth / 2 + 1) ))
        echo "depth: $depth, width: $size, layers: $layers"
        python ./create_channel.py -x $size -y $size -z $layers -r 0.2 -n "channel_d"$depth"_l"$layers".nc"
    done
done
