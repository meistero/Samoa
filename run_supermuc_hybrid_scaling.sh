# Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
# Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
# This program is licensed under the GPL, for details see the file LICENSE


#!/bin/bash

cpus=$(lscpu | grep "^CPU(s)" | grep -oE "[0-9]+" | tr "\n" " ")
output_dir=output/$(date +"%Y-%m-%d_%H-%M-%S")_Hybrid_Scaling
script_dir=$(dirname $0)

mkdir -p $output_dir
mkdir -p scripts

echo "CPU(s) detected : "$cpus
echo "Output directory: "$output_dir
echo ""
echo "Running scenarios..."

class=test
limit=02:00:00

for asagimode in 2
do
	for sections in 8 16
	do
		for concurrency in 1 2 4 8 16 32 64 128 256
		do
			processes=$(( ($concurrency - 1) / 16 + 1 ))
			threads=$(( $concurrency / $processes )) 
			nodes=$(( ($processes * $threads - 1) / 16 + 1 ))

			script="scripts/run_p"$processes"_t"$threads"_s"$sections"_a"$asagimode".sh"
			cat run_supermuc_template.sh > $script

			sed -i 's=$asagimode='$asagimode'=g' $script
			sed -i 's=$sections='$sections'=g' $script
			sed -i 's=$processes='$processes'=g' $script
			sed -i 's=$threads='$threads'=g' $script
			sed -i 's=$output_dir='$output_dir'=g' $script
			sed -i 's=$nodes='$nodes'=g' $script
			sed -i 's=$limit='$limit'=g' $script
			sed -i 's=$class='$class'=g' $script

			llsubmit $script
		done
	done
done

