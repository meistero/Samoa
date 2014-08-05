#!/bin/bash
# Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
# Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
# This program is licensed under the GPL, for details see the file LICENSE

cpus=$(lscpu | grep "^CPU(s)" | grep -oE "[0-9]+" | tr "\n" " ")
output_dir=output/Wsm_Hybrid_Scaling_$(date +"%Y-%m-%d_%H-%M-%S")
script_dir=$(dirname "$0")

mkdir -p $output_dir
mkdir -p scripts

echo "CPU(s) detected : "$cpus
echo "Output directory: "$output_dir
echo ""
echo "Compiling..."

partition=wsm

scons config=wsm.py scenario=darcy -j4 &
scons config=wsm.py scenario=swe -j4 &

wait %1 %2

echo "Running scenarios..."

class=test
limit=02:00:00
postfix=

for asagimode in 2
do
	for sections in 1
	do
		for cores in 8 16 32
		do
			processes=$(( ($cores - 1) / 32 + 1 ))
			threads=$(( $cores / $processes )) 

			script="scripts/cache/run_mac"$postfix"_p"$processes"_t"$threads"_s"$sections"_a"$asagimode".sh"
			cat "$script_dir/run_mac_template.sh" > $script

            sed -i 's=$partition='$partition'=g' $script
			sed -i 's=$asagimode='$asagimode'=g' $script
			sed -i 's=$sections='$sections'=g' $script
			sed -i 's=$processes='$processes'=g' $script
			sed -i 's=$threads='$threads'=g' $script
			sed -i 's=$output_dir='$output_dir'=g' $script
			sed -i 's=$nodes='$nodes'=g' $script
			sed -i 's=$limit='$limit'=g' $script
			sed -i 's=$class='$class'=g' $script
			sed -i 's=$postfix='$postfix'=g' $script
		    sed -i 's=-dmin 26=-dmin 23=g' $script
	        sed -i 's=-dmax 29=-dmax 29=g' $script

			sbatch $script
		done
	done
done

