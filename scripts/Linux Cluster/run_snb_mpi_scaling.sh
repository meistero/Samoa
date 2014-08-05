#!/bin/bash
# Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
# Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
# This program is licensed under the GPL, for details see the file LICENSE

cpus=$(lscpu | grep "^CPU(s)" | grep -oE "[0-9]+" | tr "\n" " ")
output_dir=output/Snb_MPI_Scaling_$(date +"%Y-%m-%d_%H-%M-%S")
script_dir=$(dirname "$0")

mkdir -p $output_dir
mkdir -p scripts

echo "CPU(s) detected : "$cpus
echo "Output directory: "$output_dir
echo ""
echo "Compiling..."

partition=snb

scons config=mac.py scenario=darcy openmp=noomp -j4 &
scons config=mac.py scenario=swe openmp=noomp -j4 &

wait %1 %2

echo "Running scenarios..."

class=test
limit=02:00:00
postfix=_noomp

for asagimode in 2
do
	for sections in 8
	do
		for cores in 1 2 4 8 16
		do
			processes=$cores
			threads=1

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

