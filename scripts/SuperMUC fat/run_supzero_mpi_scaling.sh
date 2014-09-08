# Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
# Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
# This program is licensed under the GPL, for details see the file LICENSE


#!/bin/bash

cpus=$(lscpu | grep "^CPU(s)" | grep -oE "[0-9]+" | tr "\n" " ")
output_dir=output/Fat_MPI_Scaling_$(date +"%Y-%m-%d_%H-%M-%S")
script_dir=$(dirname "$0")

mkdir -p $output_dir
mkdir -p scripts

echo "CPU(s) detected : "$cpus
echo "Output directory: "$output_dir
echo ""
echo "Compiling..."

scons config=supermuc.py scenario=darcy openmp=noomp -j4 &
scons config=supermuc.py scenario=swe openmp=noomp -j4 &

wait %1 %2

echo "Running scenarios..."

class=fattest
limit=02:00:00
postfix=_noomp

for asagimode in 2
do
	for sections in 8 16 32
	do
		for concurrency in 1 2 5 10 20 40 80 120 160
		do
			processes=$concurrency
			threads=1
			nodes=$(( ($processes * $threads - 1) / 40 + 1 ))

			script="scripts/cache/run_fat"$postfix"_p"$processes"_t"$threads"_s"$sections"_a"$asagimode"_noomp.sh"
			cat "$script_dir/run_supzero_template.sh" > $script

			sed -i 's=$asagimode='$asagimode'=g' $script
			sed -i 's=$sections='$sections'=g' $script
			sed -i 's=$processes='$processes'=g' $script
			sed -i 's=$threads='$threads'=g' $script
			sed -i 's=$output_dir='$output_dir'=g' $script
			sed -i 's=$nodes='$nodes'=g' $script
			sed -i 's=$limit='$limit'=g' $script
			sed -i 's=$class='$class'=g' $script
			sed -i 's=$postfix='$postfix'=g' $script

			llsubmit $script
		done
	done
done

