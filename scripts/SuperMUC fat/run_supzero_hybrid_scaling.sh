# Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
# Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
# This program is licensed under the GPL, for details see the file LICENSE


#!/bin/bash

cpus=$(lscpu | grep "^CPU(s)" | grep -oE "[0-9]+" | tr "\n" " ")
output_dir=output/Fat_Hybrid_Scaling_$(date +"%Y-%m-%d_%H-%M-%S")
script_dir=$(dirname "$0")

mkdir -p $output_dir
mkdir -p scripts

echo "CPU(s) detected : "$cpus
echo "Output directory: "$output_dir
echo ""
echo "Compiling..."

scons config=supermuc.py scenario=darcy openmp=tasks -j4 &
scons config=supermuc.py scenario=swe openmp=tasks -j4 &
scons config=supermuc.py scenario=darcy openmp=notasks -j4 &
scons config=supermuc.py scenario=swe openmp=notasks -j4 &

wait %1 %2 %3 %4

echo "Running scenarios..."

class=fattest
limit=02:00:00

for asagimode in 2
do
    for postfix in "" _notasks
    do
		for sections in 1 2 3 4
		do
            for cores in 10 20 40 80 120 160
            do
				processes=$(( ($cores - 1) / 10 + 1 ))
				threads=$(( $cores / $processes )) 
				nodes=$(( ($processes * $threads - 1) / 40 + 1 ))

				script="scripts/cache/run_fat"$postfix"_p"$processes"_t"$threads"_s"$sections"_a"$asagimode".sh"
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
			    sed -i 's=-dmin 26=-dmin 26=g' $script
			    sed -i 's=-dmax 29=-dmax 32=g' $script

				llsubmit $script
			done
		done
	done
done

