# Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
# Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
# This program is licensed under the GPL, for details see the file LICENSE


#!/bin/bash

cpus=$(lscpu | grep "^CPU(s)" | grep -oE "[0-9]+" | tr "\n" " ")
output_dir=output/Thin_Production_$(date +"%Y-%m-%d_%H-%M-%S")
script_dir=$(dirname "$0")

mkdir -p $output_dir
mkdir -p scripts

echo "CPU(s) detected : "$cpus
echo "Output directory: "$output_dir
echo ""
echo "Compiling..."

layers=85
#scons config=supermuc.py scenario=darcy flux_solver=upwind layers=$layers openmp=noomp -j4

echo "Running scenarios..."

class=micro
limit=48:00:00

postfix=_noomp_upwind_l$layers 

for asagimode in 2
do
	for sections in 1
	do
		for cores in 64 512
		do
			processes=$cores
			threads=1
			nodes=$(( ($processes * $threads - 1) / 16 + 1 ))
			log_layers=`echo "import numpy; print int(numpy.log2(1 + "$layers"))" | python`

			script="scripts/cache/run_thin"$postfix"_p"$processes"_t"$threads"_s"$sections"_a"$asagimode"_noomp.sh"
			cat "$script_dir/run_supermuc_template.sh" > $script

			sed -i 's=$asagimode='$asagimode'=g' $script
			sed -i 's=$sections='$sections'=g' $script
			sed -i 's=$processes='$processes'=g' $script
			sed -i 's=$threads='$threads'=g' $script
			sed -i 's=$output_dir='$output_dir'=g' $script
			sed -i 's=$nodes='$nodes'=g' $script
			sed -i 's=$limit='$limit'=g' $script
			sed -i 's=$class='$class'=g' $script
			sed -i 's=$postfix='$postfix'=g' $script
		    sed -i 's=-dmin 26=-dmin 14=g' $script
	        sed -i 's=-dmax 40=-dmax 20=g' $script
	        sed -i 's=$add_options=-tsteps -1 -tmax 1.728e8 -max_iter -1 -xmloutput -tout 86.4e4 -courant 0.5 -lbserial =g' $script

			llsubmit $script
		done
	done
done

