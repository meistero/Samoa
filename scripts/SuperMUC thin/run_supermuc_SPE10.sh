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

class=general
limit=48:00:00

layers=85
postfix=_noomp_upwind_l$layers
sections=8

#scons config=supermuc.py scenario=darcy flux_solver=upwind layers=$layers openmp=noomp -j4 &
#wait

echo "Running scenarios..."

if [ $layers = 0 ] ; then
    p_ref_th=3.0e-9
else
    p_ref_th=2.0e-8
fi

for dmax in 14 16
do
	for cores in 1024
	do
		processes=$cores
		threads=1
		nodes=$(( ($processes * $threads - 1) / 16 + 1 ))
    
		if [ $nodes -le 32 ]; then
           class=test
        else
           class=general
        fi

		script="scripts/cache/run_thin"$postfix"_p"$processes"_t"$threads"_s"$sections"_a"$asagimode"_noomp.sh"
		cat "$script_dir/run_supermuc_SPE10_template.sh" > $script

		sed -i 's=$sections='$sections'=g' $script
		sed -i 's=$processes='$processes'=g' $script
		sed -i 's=$threads='$threads'=g' $script
		sed -i 's=$output_dir='$output_dir'=g' $script
		sed -i 's=$nodes='$nodes'=g' $script
		sed -i 's=$limit='$limit'=g' $script
		sed -i 's=$class='$class'=g' $script
		sed -i 's=$postfix='$postfix'=g' $script
        sed -i 's=$dmax='$dmax'=g' $script
        sed -i 's=$p_ref_th='$p_ref_th'=g' $script

        cat $script
		#llsubmit $script
	done
done

