# Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
# Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
# This program is licensed under the GPL, for details see the file LICENSE


#!/bin/bash

cpus=$(lscpu | grep "^CPU(s)" | grep -oE "[0-9]+" | tr "\n" " ")
output_dir=output/Thin_Layers_$(date +"%Y-%m-%d_%H-%M-%S")
script_dir=$(dirname "$0")

mkdir -p scripts

echo "CPU(s) detected : "$cpus
echo "Output directory: "$output_dir
echo ""
echo "Compiling..."

limit=00:30:00

#for layers in 0 1 4 8 16 32 64 85
#do
#    scons config=supermuc_intel.py scenario=darcy flux_solver=upwind layers=$layers openmp=noomp -j4 &
#done

#wait %1 %2 %3 %4 %5 %6 %7 %8

if [ $? -ne 0 ]; then
    exit
fi

mkdir -p $output_dir

echo "Running scenarios..."

for layers in 0 1 4 8 16 32 64 85
do
    if [ $layers -ge 1 ]; then
        postfix=_noomp_upwind_l$layers
    else
        postfix=_noomp_upwind
    fi

	for sections in 8
	do
		for cores in 16
		do
            asagimode=2
			processes=$cores
			threads=1
			nodes=$(( ($processes * $threads - 1) / 16 + 1 ))
			islands=$(( ($nodes - 1) / 512 + 1 ))
			log_cores=`echo "import numpy; print int(numpy.log2("$cores"))" | python`
			log_layers=`echo "import numpy; print int(numpy.log2(1 + "$layers"))" | python`

			if [ $nodes -le 32 ]; then
               class=test
            elif [ $nodes -le 512 ]; then
               class=general
            else
               class=large
            fi

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
			sed -i 's=$islands='$islands'=g' $script
			sed -i 's=$postfix='$postfix'=g' $script
		    sed -i 's=-dmin 26=-dmin '$((22 - log_layers + log_cores))'=g' $script
	        sed -i 's=-dmax 40=-dmax '$((23 + log_cores))'=g' $script

			llsubmit $script
		done
	done
done

