# Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
# Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
# This program is licensed under the GPL, for details see the file LICENSE


#!/bin/bash

cpus=$(lscpu | grep "^CPU(s)" | grep -oE "[0-9]+" | tr "\n" " ")
output_dir=output/Thin_OpenMP_Scaling_$(date +"%Y-%m-%d_%H-%M-%S")
script_dir=$(dirname "$0")

mkdir -p $output_dir
mkdir -p scripts

echo "CPU(s) detected : "$cpus
echo "Output directory: "$output_dir
echo ""
echo "Compiling..."

scons config=supermuc_intel.py scenario=darcy mpi=no openmp=tasks -j4 &
scons config=supermuc_intel.py scenario=swe mpi=no openmp=tasks -j4 &
scons config=supermuc_intel.py scenario=darcy mpi=no openmp=notasks -j4 &
scons config=supermuc_intel.py scenario=swe mpi=no openmp=notasks -j4 &

wait %1 %2 %3 %4

echo "Running scenarios..."

class=test
limit=02:00:00

for asagimode in 2
do
	for postfix in _nompi_upwind _notasks_nompi_upwind
	do
	    for sections in 1 2 4 8 16 32
	    do
		    for cores in 1 8 16
		    do
			    processes=1
			    threads=$cores
			    nodes=$(( ($processes * $threads - 1) / 16 + 1 ))
				islands=$(( ($nodes - 1) / 512 + 1 ))

				if [ $nodes -le 32 ]; then
		           class=test
		        elif [ $nodes -le 512 ]; then
		           class=general
		        else
		           class=large
		        fi

			    script="scripts/cache/run_thin"$postfix"_p"$processes"_t"$threads"_s"$sections"_a"$asagimode".sh"
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
			    sed -i 's=-dmin 26=-dmin 22=g' $script
			    sed -i 's=-dmax 29=-dmax 29=g' $script

			    llsubmit $script
		    done
        done
	done
done

