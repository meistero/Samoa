# Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
# Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
# This program is licensed under the GPL, for details see the file LICENSE


#!/bin/bash

cpus=$(lscpu | grep "^CPU(s)" | grep -oE "[0-9]+" | tr "\n" " ")
output_dir=output/Thin_SPE10_$(date +"%Y-%m-%d_%H-%M-%S")
script_dir=$(dirname "$0")

mkdir -p scripts

echo "CPU(s) detected : "$cpus
echo "Output directory: "$output_dir
echo ""
echo "Compiling..."

limit=48:00:00
layers=85
postfix=_spe10

if [ $layers = 0 ] ; then
    p_ref_th=3.0e-9
else
    p_ref_th=2.0e-8
fi

. /etc/profile 2>/dev/null
. /etc/profile.d/modules.sh 2>/dev/null

module switch mpi.intel mpi.ibm
module unload gcc
module load gcc/4.7

scons config=supermuc_ibm.py scenario=darcy flux_solver=upwind layers=$layers exe=samoa_darcy_spe10 -j4 &
wait

if [ $? -ne 0 ]; then
    exit
fi

mkdir -p $output_dir

echo "Running scenarios..."

sections=1

for dmax in 14
do
	for cores in 128
	do
	    processes=$(( ($cores - 1) / 16 + 1 ))
	    threads=$(( $cores / $processes )) 
	    nodes=$(( ($processes * $threads - 1) / 16 + 1 ))
		islands=$(( ($nodes - 1) / 512 + 1 ))

		if [ $nodes -le 32 ]; then
           class=test
        elif [ $nodes -le 512 ]; then
           class=general
        else
           class=large
        fi

		script="scripts/cache/run_thin"$postfix"_p"$processes"_t"$threads"_s"$sections"_noomp.sh"
		cat "$script_dir/run_supermuc_SPE10_template.sh" > $script

		sed -i 's=$sections='$sections'=g' $script
		sed -i 's=$processes='$processes'=g' $script
		sed -i 's=$threads='$threads'=g' $script
		sed -i 's=$output_dir='$output_dir'=g' $script
		sed -i 's=$nodes='$nodes'=g' $script
		sed -i 's=$limit='$limit'=g' $script
		sed -i 's=$class='$class'=g' $script
        sed -i 's=$islands='$islands'=g' $script
		sed -i 's=$postfix='$postfix'=g' $script
        sed -i 's=$dmax='$dmax'=g' $script
        sed -i 's=$p_ref_th='$p_ref_th'=g' $script

        #cat $script
		llsubmit $script
	done
done

