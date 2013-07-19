# Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
# Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
# This program is licensed under the GPL, for details see the file LICENSE


#!/bin/bash

cpus=$(lscpu | grep "^CPU(s)" | grep -oE "[0-9]+" | tr "\n" " ")
output_dir=output/$(date +"%Y-%m-%d_%H-%M-%S")_OpenMP_Scaling
script_dir=$(dirname $0)

mkdir -p $output_dir
mkdir -p scripts

export KMP_AFFINITY="granularity=core,compact,1"

echo "CPU(s) detected : "$cpus
echo "Output directory: "$output_dir
echo ""
echo "Running scenarios..."

limit=02:00:00

for asagimode in 2
do
	for sections in 8
	do
		for processes in 1
		do
			for threads in 8 16 24 32 48
			do
				echo "  Running Darcy..."
				./bin/samoa_darcy_nompi -asagihints $asagimode -dmin 26 -dmax 40 -tsteps 10 -threads $threads -sections $sections > $output_dir"/darcy_p"$processes"_t"$threads"_s"$sections"_a"$asagimode".log"
				echo "  Done."

				echo "  Running SWE..."
				./bin/samoa_swe_nompi -asagihints $asagimode -dmin 8 -dmax 30 -tsteps 100 -threads $threads -sections $sections > $output_dir"/swe_p"$processes"_t"$threads"_s"$sections"_a"$asagimode".log"
				echo "  Done."
			done
		done
	done
done

