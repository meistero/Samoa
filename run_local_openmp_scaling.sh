# Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
# Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
# This program is licensed under the GPL, for details see the file LICENSE


#!/bin/bash

export KMP_AFFINITY="granularity=core,compact,1"

cpus=$(lscpu | grep "^CPU(s)" | grep -oE "[0-9]+" | tr "\n" " ")
output_dir=output/$(date +"%Y-%m-%d_%H-%M-%S")_OpenMP_Scaling
script_dir=$(dirname $0)

mkdir -p $output_dir
mkdir -p scripts

export KMP_AFFINITY="granularity=core,compact,1"

echo "CPU(s) detected : "$cpus
echo "Output directory: "$output_dir
echo ""

echo "Compiling..."

make SCENARIO=DARCY MPI=NO OPENMP=NO -j
make SCENARIO=DARCY MPI=NO OPENMP=YES -j
make SCENARIO=DARCY MPI=NO OPENMP=TASKS -j

make SCENARIO=SWE MPI=NO OPENMP=NO -j
make SCENARIO=SWE MPI=NO OPENMP=YES -j
make SCENARIO=SWE MPI=NO OPENMP=TASKS -j

echo "Running scenarios..."

for omp in _noomp _notasks ""
do
	for threads in 1 2 3 4
	do
		for asagimode in 2
		do
			for sections in 1 2
			do
				echo "  Running Darcy..."
				"./bin/samoa_darcy"$omp"_nompi" -asagihints $asagimode -dmin 18 -dmax 18 -tsteps 10 -threads $threads -sections $sections > $output_dir"/darcy"$omp"_t"$threads"_s"$sections"_a"$asagimode".log"
				echo "  Done."

				echo "  Running SWE..."
				"./bin/samoa_swe"$omp"_nompi" -asagihints $asagimode -dmin 18 -dmax 18 -tsteps 20 -threads $threads -sections $sections > $output_dir"/swe"$omp"_t"$threads"_s"$sections"_a"$asagimode".log"
				echo "  Done."
			done
		done
	done
done

