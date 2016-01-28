# Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
# Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
# This program is licensed under the GPL, for details see the file LICENSE


#!/bin/bash

export KMP_AFFINITY="granularity=core,compact,1"

cpus=$(lscpu | grep "^CPU(s)" | grep -oE "[0-9]+" | tr "\n" " ")
output_dir=output/Helium_OpenMP_Scaling_$(date +"%Y-%m-%d_%H-%M-%S")

mkdir -p $output_dir
mkdir -p scripts

export KMP_AFFINITY="granularity=core,compact,1"

echo "CPU(s) detected : "$cpus
echo "Output directory: "$output_dir
echo ""
echo "Compiling..."

scons config=supermuc.py scenario=darcy openmp=noomp -j4 &
scons config=supermuc.py scenario=swe openmp=noomp -j4 &

wait %1 %2

echo "Running scenarios..."

limit=02:00:00

for asagimode in 2
do
	for sections in 8 16 32
	do
		for concurrency in 40
		do
			processes=1
			threads=$concurrency
			nodes=$(( ($processes * $threads - 1) / 40 + 1 ))

			echo "  Running Darcy..."
			./bin/samoa_darcy_nompi -asagihints $asagimode -dmin 16 -dmax 24 -nmax 10 -threads $threads -sections $sections > $output_dir"/darcy_p"$processes"_t"$threads"_s"$sections"_a"$asagimode".log"
			echo "  Done."

			echo "  Running SWE..."
			./bin/samoa_swe_nompi -asagihints $asagimode -dmin 8 -dmax 18 -nmax 100 -threads $threads -sections $sections > $output_dir"/swe_p"$processes"_t"$threads"_s"$sections"_a"$asagimode".log"
			echo "  Done."
		done
	done
done

