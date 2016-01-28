# Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
# Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
# This program is licensed under the GPL, for details see the file LICENSE


#!/bin/bash

#@ job_name = samoa
#@ job_type = MPICH
#@ wall_clock_limit = $limit
#@ node = $nodes
#@ total_tasks = $processes
#@ island_count = 1
#@ node_usage = not_shared
#@ class = $class
#@ network.MPI = sn_all,not_shared,us
#@ initialdir = $(home)/Desktop/Samoa
#@ output = $output_dir/run$postfix_p$processes_t$threads_s$sections_a$asagimode.$(jobid).out
#@ error =  $output_dir/run$postfix_p$processes_t$threads_s$sections_a$asagimode.$(jobid).err
#@ energy_policy_tag = samoa_energy_tag
#@ minimize_time_to_solution = yes
#@ queue

. /etc/profile 2>/dev/null
. /etc/profile.d/modules.sh 2>/dev/null

export OMP_NUM_THREADS=$threads
module load mpi_pinning/hybrid_blocked
#export KMP_AFFINITY="granularity=core,compact,1"

echo "  Processes: "$processes
echo "  Threads: "$threads
echo "  Sections: "$sections
echo "  ASAGI mode: "$asagimode

echo "  Running Darcy..."
mpiexec -prepend-rank -n $processes ./bin/samoa_darcy$postfix -dmin 26 -dmax 40 -nmax 10 -asagihints $asagimode -threads $threads -sections $sections > $output_dir"/darcy"$postfix"_p"$processes"_t"$threads"_s"$sections"_a"$asagimode".log"
echo "  Done."

echo "  Running SWE..."
mpiexec -prepend-rank -n $processes ./bin/samoa_swe$postfix -dmin 8 -dmax 29 -nmax 100 -asagihints $asagimode -threads $threads -sections $sections > $output_dir"/swe"$postfix"_p"$processes"_t"$threads"_s"$sections"_a"$asagimode".log"
echo "  Done."

