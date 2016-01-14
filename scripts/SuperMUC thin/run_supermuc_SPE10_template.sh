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
#@ output = $output_dir/run$postfix_p$processes_t$threads_s$sections.$(jobid).out
#@ error =  $output_dir/run$postfix_p$processes_t$threads_s$sections.$(jobid).err
#@ energy_policy_tag = samoa_energy_tag
#@ minimize_time_to_solution = yes
#@ queue

. /etc/profile 2>/dev/null
. /etc/profile.d/modules.sh 2>/dev/null

export OMP_NUM_THREADS=$threads

echo "  Processes: "$processes
echo "  Threads: "$threads
echo "  Sections: "$sections
echo "  ASAGI mode: "$asagimode

echo "  Running SPE10..."
mpiexec -prepend-rank -n $processes ./bin/samoa_darcy$postfix -fperm "data/darcy_five_spot/spe_perm_renamed.nc" -fpor "data/darcy_five_spot/spe_por_renamed.nc" -sections $sections -threads $threads -dmin 0 -dmax $dmax -xmloutput -tout 86.4e4 -tmax 172.8e6 -epsilon 1.0e-5 -S_ref_th 1.0e2 -p_ref_th $p_ref_th -courant 1.0 -output_dir $output_dir
echo "  Done."
