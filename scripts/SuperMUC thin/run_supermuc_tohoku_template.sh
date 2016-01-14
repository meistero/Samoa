# Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
# Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
# This program is licensed under the GPL, for details see the file LICENSE


#!/bin/bash

#@ job_name = samoa
#@ job_type = parallel
#@ wall_clock_limit = $limit
#@ node = $nodes
#@ total_tasks = $processes
#@ island_count = $islands
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

module switch mpi.intel mpi.ibm
module switch gcc gcc/4.7

export OMP_NUM_THREADS=$threads

if [ $threads -ge 2 ] ; then
    export MP_SINGLE_THREAD=no
    export MP_TASK_AFFINITY=core:$threads
fi

echo "  Processes: "$processes
echo "  Threads: "$threads
echo "  Sections: "$sections

echo "  Running Tohoku..."
mpiexec -n $processes ./bin/samoa_swe$postfix -lbtime -sections $sections -threads $threads -courant 0.95 -tout 20.0 -dmin 0 -dmax $dmax -tmax 10.8e3 -fdispl "data/tohoku_static/displ.nc" -fbath "data/tohoku_static/bath.nc" -stestpoints "545735.266126 62716.4740303,935356.566012 -817289.628677,1058466.21575 765077.767857" -output_dir $output_dir
echo "  Done."
