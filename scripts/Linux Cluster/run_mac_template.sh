#!/bin/bash
# Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
# Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
# This program is licensed under the GPL, for details see the file LICENSE

#SBATCH -o /home/hpc/pr63so/di56dop/Samoa/$output_dir/run$postfix_p$processes_t$threads_s$sections_a$asagimode.%j.%N.out
#SBATCH -D /home/hpc/pr63so/di56dop/Samoa
#SBATCH -J samoa
#SBATCH --partition=$partition
#SBATCH --get-user-env
#SBATCH --ntasks=$processes
#SBATCH --cpus-per-task=$threads
#SBATCH --export=NONE
##SBATCH --mail-type=end
#SBATCH --mail-user=meistero@in.tum.de
#SBATCH --time=$limit

source /etc/profile.d/modules.sh

export OMP_NUM_THREADS=$threads
module load mpi_pinning/hybrid_blocked

echo "  Processes: "$processes
echo "  Threads: "$threads
echo "  Sections: "$sections
echo "  ASAGI mode: "$asagimode

echo "  Running Darcy..."
mpiexec -prepend-rank -n $processes ./bin/samoa_darcy$postfix -lbsplit -dmin 26 -dmax 40 -nmax 10 -asagihints $asagimode -threads $threads -sections $sections > $output_dir"/darcy"$postfix"_p"$processes"_t"$threads"_s"$sections"_a"$asagimode".log"
echo "  Done."

echo "  Running SWE..."
mpiexec -prepend-rank -n $processes ./bin/samoa_swe$postfix -lbsplit -dmin 8 -dmax 29 -nmax 100 -asagihints $asagimode -threads $threads -sections $sections > $output_dir"/swe"$postfix"_p"$processes"_t"$threads"_s"$sections"_a"$asagimode".log"
echo "  Done."

