! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

PROGRAM gridtest
	USE SFC_traversal

	implicit none

    call cfg%read_from_arguments()

    !init openmp threads
    call omp_set_num_threads(cfg%i_threads)

    !init element transformation data
    call init_transform_data()

    !init MPI
    call init_mpi()

    !run scenario selector
    call sfc_generic()

    !finalize MPI
    call finalize_mpi()

	stop
end PROGRAM gridtest
