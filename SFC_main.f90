! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

PROGRAM gridtest
	USE SFC_traversal
	USE M_kracken

	implicit none

!---------- local declarations
	logical (kind = GRID_SL), dimension(2,3,5)			:: l_boundary_data		! 2 boundary bits of 3 edges, 4 initial triangles
	real (kind = GRID_SR), dimension(2,3,5)				:: r_coords2		! 2 coordinates of 3 nodes, 4 initial triangles

	type(triangle_tree), pointer						:: p_triangle_tree  ! pointer to the triangle elements in the tree

	! Timing
	real (kind = GRID_SR)						        :: r_t1, r_t2

	! Element statistics
	integer (kind = GRID_SI)					        :: i_initial_elements, i_initial_fine_elements

	! command line arguments
	logical					                            :: l_help, l_version, mpi_flag
	integer          									:: i, i_error, i_threads, mpi_tag_upper_bound, mpi_prov_thread_support

	character(256)                                      :: arguments

!--------------------- PROGRAM START ---------------------------------------------------------------------------------------

	!define default command arguments and default values for all scenarios
    write(arguments, '(A, I0)') "-v .false. --version .false. -h .false. --help .false. -asagihints 2 -noprint .false. -sections 4 -threads ", omp_get_max_threads()

	!define additional command arguments and default values depending on the choice of the scenario
#	if defined(_DARCY)
        write(arguments, '(A, A)') trim(arguments), "  -dmin 1 -dmax 14 -tsteps -1 -tmax 2.0e1 -tout -1.0"
#	elif defined(_HEAT_EQ)
        write(arguments, '(A, A)') trim(arguments), "  -dmin 1 -dmax 16 -tsteps -1 -tmax 1.0 -tout -1.0"
#	elif defined(_SWE)
        write(arguments, '(A, A)') trim(arguments), "  -dmin 2 -dmax 14 -tsteps -1 -tmax 2.0e-3 -tout -1.0"
#	elif defined(_NUMA)
        write(arguments, '(A, A)') trim(arguments), "  -dmin 2 -dmax 14 -tsteps -1 -tmax 5 -tout -1.0"
#	elif defined(_TESTS)
        write(arguments, '(A, A)') trim(arguments), "  -dmin 18 -tsteps -20'"
#	else
#       error No scenario selected!
#	endif

	call kracken('samoa', arguments)

	!  get values
	l_help = lget('samoa_-help') .or. lget('samoa_h')
	l_version = lget('samoa_-version') .or. lget('samoa_v')
	i_min_depth = iget('samoa_dmin')
	i_max_depth = iget('samoa_dmax')
	i_max_time_steps = iget('samoa_tsteps')
	r_max_time = rget('samoa_tmax')
	r_output_time_step = rget('samoa_tout')
	l_log = lget('samoa_noprint')
	i_threads = iget('samoa_threads')
	i_sections_per_thread = iget('samoa_sections')
	i_asagi_mode = iget('samoa_asagihints')

	!if the version option was set was called, display program version
	if (l_version) then
		PRINT '(A)', " sam(oa)²: Space filling curves and Adaptive Meshes for Oceanic and Other Applications"
		PRINT '(A, I0, A, I0, A, I0)', " version ", 0, ".", 5, ".", 0
	end if

	!if the help option was set, display the list of arguments
	if (l_help) then
		PRINT '(A)',            " Usage: samoa [--help | -h] | [--version | -v] | [-asagihints <value>] [-dstart <value>] [-dmin <value>] [-dmax <value>] [-tsteps <value>] [-tmax <value>] [-tout <value>] [-threads <value>] [-sections <value>] [-noprint]"
		PRINT '(A)',            ""
		PRINT '(A)',            " Arguments:"
		PRINT '(A, I0, A)',     " 	-asagihints <value>     ASAGI mode (0: default, 1: pass thorugh, 2: nompi, 3: no mpi + small cache, 4: large grid) (value: ", i_asagi_mode, ")"
		PRINT '(A, I0, A)',     " 	-dmin <value>           minimum grid depth (value: ", i_min_depth, ")"
		PRINT '(A, I0, A)',     "	-dmax <value>           maximum grid depth (value: ", i_max_depth, ")"
		PRINT '(A, I0, A)',     "	-tsteps <value>         maximum number of time steps, less than 0: not defined (value: ", i_max_time_steps, ")"
		PRINT '(A, ES8.1, A)',  "	-tmax <value>           maximum simulation time in seconds, less than 0: not defined (value: ", r_max_time, ")"
		PRINT '(A, ES8.1, A)',  "	-tout <value>           output time step in seconds, less than 0: not defined (value: ", r_output_time_step, ")"
		PRINT '(A, I0, A)',     "	-threads <value>        number of OpenMP threads (value: ", i_threads, ")"
		PRINT '(A, I0, A)',     "	-sections <value>       number of grid sections per OpenMP thread (value: ", i_sections_per_thread, ")"
		PRINT '(A)',            "	-noprint                print log to file instead of console"
		PRINT '(A)',            "	--help, -h              display this help and exit"
		PRINT '(A)',            "	--version, -v           output version information and exit"
    end if

    !stop if the version or help command were called
	if (l_help .or. l_version) then
        stop
    end if

#	if defined(_MPI)
#       if defined(_OMP)
            call mpi_init_thread(MPI_THREAD_MULTIPLE, mpi_prov_thread_support, i_error); assert_eq(i_error, 0)
            assert_eq(MPI_THREAD_MULTIPLE, mpi_prov_thread_support)
#       else
            call mpi_init(i_error); assert_eq(i_error, 0)
#       endif

		call mpi_comm_size(MPI_COMM_WORLD, size_MPI, i_error); assert_eq(i_error, 0)
		call mpi_comm_rank(MPI_COMM_WORLD, rank_MPI, i_error); assert_eq(i_error, 0)

		call mpi_comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, mpi_tag_upper_bound, mpi_flag, i_error); assert_eq(i_error, 0)
		assert(mpi_flag)
		assert_ge(mpi_tag_upper_bound, ishft(1, 30) - 1)
#	else
		size_MPI = 1
		rank_MPI = 0
#	endif

	call omp_set_num_threads(i_threads)

! 	 	      H
	!   __________
	!	|\      /|
	!	| \  L / |
	!	|  \  /  |        Unit Square
!  H	|   \/   |H
	!	|   /\   |
	!	|  / L\  |
	!	| /    \ |
	!	|/______\|
!      	     H
	num_initial_triangles = 4

	ALLOCATE(tree(num_initial_triangles), stat=i_error); assert_eq(i_error, 0)

	i_initial_elements = 0
	i_initial_fine_elements = 0

	! The coordinates of the corners in the triangle
	! left node
	r_coords2(1,1,1) = 0.0_GRID_SR
	r_coords2(2,1,1) = 0.0_GRID_SR

	! top node
	r_coords2(1,2,1) = 0.5_GRID_SR
	r_coords2(2,2,1) = 0.5_GRID_SR

	! right node
	r_coords2(1,3,1) = 1.0_GRID_SR
	r_coords2(2,3,1) = 0.0_GRID_SR

	! "norm_vec" for triangle type 1..8
	! leg_size, hypotenuse size, area for triangles up to start_depth 40
	call initialize_triangle_constants(sqrt(dot_product(	&
		r_coords2(1:2, 1, 1) - r_coords2(1:2, 2, 1),	&
		r_coords2(1:2, 1, 1) - r_coords2(1:2, 2, 1)	)))
	p_triangle_tree => tree(1)
	!Omer
	!Add the FEM triangle type to the initial triangles
	p_triangle_tree%i_matrix_type = 5
	l_boundary_data(:,1,1) = internal_edge
	l_boundary_data(:,2,1) = dirichlet_boundary
	if (num_initial_triangles == 3) then		! 3 triangles
		l_boundary_data(:,3,1) = dirichlet_boundary
	else	! 4 triangles
		l_boundary_data(:,3,1) = internal_edge
	endif
	CALL generate_triangle_tree(p_triangle_tree, 1, i_max_depth, r_coords2(:,:,1), l_boundary_data(:,:,1), K, GREEN)
	i_initial_elements = i_initial_elements + p_triangle_tree%i_num_elements
	i_initial_fine_elements = i_initial_fine_elements + p_triangle_tree%i_num_fine_elements

	r_coords2(1,1,2) = 1.0_GRID_SR
	r_coords2(2,1,2) = 0.0_GRID_SR
	r_coords2(1,2,2) = 0.5_GRID_SR
	r_coords2(2,2,2) = 0.5_GRID_SR
	r_coords2(1,3,2) = 1.0_GRID_SR
	r_coords2(2,3,2) = 1.0_GRID_SR
	p_triangle_tree => tree(2)
	!Omer
	!Add the FEM triangle type to the initial triangles
	p_triangle_tree%i_matrix_type = 7
	l_boundary_data(:,1,2) = internal_edge
	l_boundary_data(:,2,2) = dirichlet_boundary
	l_boundary_data(:,3,2) = internal_edge
	CALL generate_triangle_tree(p_triangle_tree, 1, i_max_depth, r_coords2(:,:,2), l_boundary_data(:,:,2), V, GREEN)
	i_initial_elements = i_initial_elements + p_triangle_tree%i_num_elements
	i_initial_fine_elements = i_initial_fine_elements + p_triangle_tree%i_num_fine_elements

	r_coords2(1,1,3) = 1.0_GRID_SR
	r_coords2(2,1,3) = 1.0_GRID_SR
	r_coords2(1,2,3) = 0.5_GRID_SR
	r_coords2(2,2,3) = 0.5_GRID_SR
	r_coords2(1,3,3) = 0.0_GRID_SR
	r_coords2(2,3,3) = 1.0_GRID_SR
	p_triangle_tree => tree(3)
	!Omer
	!Add the FEM triangle type to the initial triangles
	p_triangle_tree%i_matrix_type = 1
	if (num_initial_triangles == 3) then		! 3 triangles
		l_boundary_data(:,1,3) = dirichlet_boundary
	else	! 4 triangles
		l_boundary_data(:,1,3) = internal_edge
	endif
	l_boundary_data(:,2,3) = dirichlet_boundary
	l_boundary_data(:,3,3) = internal_edge
	call generate_triangle_tree(p_triangle_tree, 1, i_max_depth, r_coords2(:,:,3), l_boundary_data(:,:,3), V, GREEN)
	i_initial_elements = i_initial_elements + p_triangle_tree%i_num_elements
	i_initial_fine_elements = i_initial_fine_elements + p_triangle_tree%i_num_fine_elements

	r_coords2(1,1,4) = 0.0_GRID_SR
	r_coords2(2,1,4) = 1.0_GRID_SR
	r_coords2(1,2,4) = 0.5_GRID_SR
	r_coords2(2,2,4) = 0.5_GRID_SR
	r_coords2(1,3,4) = 0.0_GRID_SR
	r_coords2(2,3,4) = 0.0_GRID_SR
	p_triangle_tree => tree(4)
	!Omer
	!Add the FEM triangle type to the initial triangles
	p_triangle_tree%i_matrix_type = 3
	l_boundary_data(:,1,4) = internal_edge
	l_boundary_data(:,2,4) = dirichlet_boundary
	l_boundary_data(:,3,4) = internal_edge
	call generate_triangle_tree(p_triangle_tree, 1, i_max_depth, r_coords2(:,:,4), l_boundary_data(:,:,4), -3_1, GREEN)
	i_initial_elements = i_initial_elements + p_triangle_tree%i_num_elements
	i_initial_fine_elements = i_initial_fine_elements + p_triangle_tree%i_num_fine_elements

    call sfc_generic(tree, num_initial_triangles)

	do i = 1, 4
		call delete_triangle_tree(tree(i))
	end do

	deallocate(tree, stat=i_error); assert_eq(i_error, 0)

#	if defined(_MPI)
		call mpi_barrier(MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
		call mpi_finalize(i_error); assert_eq(i_error, 0)
#	endif

	stop
end PROGRAM gridtest
