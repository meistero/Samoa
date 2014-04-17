! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_TESTS)
	MODULE Tests
		use SFC_node_traversal

		use Tests_data_types
		use Tests_basis_functions_mod
		use Tests_initialize
		use Tests_consistency
		use Tests_consistency_ringbuffer
		use Tests_node_dummy
		use Tests_memory

		use Tests_data_types

		PRIVATE
		PUBLIC tests_create, tests_run, tests_destroy

		contains

		!> Creates all required runtime objects for the scenario
		subroutine tests_create(l_log)
			implicit none
			logical 									:: l_log

			!local variables
			CHARACTER(len = 64)											:: s_format_string
			CHARACTER(len = 64)											:: s_log_name
			CHARACTER(len = 16)											:: s_date, s_time

			call date_and_time(s_date, s_time)
			write (tests_global_data%s_file_stamp, "(A, A, A8, A, A6)") "output/tests", "_", s_date, "_", s_time
			write (s_log_name, '(A, A)') TRIM(tests_global_data%s_file_stamp), ".log"

			if (l_log) then
				_log_open_file(s_log_name)
			endif
		end subroutine

		!> Destroys all required runtime objects for the scenario
		subroutine tests_destroy(l_log)
			implicit none
			logical 									:: l_log

			if (l_log) then
				_log_close_file()
			endif
		end subroutine

		!> Sets the initial values of the scenario and runs the time steps
		subroutine tests_run(triangles, num_coarse_triangles, i_triang, l_forward, i_dir, i_time_steps)
			implicit none
			type(triangle_tree), dimension(:), intent(inout), target	:: triangles				!< initial triangle array
			integer (kind = GRID_SI), intent(in)						:: num_coarse_triangles
			integer (kind = GRID_SI), intent(inout)						:: i_triang
			logical , intent(inout)						:: l_forward
			integer (kind = GRID_SI), intent(inout)						:: i_dir
			integer (kind = GRID_SI), intent(inout)						:: i_time_steps

			integer (kind = GRID_SI)									:: i
			real (kind = GRID_SR)										:: r_t1, r_t2

			_log_write(0, '(A, I0, A, I0, A, I0)') " Tests: min depth: ", i_min_depth
			_log_write(0, '(A, I0, A, ES9.2, A, ES9.2)') " Tests: time steps: ", i_time_steps
			_log_write(0, *) ""

			!interpret a negative number of time steps as infinite
			if (i_time_steps < 0) then
				i_time_steps = huge(1_GRID_SI)
			end if

			call tests_init_traversal(triangles, num_coarse_triangles, i_triang, l_forward, i_dir)

			call tests_consistency_traversal(triangles, num_coarse_triangles, i_triang, l_forward, i_dir)
			call tests_consistency_ringbuffer_traversal(triangles, num_coarse_triangles, i_triang, l_forward, i_dir)
			_log_write(0, *) ""



			_log_write(0, *) "Tests: running throughput test for cells:"

			!reset counters
			i_grid_traversals = 0
			i_traversed_elements = 0
			i_traversed_memory = 0

			CALL cpu_time(r_t1)
			do i = 1, i_time_steps
				call tests_memory_traversal_c(triangles, num_coarse_triangles, i_triang, l_forward, i_dir)
			end do
			CALL cpu_time(r_t2)

			_log_write(0, '(A, T34, F10.4, A)') "  Execution time: ", r_t2 - r_t1, " s"
			_log_write(0, '(A, T34, F10.4, A)') "  Element throughput: ", 1.0e-6 * real(i_traversed_elements) / (r_t2 - r_t1), " M/s"
			_log_write(0, '(A, T34, F10.4, A)') "  Memory throughput: ", real(i_traversed_memory) / ((1024 * 1024 * 1024) * (r_t2 - r_t1)), " GB/s"
			_log_write(0, *) ""



			_log_write(0, *) "Tests: running throughput test for cells + edges:"

			!reset counters
			i_grid_traversals = 0
			i_traversed_elements = 0
			i_traversed_memory = 0

			CALL cpu_time(r_t1)
			do i = 1, i_time_steps
				call tests_memory_traversal_ce(triangles, num_coarse_triangles, i_triang, l_forward, i_dir)
			end do
			CALL cpu_time(r_t2)

			_log_write(0, '(A, T34, F10.4, A)') "  Execution time: ", r_t2 - r_t1, " s"
			_log_write(0, '(A, T34, F10.4, A)') "  Element throughput: ", 1.0e-6 * real(i_traversed_elements) / (r_t2 - r_t1), " M/s"
			_log_write(0, '(A, T34, F10.4, A)') "  Memory throughput: ", real(i_traversed_memory) / ((1024 * 1024 * 1024) * (r_t2 - r_t1)), " GB/s"
			_log_write(0, *) ""



			_log_write(0, *) "Tests: running throughput test for cells + nodes:"

			!reset counters
			i_grid_traversals = 0
			i_traversed_elements = 0
			i_traversed_memory = 0

			CALL cpu_time(r_t1)
			do i = 1, i_time_steps
				call tests_memory_traversal_cn(triangles, num_coarse_triangles, i_triang, l_forward, i_dir)
			end do
			CALL cpu_time(r_t2)

			_log_write(0, '(A, T34, F10.4, A)') "  Execution time: ", r_t2 - r_t1, " s"
			_log_write(0, '(A, T34, F10.4, A)') "  Element throughput: ", 1.0e-6 * real(i_traversed_elements) / (r_t2 - r_t1), " M/s"
			_log_write(0, '(A, T34, F10.4, A)') "  Memory throughput: ", real(i_traversed_memory) / ((1024 * 1024 * 1024) * (r_t2 - r_t1)), " GB/s"
			_log_write(0, *) ""


			!if necessary, do a dummy traversal to synchronize edge and node streams
			if (l_forward .xor. all(p_color_edges_in%i_inc > 0)) then
				call tests_node_dummy_traversal(triangles, num_coarse_triangles, i_triang, l_forward, i_dir)
			end if



			_log_write(0, *) "Tests: running throughput test for cells + edges + nodes:"

			!reset counters
			i_grid_traversals = 0
			i_traversed_elements = 0
			i_traversed_memory = 0

			CALL cpu_time(r_t1)
			do i = 1, i_time_steps
				call tests_memory_traversal_cen(triangles, num_coarse_triangles, i_triang, l_forward, i_dir)
			end do
			CALL cpu_time(r_t2)

			_log_write(0, '(A, T34, F10.4, A)') "  Execution time: ", r_t2 - r_t1, " s"
			_log_write(0, '(A, T34, F10.4, A)') "  Element throughput: ", 1.0e-6 * real(i_traversed_elements) / (r_t2 - r_t1), " M/s"
			_log_write(0, '(A, T34, F10.4, A)') "  Memory throughput: ", real(i_traversed_memory) / ((1024 * 1024 * 1024) * (r_t2 - r_t1)), " GB/s"
			_log_write(0, *) ""



			call tests_basis_functions()

			_log_write(0, *) "Tests: done."
			_log_write(0, *) ""

			call print_memory_stats()
		end subroutine
	END MODULE Tests
#endif

