! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_HEAT_EQ)
	MODULE Heat_Eq
		use SFC_node_traversal
		use Heat_Eq_data_types

		use Heat_Eq_adapt
		use Heat_Eq_initialize
		use Heat_Eq_output
		use Heat_Eq_xml_output
		use Heat_Eq_euler_Timestep
		use Heat_Eq_midpoint_Timestep
		use Heat_Eq_heun_Timestep

		use Samoa_heat_eq

		implicit none

		PRIVATE
		PUBLIC heat_eq_create, heat_eq_run, heat_eq_destroy

		type(samoa_qr_T)												:: qr_T

		contains

		!> Creates all required runtime objects for the scenario
		subroutine heat_eq_create(grid, l_log, i_asagi_mode)
			type(t_grid), intent(inout)									:: grid
			logical 									:: l_log
			integer														:: i_asagi_mode

			!local variables
			real (kind = GRID_SR), dimension(2), target					:: default_offset
			real(kind = GRID_SR)										:: r_m
			type(t_transform_data)										:: dummy_transform_data
			CHARACTER(64)												:: s_format_string, s_log_name, s_date, s_time
            integer(kind = BYTE)                                           :: i, j

			!open log file
			call date_and_time(s_date, s_time)
			write (grid%s_file_stamp, "(A, A, A8, A, A6)") "output/heat_eq", "_", s_date, "_", s_time
			write (s_log_name, '(A, A)') TRIM(grid%s_file_stamp), ".log"

			if (l_log) then
				_log_open_file(s_log_name)
			endif

			!create a quadrature rule for the mass & stiffness matrix (integral over product of 2 derivatives of the basis functions)
			call t_qr_create_dunavant_rule(qr_T, 2 * _HEAT_EQ_ORDER)

			write (s_format_string, fmt = '(A, I0, A)') "(", _HEAT_EQ_SIZE, " (F8.5, 2X))"

			default_offset = [ 0.0_GRID_SR, 0.0_GRID_SR ]
			dummy_transform_data%plotter_data => ref_plotter_data(1)
			dummy_transform_data%custom_data%offset => default_offset
			dummy_transform_data%custom_data%scaling = 1.0_GRID_SR

			call samoa_calc_element_matrix(dummy_transform_data, qr_T%t_qr_base, sc_mass_matrix_op, gm_mass_matrix)

			_log_write(0, *) "Mass matrix: "
			_log_write(0, s_format_string) gm_mass_matrix
			_log_write(0, *) "---"

			!lump the mass matrix
#			if (_HEAT_EQ_ORDER > 0)
				do i = 1, 3
					r_m = 0.0_GRID_SR

					do j = 1, 3
						r_m = r_m + gm_mass_matrix(i, j)
					end do

					gm_mass_matrix(i, i) = r_m
				end do
#			endif

#			if (_HEAT_EQ_ORDER > 1)
				do i = 4, 6
					r_m = 0.0_GRID_SR

					do j = 4, 6
						r_m = r_m + gm_mass_matrix(i, j)
					end do

					gm_mass_matrix(i, i) = r_m
				end do
#			endif

			forall (i = 1 :_HEAT_EQ_SIZE, j = 1 : _HEAT_EQ_SIZE, i /= j)
                gm_mass_matrix(i, j) = 0.0_GRID_SR
			end forall

			_log_write(0, *) "Lumped mass matrix: "
			_log_write(0, s_format_string) gm_mass_matrix
			_log_write(0, *) "---"

			call samoa_calc_element_matrix(dummy_transform_data, qr_T%t_qr_base, sc_stiffness_matrix_op, gm_stiffness_matrix)

			_log_write(0, *) "Stiffness matrix: "
			_log_write(0, s_format_string) gm_stiffness_matrix
			_log_write(0, *) "---"

			grid%r_laser_rps = 1.0_GRID_SR
		end subroutine

		!> Destroys all required runtime objects for the scenario
		subroutine heat_eq_destroy(grid, l_log)
			type(t_grid), intent(inout)									:: grid
			logical 		:: l_log

			if (l_log) then
				_log_close_file()
			endif

			call t_qr_destroy(qr_T)
		end subroutine

		!*********************************
		! FEM Matrix operands
		!*********************************

		pure function sc_mass_matrix_op(td, qr, i, j, i_qpt) result(op)
			type(t_transform_data), intent(in)				:: td
			type(t_qr_base), intent(in)						:: qr
			integer (kind = GRID_SI), intent(in)			:: i, j
			integer (kind = GRID_SI), intent(in)			:: i_qpt
			real (kind = GRID_SR)							:: op

			!mass matrix operator phi_i * psi_j
			op = qr%r_dof_qpt_psi(i, i_qpt) * qr%r_dof_qpt_psi(j, i_qpt)
		end function

		pure function sc_stiffness_matrix_op(td, qr, i, j, i_qpt) result(op)
			type(t_transform_data), intent(in)				:: td
			type(t_qr_base), intent(in)						:: qr
			integer (kind = GRID_SI), intent(in)			:: i, j
			integer (kind = GRID_SI), intent(in)			:: i_qpt
			real (kind = GRID_SR)							:: op

			!stiffness matrix operator <J^(-T) grad (phi_i), J^(-T) grad psi_j>
			op = abs(td%plotter_data%det_jacobian) * DOT_PRODUCT( &
				samoa_barycentric_to_world_normal(td, t_qr_base_d_psi_d_lambda(qr, i, i_qpt)), &
				samoa_barycentric_to_world_normal(td, t_qr_base_d_psi_d_lambda(qr, j, i_qpt)) &
			)
		end function

		!*********************************
		! run()-method
		!*********************************

		!> Sets the initial values of the Heat_Eq and runs the time steps
		!> @param i_output_step_interval		if < 0 no output is produced, if 0 only the initial data is written out, otherwise each n-th time step is written out.
		subroutine heat_eq_run(grid)
			type(t_grid), intent(inout)									:: grid

			type(t_grid)												:: grid_temp
			integer (kind = GRID_SI)									:: i_time_step
			integer (kind = GRID_SI)									:: i_adaptions, i_adaptions_initial, i_adaptions_time_steps
			integer (kind = GRID_SI)									:: i_grid_traversals_initial
			integer (kind = GRID_DI)									:: i_traversed_elements_initial, i_traversed_memory_initial
			real (kind = GRID_SR)										:: r_t1, r_t2, r_t3, r_t4
			real (kind = GRID_SR)										:: r_time_next_output
			type(t_section_info)           	                            :: grid_info

			_log_write(0, '(A, I0, A, I0)') " Heat_Eq: min depth: ", cfg%i_min_depth, ", max depth: ", cfg%i_max_depth
			_log_write(0, '(A, I0, A, ES9.2, A, ES9.2)') " Heat_Eq: max time steps: ", i_max_time_steps, ", max sim. time: ", r_max_time, ", output step: ", cfg%r_output_time_step
			_log_write(0, *) ""

			_log_write(0, *) "Heat_Eq: setting initial values and a priori refinement.."
			_log_write(0, *) ""

			r_t1 = get_wtime()

			grid%r_time = 0.0_GRID_SR
			r_time_next_output = 0.0_GRID_SR
			grid%r_dt = 3.0 * get_cell_volume(cfg%i_max_depth)

			!init counters
			i_grid_traversals = 0
			i_traversed_elements = 0
			i_traversed_memory = 0

			i_adaptions_initial = 0
			do
				!set numerics and check for refinement
				call heat_eq_init_traversal(grid)

                grid_info = grid%get_info()
                !$omp master
				_log_write(1, "(A, I0, A, I0, A)") " Heat_Eq: ", i_adaptions_initial, " adaptions, ", grid_info%i_cells, " cells"
                !$omp end master

				if (heat_eq_init_count_issued_refinements()  .le. grid_info%i_cells / 100) then
					exit
				endif

				call heat_eq_adaption_traversal(grid, grid_temp)
				call grid_temp%move(grid)

				i_adaptions_initial = i_adaptions_initial + 1
			end do

			r_t2 = get_wtime()

			_log_write(0, *) "Heat_Eq: done."
			_log_write(0, *) ""

            grid_info = grid%get_info()
			call grid_info%print()

			!output initial grid
			if (cfg%r_output_time_step >= 0.0_GRID_SR) then
				call heat_eq_xml_output_traversal(grid)
				r_time_next_output = r_time_next_output + cfg%r_output_time_step
			end if

			!reset counters
			i_grid_traversals_initial = i_grid_traversals
			i_traversed_elements_initial = i_traversed_elements
			i_traversed_memory_initial = i_traversed_memory

			i_adaptions_time_steps = 0

 			i_time_step = 1
			i_grid_traversals = 0
			i_traversed_elements = 0
			i_traversed_memory = 0

			_log_write(0, *) "Heat_Eq: running time steps.."
			_log_write(0, *) ""

			r_t3 = get_wtime()

			do
				if ((r_max_time >= 0.0 .and. grid%r_time >= r_max_time) .or. (i_max_time_steps >= 0 .and. i_time_step >= i_max_time_steps)) then
					exit
				end if

				if (timestep_count_issued_refinements() > 0) then
					call heat_eq_adaption_traversal(grid, grid_temp)
					call grid_temp%move(grid)

					i_adaptions_time_steps = i_adaptions_time_steps + 1
				end if

				!do a timestep
				call heun_timestep_traversal(grid)
				grid%r_time = grid%r_time + grid%r_dt

                grid_info = grid%get_info()
				_log_write(1, '(A, I0, A, ES14.7, A, ES14.7, A, I0)') " Heat_Eq: time step: ", i_time_step, ", sim. time:", grid%r_time, " s, dt:", grid%r_dt, " s, cells: ", grid_info%i_cells

				!output grid
				if (cfg%r_output_time_step >= 0.0_GRID_SR .and. grid%r_time >= r_time_next_output) then
					call heat_eq_xml_output_traversal(grid)
					r_time_next_output = r_time_next_output + cfg%r_output_time_step
				end if

				i_time_step = i_time_step + 1
			end do

			r_t4 = get_wtime()

			_log_write(0, *) "Heat_Eq: done."
			_log_write(0, *) ""

			_log_write(0, *) "Initialization:"
			_log_write(0, '(A, T34, I10)') " #Traversals: ", i_grid_traversals_initial
			_log_write(0, '(A, T34, I10)') " #Adaptions: ", i_adaptions_initial
			_log_write(0, '(A, T34, F10.4, A)') " Element throughput: ", 1.0e-6 * real(i_traversed_elements_initial) / (r_t2 - r_t1), " M/s"
			_log_write(0, '(A, T34, F10.4, A)') " Memory throughput: ", real(i_traversed_memory_initial) / ((1024 * 1024 * 1024) * (r_t2 - r_t1)), " GB/s"
			_log_write(0, '(A, T34, F10.4, A)') " Initialization time:", r_t2 - r_t1, " s"
			_log_write(0, *) ""
			_log_write(0, *) "Execution:"
			_log_write(0, '(A, T34, I10)') " #Time steps: ", i_time_step - 1
			_log_write(0, '(A, T34, I10)') " #Traversals: ", i_grid_traversals
			_log_write(0, '(A, T34, I10)') " #Adaptions: ", i_adaptions_time_steps
			_log_write(0, '(A, T34, F10.4, A)') " Element throughput: ", 1.0e-6 * real(i_traversed_elements) / (r_t4 - r_t3), " M/s"
			_log_write(0, '(A, T34, F10.4, A)') " Memory throughput: ", real(i_traversed_memory) / ((1024 * 1024 * 1024) * (r_t4 - r_t3)), " GB/s"
			_log_write(0, '(A, T34, F10.4, A, F10.4, A)') " Execution time (per step):", r_t4 - r_t3, " s      (", (r_t4 - r_t3) / real(i_time_step - 1), " s)"
			_log_write(0, *) ""
			_log_write(0, '(A, T34, F10.4, A)') " Total time:", (r_t2 - r_t1) + (r_t4 - r_t3), " s"
			_log_write(0, *) "---"

            grid_info = grid%get_info()
			call grid_info%print()
		end subroutine
	END MODULE Heat_Eq
#endif

