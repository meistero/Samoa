! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"


#if defined(_DARCY)

#   define _CG                  darcy_pressure_solver
#   define _CG_mod_types        Samoa_darcy

#   define _GV_NODE_SIZE        _DARCY_P_NODE_SIZE
#   define _GV_EDGE_SIZE        _DARCY_P_EDGE_SIZE
#   define _GV_CELL_SIZE        _DARCY_P_CELL_SIZE

#   define _gm_A                darcy_gm_A
#   define _gv_x                darcy_gv_p

#   define _gv_r                darcy_gv_r
#   define _gv_d                darcy_gv_d
#   define _gv_A_d              darcy_gv_A_d
#   define _gv_trace_A          darcy_gv_mat_diagonal
#   define _gv_dirichlet        darcy_gv_is_dirichlet_boundary

#   include "../Solver/CG.f90"

	MODULE Darcy
		use Darcy_data_types
		use Darcy_initialize_pressure
		use Darcy_initialize_saturation
		use Darcy_edge_dummy
		use Darcy_vtk_output
		use Darcy_xml_output
		use Darcy_grad_p
		use Darcy_laplace_jacobi
		use Darcy_laplace_cg
		use Darcy_transport_eq
		use Darcy_permeability
		use Darcy_adapt

		use linear_solver
        use Darcy_pressure_solver

		use Samoa_darcy

		implicit none

		type t_darcy
            type(t_darcy_init_pressure_traversal)           :: init_pressure
            type(t_darcy_init_saturation_traversal)         :: init_saturation
            type(t_darcy_edge_dummy_traversal)              :: edge_dummy
            type(t_darcy_vtk_output_traversal)              :: vtk_output
            type(t_darcy_xml_output_traversal)              :: xml_output
            type(t_darcy_grad_p_traversal)                  :: grad_p
            type(t_darcy_transport_eq_traversal)            :: transport_eq
            type(t_darcy_permeability_traversal)            :: permeability
            type(t_darcy_adaption_traversal)                :: adaption
            !type(t_darcy_jacobi_solver)                    :: pressure_solver
            !type(t_darcy_cg_solver)                        :: pressure_solver
            class(t_linear_solver), allocatable             :: pressure_solver

            contains

            procedure , pass :: create => darcy_create
            procedure , pass :: run => darcy_run
            procedure , pass :: destroy => darcy_destroy
        end type

		private
		public t_darcy


		contains

		!> Creates all required runtime objects for the scenario
		subroutine darcy_create(darcy, grid, l_log, i_asagi_mode)
            class(t_darcy)                                              :: darcy
 			type(t_grid), intent(inout)									:: grid
			logical (kind = GRID_SL), intent(in)						:: l_log
			integer, intent(in)											:: i_asagi_mode

			!local variables
			type (t_transform_data)										:: dummy_transform_data
			character (len = 64)										:: s_format_string, s_log_name, s_date, s_time
			real (kind = GRID_SR)										:: r_m
			real (kind = GRID_SR), target					            :: default_offset(2)
			integer                                                     :: i_error
            type(t_darcy_pressure_solver)   :: pressure_solver
            type(t_darcy_cg_solver)         :: cg_solver
            type(t_darcy_jacobi_solver)     :: jacobi_solver

            !allocate solver
			grid%r_time = 0.0_GRID_SR
			grid%r_p0 = 1.0e6_GRID_SR          !initial pressure difference in Pa
			grid%r_epsilon = 1.0e-5_GRID_SR
			grid%r_rho = 0.2_GRID_SR
			grid%r_rel_permeability = 1.5_GRID_SR

            pressure_solver = t_darcy_pressure_solver(grid%r_epsilon * grid%r_p0)
            allocate(darcy%pressure_solver, source=cg_solver, stat=i_error); assert_eq(i_error, 0)

			!open log file
			call date_and_time(s_date, s_time)
			write (darcy%vtk_output%s_file_stamp, "(A, A, A8, A, A6)") "output/darcy", "_", s_date, "_", s_time
			write (darcy%xml_output%s_file_stamp, "(A, A, A8, A, A6)") "output/darcy", "_", s_date, "_", s_time
			write (s_log_name, '(A, A)') TRIM(darcy%vtk_output%s_file_stamp), ".log"

			if (l_log) then
				_log_open_file(s_log_name)
			endif

			!create a dunavant quadrature rule for the stiffness matrix (integral over product of 2 derivatives of the basis functions)
			call t_qr_create_dunavant_rule(qr_p, max(1, 2 * _DARCY_P_ORDER - 2))

			!create a custom quadrature rule for the cell boundary
			qr_u_boundary%i_qpts = 3
			allocate(qr_u_boundary%r_qpt_coords(2, 3), stat = i_error); assert_eq(i_error, 0)
			allocate(qr_u_boundary%r_qpt_weights(3), stat = i_error); assert_eq(i_error, 0)
			allocate(qr_u_boundary%r_qpt_normals(2, 3), stat = i_error); assert_eq(i_error, 0)

			qr_u_boundary%r_qpt_coords(:, 1) = [ 0.0_GRID_SR, 0.5_GRID_SR ]
			qr_u_boundary%r_qpt_coords(:, 2) = [ 0.5_GRID_SR, 0.5_GRID_SR ]
			qr_u_boundary%r_qpt_coords(:, 3) = [ 0.5_GRID_SR, 0.0_GRID_SR ]

			qr_u_boundary%r_qpt_normals(:, 1) = [ -1.0_GRID_SR, 0.0_GRID_SR ]
			qr_u_boundary%r_qpt_normals(:, 2) = [ sqrt(0.5_GRID_SR), sqrt(0.5_GRID_SR) ]
			qr_u_boundary%r_qpt_normals(:, 3) = [ 0.0_GRID_SR, -1.0_GRID_SR ]

			qr_u_boundary%r_qpt_weights(1) = 1.0_GRID_SR
			qr_u_boundary%r_qpt_weights(2) = sqrt(2.0_GRID_SR)
			qr_u_boundary%r_qpt_weights(3) = 1.0_GRID_SR

			call t_qr_create(qr_u_boundary)

			write (s_format_string, fmt = '(A, I0, A)') "(", _DARCY_P_SIZE, " (F8.5, 2X))"

			!compute stiffness matrices
			default_offset = [ 0.0_GRID_SR, 0.0_GRID_SR ]
			dummy_transform_data%plotter_data => ref_plotter_data(1)
			dummy_transform_data%custom_data%offset => default_offset
			dummy_transform_data%custom_data%scaling = 1.0_GRID_SR

			call samoa_calc_element_matrix(dummy_transform_data, qr_p%t_qr_base, darcy_stiffness_matrix_op, grid%stiffness_matrix)

			_log_write(2, *) "Stiffness matrix: "
			_log_write(2, s_format_string) grid%stiffness_matrix
			_log_write(2, *) "---"

			call load_permeability_data(grid, i_asagi_mode)
		end subroutine

		subroutine load_permeability_data(grid, i_asagi_mode)
			type(t_grid), target, intent(inout)		:: grid
			character (len = 64)					:: s_file_name
			integer, intent(in)						:: i_asagi_mode

			integer									:: i_error, i, j, k
			integer, pointer						:: afh

#			if defined(_ASAGI)
				afh => grid%afh_permeability
				afh = asagi_create(grid_type = GRID_FLOAT, hint = i_asagi_mode, levels = grid%i_max_depth / 2 + 1)

				do i = 0, grid%i_max_depth / 2
					do j = i, 0, -1
						write (s_file_name, fmt = '(A, I0, A, A)') "data/perm_", 2 ** j, ".nc"

                     	i_error = asagi_open(afh, trim(s_file_name), i)

						if (i_error == GRID_SUCCESS) then
							exit
						endif

						assert_gt(j, 0)
					end do

                    if (rank_MPI == 0) then
                        _log_write(1, '(A, A, A, F0.2, A, F0.2, A, F0.2, A, F0.2, A)') " Darcy: loaded '", trim(s_file_name), "', coordinate system: [", grid_min_x(afh), ", ", grid_min_y(afh), "] x [", grid_max_x(afh), ", ", grid_max_y(afh), "]"
                    end if
				end do
#			endif
		end subroutine

		!> Destroys all required runtime objects for the scenario
		subroutine darcy_destroy(darcy, grid, l_log)
            class(t_darcy)                                               :: darcy
 			type(t_grid), intent(inout)							:: grid
			logical (kind = GRID_SL)		:: l_log
			integer (kind = 1)				:: i

			call t_qr_destroy(qr_p)
			call t_qr_destroy(qr_u_boundary)

			if (l_log) then
				_log_close_file()
			endif

#			if defined(_ASAGI)
				call asagi_close(grid%afh_permeability)
#			endif
		end subroutine

		!> Sets the initial values of the scenario and runs the time steps
		subroutine darcy_run(darcy, grid, i_max_time_steps, r_max_time, r_output_step)
            class(t_darcy)                                              :: darcy
 			type(t_grid), intent(inout)									:: grid
			integer (kind = GRID_SI), intent(inout)						:: i_max_time_steps
			real (kind = GRID_SR), intent(in)							:: r_max_time
			real (kind = GRID_SR), intent(in)							:: r_output_step

			type (t_adaptive_statistics)								:: adaption_stats_initial, adaption_stats_time_steps
			type (t_statistics)									        :: cg_stats_initial, cg_stats_time_steps
			type (t_statistics)									        :: grid_stats_initial, grid_stats_time_steps
            integer (kind = GRID_SI)									:: i_lse_iterations, i_lse_iterations_initial
			double precision										    :: r_t1, r_t2, r_t3, r_t4
			real (kind = GRID_SR)										:: r_time_next_output
			type(t_section_info)           	                            :: grid_info

			!init parameters
			r_time_next_output = 0.0_GRID_SR

            !$omp master
            if (rank_MPI == 0) then
                _log_write(0, *) "Darcy: setting initial values and solving initial system.."
                _log_write(0, *) ""
            end if
            !$omp end master

			r_t1 = omp_get_wtime()

			!set pressure initial condition
			call darcy%init_pressure%traverse(grid)

			do
				!reset saturation to initial condition, compute permeability and set refinement flag
				call darcy%init_saturation%traverse(grid)

				!solve pressure equation
				i_lse_iterations = darcy%pressure_solver%solve(grid)

                grid_info = grid%get_capacity(.false.)

                !$omp master
                if (rank_MPI == 0) then
                    _log_write(1, "(A, I0, A, I0, A, I0, A)") " Darcy: ", darcy%adaption%stats%i_traversals, " adaptions, ", i_lse_iterations, " iterations, ", grid_info%i_cells, " cells"
                end if
                !$omp end master

                grid_info = grid%get_capacity(.true.)
				if (darcy%init_saturation%i_refinements_issued .le. grid_info%i_cells / 100) then
					exit
				endif

				!refine grid
				call darcy%adaption%traverse(grid)
			end do

			r_t2 = omp_get_wtime()

            !$omp master
            if (rank_MPI == 0) then
                _log_write(0, *) "Darcy: done."
                _log_write(0, *) ""

                call grid_info%print()
            end if
			!$omp end master

			!output initial grid
			if (r_output_step >= 0.0_GRID_SR) then
				call darcy%grad_p%traverse(grid)
				call darcy%xml_output%traverse(grid)
				r_time_next_output = r_time_next_output + r_output_step
			end if

            !$omp master
			!copy counters
            adaption_stats_initial = darcy%adaption%stats
            cg_stats_initial = darcy%pressure_solver%stats
			grid_stats_initial = grid%stats
			r_asagi_time_initial = r_asagi_time

            if (rank_MPI == 0) then
                _log_write(0, *) "Darcy: running time steps.."
                _log_write(0, *) ""
            end if
			!$omp end master

			r_t3 = omp_get_wtime()

			do
				if ((r_max_time >= 0.0 .and. grid%r_time >= r_max_time) .or. (i_max_time_steps >= 0 .and. darcy%transport_eq%stats%i_traversals >= i_max_time_steps)) then
					exit
				end if

				!refine grid if necessary
				!if (darcy%permeability%i_refinements_issued > 0) then
					!refine grid
					call darcy%adaption%traverse(grid)

					!recompute permeability + refinement flag
					call darcy%permeability%traverse(grid)
				!end if

				!solve pressure equation
				i_lse_iterations = darcy%pressure_solver%solve(grid)

				!compute velocity field
				call darcy%grad_p%traverse(grid)

				!transport equation time step
				call darcy%transport_eq%traverse(grid)

				!compute permeability field + refinement flag
				call darcy%permeability%traverse(grid)

                grid_info = grid%get_capacity(.false.)

                !$omp master
                if (rank_MPI == 0) then
                    _log_write(1, '(A, I0, A, ES14.7, A, ES14.7, A, I0, A, I0)') " Darcy: time step: ", darcy%transport_eq%stats%i_traversals, ", sim. time:", grid%r_time, " s, dt:", grid%r_dt, " s, cells: ", grid_info%i_cells, ", LSE iterations: ", i_lse_iterations
                end if
                !$omp end master

				!output grid
				if (r_output_step >= 0.0_GRID_SR .and. grid%r_time >= r_time_next_output) then
					call darcy%xml_output%traverse(grid)
					r_time_next_output = r_time_next_output + r_output_step
				end if
			end do

			r_t4 = omp_get_wtime()

            grid_info = grid%get_capacity(.true.)

            !$omp master
            adaption_stats_time_steps = darcy%adaption%stats - adaption_stats_initial
            cg_stats_time_steps = darcy%pressure_solver%stats - cg_stats_initial
			grid_stats_time_steps = grid%stats - grid_stats_initial
			r_asagi_time = r_asagi_time - r_asagi_time_initial

            call darcy%init_saturation%stats%reduce()
            call adaption_stats_initial%reduce()
            call cg_stats_initial%reduce()
            call grid_stats_initial%reduce()

            call darcy%transport_eq%stats%reduce()
            call adaption_stats_time_steps%reduce()
            call cg_stats_time_steps%reduce()
            call grid_stats_time_steps%reduce()

            call reduce(r_asagi_time_initial, MPI_MAX)
            call reduce(r_asagi_time, MPI_MAX)

            if (rank_MPI == 0) then
                _log_write(0, *) "Darcy: done."
                _log_write(0, *) ""

                _log_write(0, *) "Initialization:"
                _log_write(0, '(A, T34, A)') " Init: ", trim(darcy%init_saturation%stats%to_string())
                _log_write(0, '(A, T34, A)') " Adaptions: ", trim(adaption_stats_initial%to_string())
                _log_write(0, '(A, T34, A)') " CG: ", trim(cg_stats_initial%to_string())
                _log_write(0, '(A, T34, A)') " Grid: ", trim(grid_stats_initial%to_string())
                _log_write(0, '(A, T34, F10.4, A)') " Element throughput: ", 1.0e-6 * real(grid_stats_initial%i_traversed_cells, GRID_SR) / (r_t2 - r_t1), " M/s"
                _log_write(0, '(A, T34, F10.4, A)') " Memory throughput: ", real(grid_stats_initial%i_traversed_memory, GRID_SR) / ((1024 * 1024 * 1024) * (r_t2 - r_t1)), " GB/s"
                _log_write(0, '(A, T34, F10.4, A)') " Asagi time:", r_asagi_time_initial, " s"
                _log_write(0, '(A, T34, F10.4, A)') " Initialization time:", r_t2 - r_t1, " s"
                _log_write(0, *) ""
                _log_write(0, *) "Execution:"
                _log_write(0, '(A, T34, A)') " Time steps: ", trim(darcy%transport_eq%stats%to_string())
                _log_write(0, '(A, T34, A)') " Adaptions: ", trim(adaption_stats_time_steps%to_string())
                _log_write(0, '(A, T34, A)') " CG: ", trim(cg_stats_time_steps%to_string())
                _log_write(0, '(A, T34, A)') " Grid: ", trim(grid_stats_time_steps%to_string())
                _log_write(0, '(A, T34, F10.4, A)') " Element throughput: ", 1.0e-6 * real(grid_stats_time_steps%i_traversed_cells, GRID_SR) / (r_t4 - r_t3), " M/s"
                _log_write(0, '(A, T34, F10.4, A)') " Memory throughput: ", real(grid_stats_time_steps%i_traversed_memory, GRID_SR) / ((1024 * 1024 * 1024) * (r_t4 - r_t3)), " GB/s"
                _log_write(0, '(A, T34, F10.4, A)') " Asagi time:", r_asagi_time, " s"
                _log_write(0, '(A, T34, F10.4, A)') " Execution time:", r_t4 - r_t3, " s"
                _log_write(0, *) ""
                _log_write(0, '(A, T34, F10.4, A)') " Total time:", (r_t2 - r_t1) + (r_t4 - r_t3), " s"
                _log_write(0, *) "---"

                call grid_info%print()
            end if
            !$omp end master
		end subroutine

		!*********************************
		! FEM Matrix operands
		!*********************************

		pure function darcy_stiffness_matrix_op(td, qr, i, j, i_qpt) result(op)
			type(t_transform_data), intent(in)				:: td
			type(t_qr_base), intent(in)						:: qr
			integer (kind = GRID_SI), intent(in)			:: i, j
			integer (kind = GRID_SI), intent(in)			:: i_qpt
			real (kind = GRID_SR)							:: op

			!stiffness matrix operator dot(grad psi_i, grad psi_j)
			op = abs(td%plotter_data%det_jacobian) * DOT_PRODUCT( &
				samoa_barycentric_to_world_normal(td, t_qr_base_d_psi_d_lambda(qr, i, i_qpt)), &
				samoa_barycentric_to_world_normal(td, t_qr_base_d_psi_d_lambda(qr, j, i_qpt)) &
			)
		end function
	END MODULE Darcy
#endif

