! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_DARCY)
	MODULE Darcy_xml_output
		use LIB_VTK_IO

		use SFC_edge_traversal
		use Darcy_grad_p

		use Samoa_darcy

		!> Output point data
		type t_output_point_data
			real (kind = GRID_SR), allocatable				    :: coords(:, :)
			real (kind = GRID_SR), allocatable			        :: p(:)
			real (kind = GRID_SR), allocatable			        :: S(:)
			real (kind = GRID_SR), allocatable			        :: rhs(:)

			contains

			procedure, pass :: create => point_data_create
			procedure, pass :: destroy => point_data_destroy
		end type

		!> Output cell data
		type t_output_cell_data
			real (kind = GRID_SR), allocatable			        :: permeability(:, :)
			real (kind = GRID_SR), allocatable			        :: porosity(:)
			real (kind = GRID_SR), allocatable					:: u(:, :)
            integer (kind = GRID_SI), allocatable			    :: rank(:)
            integer (kind = GRID_SI), allocatable			    :: section_index(:)
			integer (BYTE), allocatable					        :: depth(:)
			integer (BYTE), allocatable					        :: refinement(:)

			contains

			procedure, pass :: create => cell_data_create
			procedure, pass :: destroy => cell_data_destroy
		end type

        type num_traversal_data
            type(t_output_point_data)                               :: point_data
            type(t_output_cell_data)                                :: cell_data
            integer (kind = GRID_SI), allocatable			        :: i_connectivity(:)
            character(len=64)							            :: s_file_stamp

            integer (kind = GRID_SI)								:: i_output_iteration = 0
            integer (kind = GRID_SI)								:: i_point_data_index
            integer (kind = GRID_SI)								:: i_cell_data_index
        end type

        interface node_first_touch_op
            module procedure node_first_touch_op_scalar
            module procedure node_first_touch_op_array
        end interface

        type(darcy_gv_p)										    :: gv_p
        type(darcy_gv_rhs)										    :: gv_rhs
        type(darcy_gv_saturation)								    :: gv_saturation
        type(darcy_gv_r)										    :: gv_r

#		define	_GT_NAME							t_darcy_xml_output_traversal

#		if (_DARCY_P_EDGE_SIZE > 0)
#			define _GT_EDGES
#		endif

#		define	_GT_NODES

#		define	_GT_PRE_TRAVERSAL_GRID_OP			pre_traversal_grid_op
#		define	_GT_POST_TRAVERSAL_GRID_OP			post_traversal_grid_op
#		define	_GT_PRE_TRAVERSAL_OP				pre_traversal_op
#		define	_GT_POST_TRAVERSAL_OP				post_traversal_op

#		define	_GT_ELEMENT_OP						element_op

#		define  _GT_NODE_FIRST_TOUCH_OP			    node_first_touch_op

#		include "SFC_generic_traversal_ringbuffer.f90"

		subroutine pre_traversal_grid_op(traversal, grid)
			type(t_darcy_xml_output_traversal), intent(inout)		:: traversal
			type(t_grid), intent(inout)							    :: grid

            if (rank_MPI == 0) then
                _log_write(1, '(A, I0, A, I0)') " Darcy: output step: ", traversal%i_output_iteration
            end if

            call scatter(traversal%s_file_stamp, traversal%children%s_file_stamp)
            call scatter(traversal%i_output_iteration, traversal%children%i_output_iteration)
		end subroutine

		subroutine post_traversal_grid_op(traversal, grid)
			type(t_darcy_xml_output_traversal), intent(inout)		:: traversal
			type(t_grid), intent(inout)							    :: grid

			character (len = 64)							:: s_file_name
            integer                                         :: i_error
			integer(4)										:: i_rank, i_section, e_io
			logical                                         :: l_exists
            type(t_vtk_writer)                              :: vtk

#           if defined(_MPI)
                call mpi_barrier(MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
#           endif

#           if defined(_QUAD_PRECISION)
#               warning VTK output does not work for quad precision
#           else
                if (rank_MPI == 0) then
                    write (s_file_name, "(A, A, I0, A, I0, A)") TRIM(traversal%s_file_stamp), "_", traversal%i_output_iteration, ".pvtu"

                    e_io = vtk%VTK_INI_XML('ascii', s_file_name, 'PUnstructuredGrid')
                        e_io = vtk%VTK_DAT_XML('pfield', 'OPEN')
                            e_io = vtk%VTK_VAR_XML('time', 1.0_GRID_SR, 1)
                        e_io = vtk%VTK_DAT_XML('pfield', 'CLOSE')

                        e_io = vtk%VTK_DAT_XML('pnode', 'OPEN')
                            e_io = vtk%VTK_VAR_XML('pressure', 1.0_GRID_SR, 1)
                            e_io = vtk%VTK_VAR_XML('saturation', 1.0_GRID_SR, 1)
                            e_io = vtk%VTK_VAR_XML('rhs', 1.0_GRID_SR, 1)
                        e_io = vtk%VTK_DAT_XML('pnode', 'CLOSE')

                        e_io = vtk%VTK_DAT_XML('pcell', 'OPEN')
                            e_io = vtk%VTK_VAR_XML('permeability', 1.0_GRID_SR, 3)
                            e_io = vtk%VTK_VAR_XML('porosity', 1.0_GRID_SR, 1)
                            e_io = vtk%VTK_VAR_XML('velocity', 1.0_GRID_SR, 3)
                            e_io = vtk%VTK_VAR_XML('rank', 1_GRID_SI, 1)
                            e_io = vtk%VTK_VAR_XML('section index', 1_GRID_SI, 1)
                            e_io = vtk%VTK_VAR_XML('depth', 1_1, 1)
                            e_io = vtk%VTK_VAR_XML('refinement flag', 1_1, 1)
                        e_io = vtk%VTK_DAT_XML('pcell', 'CLOSE')

                        e_io = vtk%VTK_GEO_XML(1.0_GRID_SR)

                        do i_rank = 0, size_MPI
                            do i_section = 1, huge(1)
                                write (s_file_name, "(A, A, I0, A, I0, A, I0, A)") trim(traversal%s_file_stamp), "_", traversal%i_output_iteration, "_r", i_rank, "_s", i_section, ".vtu"
                                inquire(file = s_file_name, exist = l_exists)

                                if (l_exists) then
                                    write(s_file_name, "(A)") trim(s_file_name(scan(s_file_name, "/\", .true.) + 1 : len(s_file_name)))
                                    e_io = vtk%VTK_GEO_XML(s_file_name)
                                else
                                    exit
                                end if
                            end do
                        end do

                    e_io = vtk%VTK_END_XML()
                end if
#           endif

            traversal%i_output_iteration = traversal%i_output_iteration + 1
		end subroutine

		subroutine pre_traversal_op(traversal, section)
 			type(t_darcy_xml_output_traversal), intent(inout)			:: traversal
 			type(t_grid_section), intent(inout)							:: section

			type(t_section_info)                                        :: grid_info
			integer (kind = GRID_SI)									:: i_error, i_cells, i_points

            grid_info = section%get_info()
			i_cells = max(1, _DARCY_LAYERS) * grid_info%i_cells

			i_points = (_DARCY_LAYERS + 1) * (grid_info%i_nodes + sum(grid_info%i_boundary_nodes))

#           if (_DARCY_LAYERS > 0)
                allocate(traversal%i_connectivity(6 * i_cells), stat = i_error); assert_eq(i_error, 0)
#           else
                allocate(traversal%i_connectivity(3 * i_cells), stat = i_error); assert_eq(i_error, 0)
#           endif
			call traversal%point_data%create(i_points)
			call traversal%cell_data%create(i_cells)

			traversal%i_cell_data_index = 1
			traversal%i_point_data_index = 1
		end subroutine

		subroutine post_traversal_op(traversal, section)
 			type(t_darcy_xml_output_traversal), intent(inout)			:: traversal
			type(t_grid_section), intent(inout)							:: section

			integer (kind = GRID_SI), dimension(:), allocatable			:: i_offsets
			integer (1), dimension(:), allocatable						:: i_types
            type(t_vtk_writer)                                          :: vtk

			type(t_section_info)                                        :: grid_info
			integer (kind = GRID_SI)									:: i_error, i_cells, i_points
			integer(4)													:: e_io, i

            grid_info = section%get_info()

			i_cells = max(1, _DARCY_LAYERS) * grid_info%i_cells
			i_points = (_DARCY_LAYERS + 1) * (grid_info%i_nodes + sum(grid_info%i_boundary_nodes))

			allocate(i_offsets(i_cells), stat = i_error); assert_eq(i_error, 0)
			allocate(i_types(i_cells), stat = i_error); assert_eq(i_error, 0)

#           if (_DARCY_LAYERS > 0)
                i_types = 13_1

                forall (i = 1 : i_cells)
                    i_offsets(i) = 6_GRID_SI * i
                end forall
#           else
                i_types = 5_1

                forall (i = 1 : i_cells)
                    i_offsets(i) = 3_GRID_SI * i
                end forall
#           endif

			write (traversal%s_file_stamp, "(A, A, I0, A, I0, A, I0, A)") TRIM(traversal%s_file_stamp), "_", traversal%i_output_iteration, "_r", rank_MPI, "_s", section%index, ".vtu"

#           if defined(_QUAD_PRECISION)
#               warning VTK output does not work for quad precision
#           else
                e_io = vtk%VTK_INI_XML('binary', traversal%s_file_stamp, 'UnstructuredGrid')
                    e_io = vtk%VTK_DAT_XML('field', 'OPEN')
                        e_io = vtk%VTK_VAR_XML(1, 'time', [section%r_time])
                    e_io = vtk%VTK_DAT_XML('field', 'CLOSE')

                    e_io = vtk%VTK_GEO_XML(i_points, i_cells, traversal%point_data%coords(:, 1), traversal%point_data%coords(:, 2), traversal%point_data%coords(:, 3))

                    e_io = vtk%VTK_CON_XML(i_cells, traversal%i_connectivity, i_offsets, i_types)

                    e_io = vtk%VTK_DAT_XML('node', 'OPEN')
                        e_io = vtk%VTK_VAR_XML(i_points, 'pressure', traversal%point_data%p)
                        e_io = vtk%VTK_VAR_XML(i_points, 'saturation', traversal%point_data%S)
                        e_io = vtk%VTK_VAR_XML(i_points, 'rhs',  traversal%point_data%rhs)
                    e_io = vtk%VTK_DAT_XML('node', 'CLOSE')

                    e_io = vtk%VTK_DAT_XML('cell', 'OPEN')
                        e_io = vtk%VTK_VAR_XML(i_cells, 'permeability', traversal%cell_data%permeability(:, 1), traversal%cell_data%permeability(:, 2), traversal%cell_data%permeability(:, 3))
                        e_io = vtk%VTK_VAR_XML(i_cells, 'porosity', traversal%cell_data%porosity)
                        e_io = vtk%VTK_VAR_XML(i_cells, 'velocity', traversal%cell_data%u(:, 1), traversal%cell_data%u(:, 2), traversal%cell_data%u(:, 3))
                        e_io = vtk%VTK_VAR_XML(i_cells, 'rank', traversal%cell_data%rank)
                        e_io = vtk%VTK_VAR_XML(i_cells, 'section index', traversal%cell_data%section_index)
                        e_io = vtk%VTK_VAR_XML(i_cells, 'depth', traversal%cell_data%depth)
                        e_io = vtk%VTK_VAR_XML(i_cells, 'refinement flag', traversal%cell_data%refinement)
                    e_io = vtk%VTK_DAT_XML('cell', 'CLOSE')

                    e_io = vtk%VTK_GEO_XML()
                e_io = vtk%VTK_END_XML()
#           endif

			deallocate(i_offsets, stat = i_error); assert_eq(i_error, 0)
			deallocate(i_types, stat = i_error); assert_eq(i_error, 0)

			deallocate(traversal%i_connectivity, stat = i_error); assert_eq(i_error, 0)

			call traversal%point_data%destroy()
			call traversal%cell_data%destroy()

			traversal%i_output_iteration = traversal%i_output_iteration + 1
		end subroutine

		!******************
		!Geometry operators
		!******************

		subroutine element_op(traversal, section, element)
 			type(t_darcy_xml_output_traversal), intent(inout)	:: traversal
			type(t_grid_section), intent(inout)					:: section
			type(t_element_base), intent(inout)					:: element

			!local variables

			integer (kind = GRID_SI)							:: i, layer

			!local variables

			real (kind = GRID_SR)               :: p(_DARCY_LAYERS + 1, 3)
			real (kind = GRID_SR)               :: rhs(_DARCY_LAYERS + 1, 3)
			real (kind = GRID_SR)				:: saturation(_DARCY_LAYERS + 1, 3)

			real (kind = GRID_SR)               :: r_point_data_indices(_DARCY_LAYERS + 1, 3)		!< point data indices
			integer (kind = GRID_SI)            :: point_data_indices(3)		!< point data indices

            real (kind = GRID_SR)				:: edge_length, dz
            real (kind = SR)                    :: u_w(3), u_n(3)

#           if (_DARCY_LAYERS > 0)
                real (kind = GRID_SR)			:: g_local(3)
#           else
                real (kind = GRID_SR)			:: g_local(2)
#           endif

			call gv_p%read_from_element(element, p)

			call gv_rhs%read_from_element(element, rhs)
			call gv_r%read_from_element(element, r_point_data_indices)
			call gv_saturation%read_from_element(element, saturation)

            !rotate g so it points in the right direction (no scaling!)
            g_local = g
            g_local(1:2) = samoa_world_to_barycentric_normal(element%transform_data, g_local(1:2))
            g_local(1:2) = g_local(1:2) / (element%transform_data%custom_data%scaling * sqrt(abs(element%transform_data%plotter_data%det_jacobian)))

#           if (_DARCY_LAYERS > 0)
                edge_length = element%cell%geometry%get_leg_size()
                dz = cfg%dz

                do layer = 1, _DARCY_LAYERS
                    traversal%cell_data%rank(traversal%i_cell_data_index) = rank_MPI
                    traversal%cell_data%section_index(traversal%i_cell_data_index) = section%index
                    traversal%cell_data%permeability(traversal%i_cell_data_index, 1:2) = element%cell%data_pers%base_permeability(layer, 1) / 9.869233e-16_SR * (cfg%scaling ** 2)
                    traversal%cell_data%permeability(traversal%i_cell_data_index, 3) = element%cell%data_pers%base_permeability(layer, 2) / 9.869233e-16_SR * (cfg%scaling ** 2)
                    traversal%cell_data%porosity(traversal%i_cell_data_index) = element%cell%data_pers%porosity(layer)

                    call compute_velocity_1D(edge_length, 1.0_SR, element%cell%data_pers%base_permeability(layer, 1), 0.5_SR * (p(layer, 2) + p(layer + 1, 2)), 0.5_SR * (p(layer, 1) + p(layer + 1, 1)), u_w(1), u_n(1), g_local(1))
                    call compute_velocity_1D(edge_length, 1.0_SR, element%cell%data_pers%base_permeability(layer, 1), 0.5_SR * (p(layer, 2) + p(layer + 1, 2)), 0.5_SR * (p(layer, 3) + p(layer + 1, 3)), u_w(2), u_n(2), g_local(2))
                    call compute_velocity_1D(dz, 1.0_SR, element%cell%data_pers%base_permeability(layer, 2), &
                    (0.25_SR * p(layer, 1) + 0.5_SR * p(layer, 2) + 0.25_SR * p(layer, 3)), &
                    (0.25_SR * p(layer + 1, 1) + 0.5_SR * p(layer + 1, 2) + 0.25_SR * p(layer + 1, 3)), u_w(3), u_n(3), g_local(3))

                    u_w(1:2) = samoa_barycentric_to_world_normal(element%transform_data, u_w(1:2))
                    u_w(1:2) = u_w(1:2) * (element%transform_data%custom_data%scaling * sqrt(abs(element%transform_data%plotter_data%det_jacobian)))

                    traversal%cell_data%u(traversal%i_cell_data_index, :) = cfg%scaling * u_w

                    traversal%cell_data%depth(traversal%i_cell_data_index) = element%cell%geometry%i_depth
                    traversal%cell_data%refinement(traversal%i_cell_data_index) = element%cell%geometry%refinement

                    point_data_indices = int(r_point_data_indices(layer, :), kind=GRID_SI)
                    traversal%i_connectivity(6 * traversal%i_cell_data_index - 5 : 6 * traversal%i_cell_data_index - 3) = point_data_indices - 1

                    point_data_indices = int(r_point_data_indices(layer + 1, :), kind=GRID_SI)
                    traversal%i_connectivity(6 * traversal%i_cell_data_index - 2 : 6 * traversal%i_cell_data_index) = point_data_indices - 1

                    traversal%i_cell_data_index = traversal%i_cell_data_index + 1
                end do

                do layer = 1, _DARCY_LAYERS + 1
                    point_data_indices(:) = int(r_point_data_indices(layer, :), kind=GRID_SI)

                    forall (i = 1 : 3)
                        traversal%point_data%coords(point_data_indices(i), 1:2) = cfg%scaling * samoa_barycentric_to_world_point(element%transform_data, samoa_basis_p_get_dof_coords(i)) + cfg%offset
                        traversal%point_data%coords(point_data_indices(i), 3) = cfg%scaling * real(layer - 1, SR) * cfg%dz
                        traversal%point_data%p(point_data_indices(i)) = p(layer, i) / (cfg%scaling * 6894.75729_SR)
                        traversal%point_data%rhs(point_data_indices(i)) = rhs(layer, i) * (cfg%scaling ** 3)
                        traversal%point_data%S(point_data_indices(i)) = saturation(layer, i)
                    end forall
                end do
#           else
                edge_length = element%cell%geometry%get_leg_size()

                point_data_indices(:) = int(r_point_data_indices(1, :), kind=GRID_SI)

                traversal%cell_data%rank(traversal%i_cell_data_index) = rank_MPI
                traversal%cell_data%section_index(traversal%i_cell_data_index) = section%index
                traversal%cell_data%permeability(traversal%i_cell_data_index, :) = element%cell%data_pers%base_permeability / 9.869233e-16_SR * (cfg%scaling ** 2)
                traversal%cell_data%porosity(traversal%i_cell_data_index) = element%cell%data_pers%porosity

                call compute_velocity_1D(edge_length, 1.0_SR, element%cell%data_pers%base_permeability, p(1, 2), p(1, 1), u_w(1), u_n(1), g_local(1))
                call compute_velocity_1D(edge_length, 1.0_SR, element%cell%data_pers%base_permeability, p(1, 2), p(1, 3), u_w(2), u_n(2), g_local(2))

                u_w(1:2) = samoa_barycentric_to_world_normal(element%transform_data, u_w(1:2))
                u_w(1:2) = u_w(1:2) * (element%transform_data%custom_data%scaling * sqrt(abs(element%transform_data%plotter_data%det_jacobian)))
                u_w(3) = 0.0_SR

                traversal%cell_data%u(traversal%i_cell_data_index, :) = cfg%scaling * u_w

                traversal%cell_data%depth(traversal%i_cell_data_index) = element%cell%geometry%i_depth
                traversal%cell_data%refinement(traversal%i_cell_data_index) = element%cell%geometry%refinement

                forall (i = 1 : 3)
                    traversal%point_data%coords(point_data_indices(i), 1:2) = cfg%scaling * samoa_barycentric_to_world_point(element%transform_data, samoa_basis_p_get_dof_coords(i)) + cfg%offset
                    traversal%point_data%coords(point_data_indices(i), 3) = 0.0_SR
                    traversal%point_data%p(point_data_indices(i)) = p(1, i) / (cfg%scaling * 6894.75729_SR)
                    traversal%point_data%rhs(point_data_indices(i)) = rhs(1, i) * (cfg%scaling ** 2)
                    traversal%point_data%S(point_data_indices(i)) = saturation(1, i)
                end forall

                traversal%i_connectivity(3 * traversal%i_cell_data_index - 2 : 3 * traversal%i_cell_data_index) = point_data_indices(1:3) - 1

                traversal%i_cell_data_index = traversal%i_cell_data_index + 1
#           endif
		end subroutine

		subroutine node_first_touch_op_array(traversal, section, nodes)
 			type(t_darcy_xml_output_traversal), intent(inout)	:: traversal
 			type(t_grid_section), intent(in)				    :: section
			type(t_node_data), intent(inout)				    :: nodes(:)

			integer (kind = GRID_SI)						    :: j

            do j = 1, size(nodes)
                call node_first_touch_op_scalar(traversal, section, nodes(j))
            end do
		end subroutine

		subroutine node_first_touch_op_scalar(traversal, section, node)
 			type(t_darcy_xml_output_traversal), intent(inout)	:: traversal
 			type(t_grid_section), intent(in)				    :: section
			type(t_node_data), intent(inout)				    :: node

			integer (kind = GRID_SI) :: i

            do i = 1, _DARCY_LAYERS + 1
                call pre_dof_op(traversal%i_point_data_index, node%data_pers%r(i))
            end do
		end subroutine

		elemental subroutine pre_dof_op(i_point_data_index, r)
 			integer(kind = GRID_SI), intent(inout)			:: i_point_data_index
			real(kind = GRID_SR), intent(out)				:: r

			r = real(i_point_data_index, GRID_SR)
			i_point_data_index = i_point_data_index + 1
		end subroutine


		subroutine point_data_create(point_data, i_points)
			class(t_output_point_data), intent(inout)		:: point_data
			integer(kind = GRID_SI), intent(in)			    :: i_points

			integer (kind = GRID_SI)						:: i_error

			allocate(point_data%coords(i_points, 3), stat = i_error); assert_eq(i_error, 0)
			allocate(point_data%p(i_points), stat = i_error); assert_eq(i_error, 0)
			allocate(point_data%S(i_points), stat = i_error); assert_eq(i_error, 0)
			allocate(point_data%rhs(i_points), stat = i_error); assert_eq(i_error, 0)
		end subroutine

		subroutine point_data_destroy(point_data)
			class(t_output_point_data), intent(inout)		:: point_data

			integer (kind = GRID_SI)						:: i_error

			deallocate(point_data%coords, stat = i_error); assert_eq(i_error, 0)
			deallocate(point_data%p, stat = i_error); assert_eq(i_error, 0)
			deallocate(point_data%S, stat = i_error); assert_eq(i_error, 0)
			deallocate(point_data%rhs, stat = i_error); assert_eq(i_error, 0)
		end subroutine

		subroutine cell_data_create(cell_data, i_cells)
			class(t_output_cell_data), intent(inout)		:: cell_data
			integer(kind = GRID_SI), intent(in)			    :: i_cells

			integer (kind = GRID_SI)						:: i_error

			allocate(cell_data%permeability(i_cells, 3), stat = i_error); assert_eq(i_error, 0)
			allocate(cell_data%porosity(i_cells), stat = i_error); assert_eq(i_error, 0)
			allocate(cell_data%u(i_cells, 3), stat = i_error); assert_eq(i_error, 0)
			allocate(cell_data%rank(i_cells), stat = i_error); assert_eq(i_error, 0)
			allocate(cell_data%section_index(i_cells), stat = i_error); assert_eq(i_error, 0)
			allocate(cell_data%depth(i_cells), stat = i_error); assert_eq(i_error, 0)
			allocate(cell_data%refinement(i_cells), stat = i_error); assert_eq(i_error, 0)
		end subroutine

		subroutine cell_data_destroy(cell_data)
			class(t_output_cell_data), intent(inout)		:: cell_data

			integer (kind = GRID_SI)						:: i_error

			deallocate(cell_data%permeability, stat = i_error); assert_eq(i_error, 0)
			deallocate(cell_data%porosity, stat = i_error); assert_eq(i_error, 0)
			deallocate(cell_data%u, stat = i_error); assert_eq(i_error, 0)
			deallocate(cell_data%rank, stat = i_error); assert_eq(i_error, 0)
			deallocate(cell_data%section_index, stat = i_error); assert_eq(i_error, 0)
			deallocate(cell_data%depth, stat = i_error); assert_eq(i_error, 0)
			deallocate(cell_data%refinement, stat = i_error); assert_eq(i_error, 0)
		end subroutine
	END MODULE
#endif
