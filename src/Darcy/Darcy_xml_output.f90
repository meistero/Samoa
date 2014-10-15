! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_DARCY)
	MODULE Darcy_xml_output
		use LIB_VTK_IO

		use SFC_edge_traversal

		use Samoa_darcy

		!> Output point data
		type t_output_point_data
			real (kind = GRID_SR)				    :: coords(2)
			real (kind = GRID_SR)			        :: p
			real (kind = GRID_SR)			        :: S
		end type

		!> Output cell data
		type t_output_cell_data
			real (kind = GRID_SR)			        :: permeability
			real (kind = GRID_SR)			        :: porosity
			real (kind = GRID_SR)					:: u(2)
            integer (kind = GRID_SI)			    :: rank
            integer (kind = GRID_SI)			    :: section_index
			integer (BYTE)					        :: depth
			integer (BYTE)					        :: refinement
		end type

		logical, parameter											:: l_second_order = (_DARCY_P_ORDER > 1)

        type num_traversal_data
            type(t_output_point_data), allocatable		            :: point_data(:)
            type(t_output_cell_data), allocatable			        :: cell_data(:)
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
        type(darcy_gv_u)										    :: gv_u
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
                        e_io = vtk%VTK_DAT_XML('pnode', 'OPEN')
                            e_io = vtk%VTK_VAR_XML('pressure', 1.0_GRID_SR, 1)
                            e_io = vtk%VTK_VAR_XML('saturation', 1.0_GRID_SR, 1)
                        e_io = vtk%VTK_DAT_XML('pnode', 'CLOSE')

                        e_io = vtk%VTK_DAT_XML('pcell', 'OPEN')
                            e_io = vtk%VTK_VAR_XML('permeability', 1.0_GRID_SR, 1)
                            e_io = vtk%VTK_VAR_XML('porosity', 1.0_GRID_SR, 1)
                            e_io = vtk%VTK_VAR_XML('velocity', 1.0_GRID_SR, 3)
                            e_io = vtk%VTK_VAR_XML('rank', 1_GRID_SI, 1)
                            e_io = vtk%VTK_VAR_XML('section index', 1_GRID_SI, 1)
                            e_io = vtk%VTK_VAR_XML('depth', 1_1, 1)
                            e_io = vtk%VTK_VAR_XML('refinement flag', 1_1, 1)
                        e_io = vtk%VTK_DAT_XML('pcell', 'CLOSE')

                        e_io = vtk%VTK_GEO_XML(1.0_GRID_SR)

                        do i_rank = 0, size_MPI
                            do i_section = 1, omp_get_max_threads() * cfg%i_sections_per_thread * 2
                                write (s_file_name, "(A, A, I0, A, I0, A, I0, A)") trim(traversal%s_file_stamp), "_", traversal%i_output_iteration, "_r", i_rank, "_s", i_section, ".vtu"
                                inquire(file = s_file_name, exist = l_exists)

                                if (l_exists) then
                                    write(s_file_name, "(A)") trim(s_file_name(scan(s_file_name, "/\", .true.) + 1 : len(s_file_name)))
                                    e_io = vtk%VTK_GEO_XML(s_file_name)
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
			i_cells = grid_info%i_cells

			if (l_second_order) then
				i_points = grid_info%i_crossed_edges + grid_info%i_color_edges + grid_info%i_nodes + sum(grid_info%i_boundary_edges + grid_info%i_boundary_nodes)
				allocate(traversal%i_connectivity(6 * i_cells), stat = i_error); assert_eq(i_error, 0)
			else
				i_points = grid_info%i_nodes + sum(grid_info%i_boundary_nodes)
				allocate(traversal%i_connectivity(3 * i_cells), stat = i_error); assert_eq(i_error, 0)
			end if

			allocate(traversal%point_data(i_points), stat = i_error); assert_eq(i_error, 0)
			allocate(traversal%cell_data(i_cells), stat = i_error); assert_eq(i_error, 0)

			traversal%i_cell_data_index = 1
			traversal%i_point_data_index = 1
		end subroutine

		subroutine post_traversal_op(traversal, section)
 			type(t_darcy_xml_output_traversal), intent(inout)			:: traversal
			type(t_grid_section), intent(inout)							:: section

			integer (kind = GRID_SI), dimension(:), allocatable			:: i_offsets
			integer (1), dimension(:), allocatable						:: i_types
			real (kind = GRID_SR), dimension(:), allocatable			:: r_empty
            type(t_vtk_writer)                                          :: vtk

			type(t_section_info)                                        :: grid_info
			integer (kind = GRID_SI)									:: i_error, i_cells, i_points
			integer(4)													:: e_io, i

            grid_info = section%get_info()
			i_cells = grid_info%i_cells

			if (l_second_order) then
				i_points = grid_info%i_crossed_edges + grid_info%i_color_edges + grid_info%i_nodes + sum(grid_info%i_boundary_edges + grid_info%i_boundary_nodes)
			else
				i_points = grid_info%i_nodes + sum(grid_info%i_boundary_nodes)
			end if

			allocate(i_offsets(i_cells), stat = i_error); assert_eq(i_error, 0)
			allocate(i_types(i_cells), stat = i_error); assert_eq(i_error, 0)
			allocate(r_empty(max(i_cells, i_points)), stat = i_error); assert_eq(i_error, 0)

			r_empty(:) = 0.0_GRID_SR

			if (l_second_order) then
				i_types = 22_1

				forall (i = 1 : i_cells)
					i_offsets(i) = 6_GRID_SI * i
				end forall
			else
				i_types = 5_1

				forall (i = 1 : i_cells)
					i_offsets(i) = 3_GRID_SI * i
				end forall
			end if

			write (traversal%s_file_stamp, "(A, A, I0, A, I0, A, I0, A)") TRIM(traversal%s_file_stamp), "_", traversal%i_output_iteration, "_r", rank_MPI, "_s", section%index, ".vtu"

#           if defined(_QUAD_PRECISION)
#               warning VTK output does not work for quad precision
#           else
                e_io = vtk%VTK_INI_XML('binary', traversal%s_file_stamp, 'UnstructuredGrid')
                    e_io = vtk%VTK_GEO_XML(i_points, i_cells, traversal%point_data%coords(1), traversal%point_data%coords(2), r_empty(1:i_points))

                    e_io = vtk%VTK_CON_XML(i_cells, traversal%i_connectivity, i_offsets, i_types)

                    e_io = vtk%VTK_DAT_XML('node', 'OPEN')
                        e_io = vtk%VTK_VAR_XML(i_points, 'pressure', traversal%point_data%p)
                        e_io = vtk%VTK_VAR_XML(i_points, 'saturation', traversal%point_data%S)
                    e_io = vtk%VTK_DAT_XML('node', 'CLOSE')

                    e_io = vtk%VTK_DAT_XML('cell', 'OPEN')
                        e_io = vtk%VTK_VAR_XML(i_cells, 'permeability', traversal%cell_data%permeability)
                        e_io = vtk%VTK_VAR_XML(i_cells, 'porosity', traversal%cell_data%porosity)
                        e_io = vtk%VTK_VAR_XML(i_cells, 'velocity', traversal%cell_data%u(1), traversal%cell_data%u(2), r_empty(1:i_cells))
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
			deallocate(r_empty, stat = i_error); assert_eq(i_error, 0)

			deallocate(traversal%i_connectivity, stat = i_error); assert_eq(i_error, 0)
			deallocate(traversal%cell_data, stat = i_error); assert_eq(i_error, 0)
			deallocate(traversal%point_data, stat = i_error); assert_eq(i_error, 0)

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

			integer (kind = GRID_SI)							:: i

			!local variables

			real (kind = GRID_SR), dimension(_DARCY_P_SIZE)					:: p
			real (kind = GRID_SR), dimension(2, _DARCY_U_SIZE)				:: u
			real (kind = GRID_SR), dimension(_DARCY_FLOW_SIZE)				:: saturation

			real (kind = GRID_SR), dimension(_DARCY_P_SIZE)					:: r_point_data_indices		!< point data indices
			integer (kind = GRID_SI), dimension(_DARCY_P_SIZE)				:: point_data_indices		!< point data indices

			call gv_p%read(element, p)
			call gv_r%read(element, r_point_data_indices)
			call gv_u%read_from_element(element, u)
			call gv_saturation%read(element, saturation)

			point_data_indices(:) = int(r_point_data_indices(:), kind=GRID_SI)
			p = samoa_basis_p_dofs_to_values(p) / 6.89e3_GRID_SR
			u(1, :) = samoa_basis_u_dofs_to_values(u(1, :))
			u(2, :) = samoa_basis_u_dofs_to_values(u(2, :))
			saturation = samoa_basis_flow_dofs_to_values(saturation)

			traversal%cell_data(traversal%i_cell_data_index)%rank = rank_MPI
			traversal%cell_data(traversal%i_cell_data_index)%section_index = section%index
			traversal%cell_data(traversal%i_cell_data_index)%permeability = element%cell%data_pers%permeability / 9.869233e-16_SR
			traversal%cell_data(traversal%i_cell_data_index)%porosity = element%cell%data_pers%porosity
			traversal%cell_data(traversal%i_cell_data_index)%u = cfg%scaling * u(:, 1)
			traversal%cell_data(traversal%i_cell_data_index)%depth = element%cell%geometry%i_depth
			traversal%cell_data(traversal%i_cell_data_index)%refinement = element%cell%geometry%refinement

			if (l_second_order) then
				forall (i = 1 : 6)
					traversal%point_data(point_data_indices(i))%coords = cfg%scaling * samoa_barycentric_to_world_point(element%transform_data, samoa_basis_p_get_dof_coords(i)) + cfg%offset
					traversal%point_data(point_data_indices(i))%p = p(i)
					traversal%point_data(point_data_indices(i))%S = saturation(i)
				end forall

				traversal%i_connectivity(6 * traversal%i_cell_data_index - 5 : 6 * traversal%i_cell_data_index) = point_data_indices([ 1, 2, 3, 6, 4, 5 ]) - 1
			else
				forall (i = 1 : 3)
					traversal%point_data(point_data_indices(i))%coords = cfg%scaling * samoa_barycentric_to_world_point(element%transform_data, samoa_basis_p_get_dof_coords(i)) + cfg%offset
					traversal%point_data(point_data_indices(i))%p = p(i)
					traversal%point_data(point_data_indices(i))%S = saturation(i)
				end forall

				traversal%i_connectivity(3 * traversal%i_cell_data_index - 2 : 3 * traversal%i_cell_data_index) = point_data_indices(1:3) - 1
			end if

			traversal%i_cell_data_index = traversal%i_cell_data_index + 1
		end subroutine

		subroutine node_first_touch_op_array(traversal, section, nodes)
 			type(t_darcy_xml_output_traversal), intent(inout)	:: traversal
 			type(t_grid_section), intent(in)				    :: section
			type(t_node_data), intent(inout)				    :: nodes(:)

			integer (kind = GRID_SI)						    :: i, j

            do j = 1, size(nodes)
                do i = 1, _DARCY_P_NODE_SIZE
                    call pre_dof_op(traversal%i_point_data_index, nodes(j)%data_pers%r(i))
                end do
            end do
		end subroutine

		subroutine node_first_touch_op_scalar(traversal, section, node)
 			type(t_darcy_xml_output_traversal), intent(inout)	:: traversal
 			type(t_grid_section), intent(in)				    :: section
			type(t_node_data), intent(inout)				    :: node

			integer (kind = GRID_SI)						    :: i

			do i = 1, _DARCY_P_NODE_SIZE
				call pre_dof_op(traversal%i_point_data_index, node%data_pers%r(i))
			end do
		end subroutine

		elemental subroutine pre_dof_op(i_point_data_index, r)
 			integer(kind = GRID_SI), intent(inout)			:: i_point_data_index
			real(kind = GRID_SR), intent(out)				:: r

			r = real(i_point_data_index, GRID_SR)
			i_point_data_index = i_point_data_index + 1
		end subroutine
	END MODULE
#endif
