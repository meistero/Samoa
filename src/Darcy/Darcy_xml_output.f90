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
			real (kind = GRID_SR), allocatable				    :: coords(:, :)
			real (kind = GRID_SR), allocatable			        :: p(:)
			real (kind = GRID_SR), allocatable			        :: S(:)
			real (kind = GRID_SR), allocatable			        :: rhs(:)

            integer (kind = GRID_SI)                            :: index

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
            integer (kind = GRID_SI), allocatable			    :: connectivity(:)

            integer (kind = GRID_SI)                            :: index

			contains

			procedure, pass :: create => cell_data_create
			procedure, pass :: destroy => cell_data_destroy
		end type

        type num_traversal_data
            type(t_output_point_data), pointer                      :: point_data => null()
            type(t_output_cell_data), pointer                       :: cell_data => null()
            character(len=256)							            :: s_file_stamp
            integer (kind = GRID_SI)								:: i_output_iteration = 0

            contains

            procedure, private, pass :: assign_num_traversal_data

            generic :: assignment(=) => assign_num_traversal_data
        end type

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

#		define	_GT_ELEMENT_OP						element_op

#		define  _GT_NODE_FIRST_TOUCH_OP			    node_first_touch_op

#		include "SFC_generic_traversal_ringbuffer.f90"

		subroutine point_data_create(point_data, i_points)
			class(t_output_point_data), intent(inout)		:: point_data
			integer(kind = GRID_SI), intent(in)			    :: i_points

			integer (kind = GRID_SI)						:: i_error

			allocate(point_data%coords(i_points, 3), stat = i_error); assert_eq(i_error, 0)
			allocate(point_data%p(i_points), stat = i_error); assert_eq(i_error, 0)
			allocate(point_data%S(i_points), stat = i_error); assert_eq(i_error, 0)
			allocate(point_data%rhs(i_points), stat = i_error); assert_eq(i_error, 0)

			point_data%index = 1
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

#           if (_DARCY_LAYERS > 0)
                allocate(cell_data%connectivity(6 * i_cells), stat = i_error); assert_eq(i_error, 0)
#           else
                allocate(cell_data%connectivity(3 * i_cells), stat = i_error); assert_eq(i_error, 0)
#           endif

			cell_data%index = 1
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
			deallocate(cell_data%connectivity, stat = i_error); assert_eq(i_error, 0)
		end subroutine

        subroutine assign_num_traversal_data(data_1, data_2)
            class(num_traversal_data), intent(inout)  :: data_1
            type(num_traversal_data), intent(in)      :: data_2

            assert(associated(data_2%point_data))
            assert(associated(data_2%cell_data))

            data_1%point_data => data_2%point_data
            data_1%cell_data => data_2%cell_data

            data_1%s_file_stamp = data_2%s_file_stamp
            data_1%i_output_iteration = data_2%i_output_iteration
		end subroutine

		subroutine pre_traversal_grid_op(traversal, grid)
			type(t_darcy_xml_output_traversal), intent(inout)		:: traversal
			type(t_grid), intent(inout)							    :: grid

			type(t_grid_info)                                       :: grid_info
            type(t_section_info)         	                        :: section_info
			integer (kind = GRID_SI)								:: i_section, i_cells, i_points, i_error, i

            if (rank_MPI == 0) then
                _log_write(1, '(A, I0, A, I0)') " Darcy: output step: ", traversal%i_output_iteration
            end if

            !we can not use grid%get_info in a sequential OpenMP environment, hence reduce the grid info manually

            grid_info = t_grid_info()

            do i_section = 1, size(grid%sections%elements_alloc)
                section_info = grid%sections%elements_alloc(i_section)%get_info()
                grid_info = grid_info + section_info%t_grid_info
            end do

			i_cells = max(1, _DARCY_LAYERS) * grid_info%i_cells
			i_points = (_DARCY_LAYERS + 1) * (grid_info%i_nodes + sum(grid_info%i_boundary_nodes))

            allocate(traversal%point_data, stat = i_error); assert_eq(i_error, 0)
            allocate(traversal%cell_data, stat = i_error); assert_eq(i_error, 0)

			call traversal%point_data%create(i_points)
			call traversal%cell_data%create(i_cells)

            do i = 1, size(traversal%sections)
                traversal%sections(i)%point_data => traversal%point_data
                traversal%sections(i)%cell_data => traversal%cell_data
            end do

            call scatter(traversal%s_file_stamp, traversal%sections%s_file_stamp)
            call scatter(traversal%i_output_iteration, traversal%sections%i_output_iteration)
		end subroutine

		subroutine post_traversal_grid_op(traversal, grid)
			type(t_darcy_xml_output_traversal), intent(inout)		:: traversal
			type(t_grid), intent(inout)							    :: grid

            integer (kind = GRID_SI)								:: i_error
            call write_vtu_file(traversal, grid)

            if (rank_MPI == 0) then
                call write_pvtu_file(traversal, grid)
            end if

			call traversal%point_data%destroy()
			call traversal%cell_data%destroy()

            deallocate(traversal%point_data, stat = i_error); assert_eq(i_error, 0)
            deallocate(traversal%cell_data, stat = i_error); assert_eq(i_error, 0)

            traversal%i_output_iteration = traversal%i_output_iteration + 1
        end subroutine

        subroutine write_vtu_file(traversal, grid)
			type(t_darcy_xml_output_traversal), intent(inout)		:: traversal
			type(t_grid), intent(inout)							    :: grid

			integer (kind = GRID_SI)								:: i_cells, i_points, i_error, i
            integer (kind = GRID_SI), allocatable                   :: i_offsets(:)
            integer (kind = 1), allocatable                         :: i_types(:)

			character (len = 256)							        :: s_file_name
            type(t_vtk_writer)                                      :: vtk
			integer(4)												:: e_io

            !we can figurte the actual number of cells and points out by using the global counters
			i_cells = traversal%cell_data%index - 1_SI
			i_points = traversal%point_data%index - 1_SI

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

			write (s_file_name, "(A, A, I0, A, I0, A, I0, A)") trim(traversal%s_file_stamp), "_", traversal%i_output_iteration, "_r", rank_MPI, ".vtu"

#           if defined(_QUAD_PRECISION)
#               warning VTK output does not work for quad precision
#           else
                e_io = vtk%VTK_INI_XML('binary', s_file_name, 'UnstructuredGrid')
                    e_io = vtk%VTK_DAT_XML('field', 'OPEN')
                        e_io = vtk%VTK_VAR_XML(1, 'time', [grid%r_time])
                    e_io = vtk%VTK_DAT_XML('field', 'CLOSE')

                    e_io = vtk%VTK_GEO_XML(i_points, i_cells, traversal%point_data%coords(1:i_points, 1), traversal%point_data%coords(1:i_points, 2), traversal%point_data%coords(1:i_points, 3))

                    e_io = vtk%VTK_CON_XML(i_cells, traversal%cell_data%connectivity, i_offsets, i_types)

                    e_io = vtk%VTK_DAT_XML('node', 'OPEN')
                        e_io = vtk%VTK_VAR_XML(i_points, 'pressure', traversal%point_data%p(1:i_points))
                        e_io = vtk%VTK_VAR_XML(i_points, 'saturation', traversal%point_data%S(1:i_points))
                        e_io = vtk%VTK_VAR_XML(i_points, 'rhs',  traversal%point_data%rhs(1:i_points))
                    e_io = vtk%VTK_DAT_XML('node', 'CLOSE')

                    e_io = vtk%VTK_DAT_XML('cell', 'OPEN')
                        e_io = vtk%VTK_VAR_XML(i_cells, 'permeability', traversal%cell_data%permeability(:, 1), traversal%cell_data%permeability(:, 2), traversal%cell_data%permeability(:, 3))
                        e_io = vtk%VTK_VAR_XML(i_cells, 'porosity', traversal%cell_data%porosity)
                        e_io = vtk%VTK_VAR_XML(i_cells, 'flux', traversal%cell_data%u(:, 1), traversal%cell_data%u(:, 2), traversal%cell_data%u(:, 3))
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
        end subroutine

        subroutine write_pvtu_file(traversal, grid)
			type(t_darcy_xml_output_traversal), intent(inout)		:: traversal
			type(t_grid), intent(inout)							    :: grid

			character (len = 256)							        :: s_file_name
            type(t_vtk_writer)                                      :: vtk
			integer(4)												:: e_io

			integer (kind = GRID_SI) :: i_rank

#           if defined(_QUAD_PRECISION)
#               warning VTK output does not work for quad precision
#           else
                write (s_file_name, "(A, A, I0, A, I0, A)") trim(traversal%s_file_stamp), "_", traversal%i_output_iteration, ".pvtu"

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
                        e_io = vtk%VTK_VAR_XML('flux', 1.0_GRID_SR, 3)
                        e_io = vtk%VTK_VAR_XML('rank', 1_GRID_SI, 1)
                        e_io = vtk%VTK_VAR_XML('section index', 1_GRID_SI, 1)
                        e_io = vtk%VTK_VAR_XML('depth', 1_1, 1)
                        e_io = vtk%VTK_VAR_XML('refinement flag', 1_1, 1)
                    e_io = vtk%VTK_DAT_XML('pcell', 'CLOSE')

                    e_io = vtk%VTK_GEO_XML(1.0_GRID_SR)

                    do i_rank = 0, size_MPI - 1
                        write(s_file_name, '(A, "_", I0, "_r", I0, ".vtu")') trim(remove_dir_from_path(traversal%s_file_stamp)), traversal%i_output_iteration, i_rank

                        e_io = vtk%VTK_GEO_XML(s_file_name)
                    end do

                e_io = vtk%VTK_END_XML()
#           endif
		end subroutine

        function remove_dir_from_path(in_file_name) result(out_file_name)
            character(*)            :: in_file_name
            character(len = 256)    :: out_file_name

            out_file_name = trim(in_file_name(scan(in_file_name, "/\", .true.) + 1 :))
        end function

		!******************
		!Geometry operators
		!******************

		subroutine element_op(traversal, section, element)
 			type(t_darcy_xml_output_traversal), intent(inout)	:: traversal
			type(t_grid_section), intent(inout)					:: section
			type(t_element_base), intent(inout)					:: element

#           if (_DARCY_LAYERS > 0)
                real(kind = GRID_SR)        :: p(_DARCY_LAYERS + 1, 3)
                real(kind = GRID_SR)        :: rhs(_DARCY_LAYERS + 1, 3)
                real(kind = GRID_SR)        :: saturation(_DARCY_LAYERS + 1, 3)

                real (kind = GRID_SR)       :: r_point_data_indices(_DARCY_LAYERS + 1, 3)   !< point data indices
                integer (kind = GRID_SI)    :: point_data_indices(_DARCY_LAYERS + 1, 3)		!< point data indices
#           else
                real(kind = GRID_SR)        :: p(3)
                real(kind = GRID_SR)        :: rhs(3)
                real(kind = GRID_SR)        :: saturation(3)

                real (kind = GRID_SR)       :: r_point_data_indices(3)		!< point data indices
                integer (kind = GRID_SI)    :: point_data_indices(3)		!< point data indices
#           endif

			call gv_p%read_from_element(element, p)

			call gv_rhs%read_from_element(element, rhs)
			call gv_r%read_from_element(element, r_point_data_indices)
			call gv_saturation%read_from_element(element, saturation)

			point_data_indices = int(r_point_data_indices, kind=GRID_SI)

            call write_element_data(section%index, element, p, rhs, saturation, element%cell%data_pers%base_permeability, element%cell%data_pers%porosity, traversal%cell_data, traversal%point_data, point_data_indices)
		end subroutine


        subroutine write_element_data(section_index, element, p, rhs, saturation, base_permeability, porosity, cell_data, point_data, point_data_indices)
			type(t_element_base), intent(inout)	        :: element
            type(t_output_cell_data), intent(inout)     :: cell_data
            type(t_output_point_data), intent(inout)    :: point_data
            integer  (kind = GRID_SI), intent(in)       :: section_index

#           if (_DARCY_LAYERS > 0)
                real (kind = GRID_SR), intent(in)	    :: p(:, :)
                real (kind = GRID_SR), intent(in)	    :: rhs(:, :)
                real (kind = GRID_SR), intent(in)	    :: base_permeability(:,:)
                real (kind = GRID_SR), intent(in)	    :: porosity(:)
                real (kind = GRID_SR), intent(inout)	:: saturation(:, :)
                integer (kind = GRID_SI), intent(in)    :: point_data_indices(:, :)

                real (kind = SR)                :: lambda_w(_DARCY_LAYERS + 1, 3), lambda_n(_DARCY_LAYERS + 1, 3)
                real (kind = SR)                :: u_w(_DARCY_LAYERS, 7), u_n(_DARCY_LAYERS, 7)
                real (kind = GRID_SR)			:: g_local(3), flux_w(_DARCY_LAYERS, 3), flux_n(_DARCY_LAYERS, 3), flux_t(3)
#           else
                real (kind = GRID_SR), intent(in)	    :: p(:)
                real (kind = GRID_SR), intent(in)	    :: rhs(:)
                real (kind = GRID_SR), intent(in)	    :: base_permeability
                real (kind = GRID_SR), intent(in)	    :: porosity
                real (kind = GRID_SR), intent(inout)	:: saturation(:)
                integer (kind = GRID_SI), intent(in)    :: point_data_indices(:)

                real (kind = SR)                :: lambda_w(3), lambda_n(3)
                real (kind = SR)                :: u_w(2), u_n(2)
                real (kind = GRID_SR)			:: g_local(2), flux_w(2), flux_n(2), flux_t(2)
#           endif

            integer                         :: i, layer
            real (kind = GRID_SR)			:: edge_length, dz, pos_tmp(3)
            integer  (kind = GRID_SI)       :: i_cell_data_index

            !$omp critical(cell_data_index)
                i_cell_data_index = cell_data%index
                cell_data%index = cell_data%index + max(1_SI, _DARCY_LAYERS)
            !$omp end critical(cell_data_index)

            lambda_w = l_w(saturation)
            lambda_n = l_n(saturation)

            !rotate g so it points in the right direction (no scaling!)

#           if (_DARCY_LAYERS > 0)
                g_local = cfg%g
                g_local(1:2) = samoa_world_to_barycentric_normal(element%transform_data, cfg%g(1:2))
#           else
                g_local(1:2) = samoa_world_to_barycentric_normal(element%transform_data, cfg%g(1:2))
#           endif

            flux_w = 0.0_SR
            flux_n = 0.0_SR

#           if (_DARCY_LAYERS > 0)
                edge_length = element%cell%geometry%get_leg_size()
                dz = cfg%dz

                !compute fluxes

                call compute_base_fluxes_3D(p, base_permeability, edge_length, edge_length, dz, 1.0_SR, 1.0_SR, 1.0_SR, g_local, u_w, u_n)
                call compute_flux_vector_3D(saturation, u_w, u_n, flux_w, flux_n)

                do layer = 1, _DARCY_LAYERS
                    cell_data%rank(i_cell_data_index) = rank_MPI
                    cell_data%section_index(i_cell_data_index) = section_index
                    !convert the permeability back to millidarcy
                    cell_data%permeability(i_cell_data_index, 1:2) = base_permeability(layer, 1) / _MDY
                    cell_data%permeability(i_cell_data_index, 3) = base_permeability(layer, 2) / _MDY

                    !the grid porosity also contains residual saturation which has to be removed for output
                    cell_data%porosity(i_cell_data_index) = porosity(layer) / (1.0_SR - cfg%S_wr - cfg%S_nr)

                    flux_t = 2.0_SR * (flux_w(layer, :) + flux_n(layer, :))
                    flux_t(1:2) = samoa_barycentric_to_world_normal(element%transform_data, flux_t(1:2))

                    cell_data%u(i_cell_data_index, :) = flux_t / (_M / _S)  !return the flux in m/s

                    cell_data%depth(i_cell_data_index) = element%cell%geometry%i_depth
                    cell_data%refinement(i_cell_data_index) = element%cell%geometry%refinement

                    if (element%transform_data%plotter_data%orientation > 0) then
                        cell_data%connectivity(6 * i_cell_data_index - 5) = point_data_indices(layer, 1) - 1
                        cell_data%connectivity(6 * i_cell_data_index - 4) = point_data_indices(layer, 2) - 1
                        cell_data%connectivity(6 * i_cell_data_index - 3) = point_data_indices(layer, 3) - 1
                        cell_data%connectivity(6 * i_cell_data_index - 2) = point_data_indices(layer + 1, 1) - 1
                        cell_data%connectivity(6 * i_cell_data_index - 1) = point_data_indices(layer + 1, 2) - 1
                        cell_data%connectivity(6 * i_cell_data_index - 0) = point_data_indices(layer + 1, 3) - 1
                    else
                        cell_data%connectivity(6 * i_cell_data_index - 5) = point_data_indices(layer, 3) - 1
                        cell_data%connectivity(6 * i_cell_data_index - 4) = point_data_indices(layer, 2) - 1
                        cell_data%connectivity(6 * i_cell_data_index - 3) = point_data_indices(layer, 1) - 1
                        cell_data%connectivity(6 * i_cell_data_index - 2) = point_data_indices(layer + 1, 3) - 1
                        cell_data%connectivity(6 * i_cell_data_index - 1) = point_data_indices(layer + 1, 2) - 1
                        cell_data%connectivity(6 * i_cell_data_index - 0) = point_data_indices(layer + 1, 1) - 1
                    end if

                    i_cell_data_index = i_cell_data_index + 1
                end do

                do layer = 1, _DARCY_LAYERS + 1
                    do i = 1, 3
                        pos_tmp(1:2) = cfg%scaling * samoa_barycentric_to_world_point(element%transform_data, samoa_basis_p_get_dof_coords(i)) + cfg%offset(1:2)
                        pos_tmp(3) = cfg%scaling * real(layer - 1, SR) * cfg%dz + cfg%offset(3)

                        !these are potential concurrent write accesses as other threads are able to access the data at the same time
                        point_data%coords(point_data_indices(layer, i), :) = pos_tmp
                        point_data%p(point_data_indices(layer, i)) = p(layer, i) / _PPSI                !return the pressure in ppsi
                        point_data%rhs(point_data_indices(layer, i)) = rhs(layer, i) / ((_M ** 3) / _S) !return the rhs in m^3 / s
                        point_data%S(point_data_indices(layer, i)) = saturation(layer, i)
                    end do
                end do
#           else
                edge_length = element%cell%geometry%get_leg_size()

                cell_data%rank(i_cell_data_index) = rank_MPI
                cell_data%section_index(i_cell_data_index) = section_index

                !convert the permeability back to millidarcy
                cell_data%permeability(i_cell_data_index, :) = base_permeability / _MDY

                !the grid porosity also contains residual saturation which has to be removed for output
                cell_data%porosity(i_cell_data_index) = porosity / (1.0_SR - cfg%S_wr - cfg%S_nr)

                !compute fluxes

                call compute_base_fluxes_2D(p, base_permeability, edge_length, edge_length, 1.0_SR, 1.0_SR, g_local(1:2), u_w, u_n)
                call compute_flux_vector_2D(saturation, u_w, u_n, flux_w, flux_n)

                flux_t = 2.0_SR * (flux_w + flux_n)
                flux_t = samoa_barycentric_to_world_normal(element%transform_data, flux_t)

                cell_data%u(i_cell_data_index, 1:2) = flux_t / (_M / _S)    !return the flux in m/s
                cell_data%u(i_cell_data_index, 3) = 0.0_SR

                cell_data%depth(i_cell_data_index) = element%cell%geometry%i_depth
                cell_data%refinement(i_cell_data_index) = element%cell%geometry%refinement

                do i = 1, 3
                    pos_tmp(1:2) = cfg%scaling * samoa_barycentric_to_world_point(element%transform_data, samoa_basis_p_get_dof_coords(i)) + cfg%offset(1:2)
                    pos_tmp(3) = cfg%offset(3)

                    !these are potential concurrent write accesses as other threads are able to access the data at the same time
                    point_data%coords(point_data_indices(i), :) = pos_tmp
                    point_data%p(point_data_indices(i)) = p(i) / _PPSI                  !return the pressure in ppsi
                    point_data%rhs(point_data_indices(i)) = rhs(i) / ((_M ** 2) / _S)   !return the rhs in m^2 / s
                    point_data%S(point_data_indices(i)) = saturation(i)
                end do

                if (element%transform_data%plotter_data%orientation > 0) then
                    cell_data%connectivity(3 * i_cell_data_index - 2) = point_data_indices(1) - 1
                    cell_data%connectivity(3 * i_cell_data_index - 1) = point_data_indices(2) - 1
                    cell_data%connectivity(3 * i_cell_data_index - 0) = point_data_indices(3) - 1
                else
                    cell_data%connectivity(3 * i_cell_data_index - 2) = point_data_indices(3) - 1
                    cell_data%connectivity(3 * i_cell_data_index - 1) = point_data_indices(2) - 1
                    cell_data%connectivity(3 * i_cell_data_index - 0) = point_data_indices(1) - 1
                end if
#           endif
		end subroutine

		subroutine node_first_touch_op(traversal, section, node)
 			type(t_darcy_xml_output_traversal), intent(inout)	:: traversal
 			type(t_grid_section), intent(in)				    :: section
			type(t_node_data), intent(inout)				    :: node

            integer (kind = GRID_SI) :: i, i_point_data_index

            !$omp critical(point_data_index)
                i_point_data_index = traversal%point_data%index
                traversal%point_data%index = traversal%point_data%index + _DARCY_LAYERS + 1
            !$omp end critical(point_data_index)

            do i = 1, _DARCY_LAYERS + 1
                call pre_dof_op(i_point_data_index + i - 1, node%data_pers%r(i))
            end do
		end subroutine

		subroutine pre_dof_op(i_point_data_index, r)
 			integer(kind = GRID_SI), intent(in)			    :: i_point_data_index
			real(kind = GRID_SR), intent(out)				:: r

            r = real(i_point_data_index, GRID_SR)
		end subroutine
	END MODULE
#endif
