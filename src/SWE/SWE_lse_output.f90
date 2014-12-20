! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_SWE)
	MODULE SWE_lse_output
		use LIB_VTK_IO

		use SFC_edge_traversal

		use Samoa_swe

        type num_traversal_data
            real (kind = GRID_SR), allocatable                      :: A(:, :)
            real (kind = GRID_SR), allocatable                      :: qp(:)
            real (kind = GRID_SR), allocatable                      :: rhs(:)
            integer (kind = GRID_SI)								:: i_output_iteration = 0
            integer (kind = GRID_SI)								:: i_point_data_index
        end type

        interface node_first_touch_op
            module procedure node_first_touch_op_scalar
            module procedure node_first_touch_op_array
        end interface

        type(swe_gm_A)										    :: gm_A
        type(swe_gv_qp)										    :: gv_qp
        type(swe_gv_r)										    :: gv_r
        type(swe_gv_rhs)										 :: gv_rhs

#		define	_GT_NAME							t_swe_lse_output_traversal

!#		if (_DARCY_P_EDGE_SIZE > 0)
!#			define _GT_EDGES
!#		endif

#		define	_GT_NODES
#       define  _GT_EDGES
#		define _GT_EDGES_TEMP

#		define	_GT_PRE_TRAVERSAL_GRID_OP			pre_traversal_grid_op
#		define	_GT_POST_TRAVERSAL_GRID_OP			post_traversal_grid_op
#		define	_GT_PRE_TRAVERSAL_OP				pre_traversal_op
#		define	_GT_POST_TRAVERSAL_OP				post_traversal_op

#		define	_GT_ELEMENT_OP						element_op

#		define  _GT_NODE_FIRST_TOUCH_OP			    node_first_touch_op

#		include "SFC_generic_traversal_ringbuffer.f90"

		subroutine pre_traversal_grid_op(traversal, grid)
			type(t_swe_lse_output_traversal), intent(inout)		:: traversal
			type(t_grid), intent(inout)							    :: grid

			_log_write(1, '(A, I0, A, I0)') " SWE: LSE output step: ", traversal%i_output_iteration

            call scatter(traversal%i_output_iteration, traversal%children%i_output_iteration)
		end subroutine

		subroutine post_traversal_grid_op(traversal, grid)
			type(t_swe_lse_output_traversal), intent(inout)		:: traversal
			type(t_grid), intent(inout)							    :: grid


            traversal%i_output_iteration = traversal%i_output_iteration + 1
		end subroutine

		subroutine pre_traversal_op(traversal, section)
 			type(t_swe_lse_output_traversal), intent(inout)			:: traversal
 			type(t_grid_section), intent(inout)							:: section

			type(t_section_info)                                        :: grid_info
			integer (kind = GRID_SI)									:: i_error, i_points

			traversal%i_point_data_index = 1

            grid_info = section%get_info()
			i_points = grid_info%i_nodes + sum(grid_info%i_boundary_nodes)

			allocate(traversal%A(i_points, i_points), stat=i_error); assert_eq(i_error, 0)
			allocate(traversal%qp(i_points), stat=i_error); assert_eq(i_error, 0)
			allocate(traversal%rhs(i_points), stat=i_error); assert_eq(i_error, 0)

			traversal%A(:,:) = 0.0d0
			traversal%qp(:) = 0.0d0
			traversal%rhs(:) = 0.0d0
		end subroutine

		subroutine post_traversal_op(traversal, section)
 			type(t_swe_lse_output_traversal), intent(inout)			:: traversal
			type(t_grid_section), intent(inout)							:: section

			integer (kind = GRID_SI)    :: i_points
            integer                     :: i_error, i, j, fileunit
            real (kind = GRID_SR)       :: nrm
            character(256)              :: filename

            fileunit = 12345 + section%index * size_MPI + rank_MPI
            write(filename, fmt = '("output/lse_r", I0, "_s", I0, "_", I0, ".txt")') rank_MPI, section%index, traversal%i_output_iteration

            open(unit=fileunit, file=filename, status='replace', access='sequential', iostat=i_error); assert_eq(i_error, 0)

			i_points = traversal%i_point_data_index - 1

            write(fileunit, '("A:")')
            do i = 1, i_points
                do j = 1, i_points
                    write(fileunit, '(ES23.15, 1X, $)') traversal%A(j, i)
                end do

                write(fileunit, '()')
            end do
            write(fileunit, '()')

            write(fileunit, '("qp:")')
            do i = 1, i_points
                write(fileunit, '(ES23.15)') traversal%qp(i)
            end do
            write(fileunit, '()')

            write(fileunit, '("rhs:")')
            do i = 1, i_points
                write(fileunit, '(ES23.15)') traversal%rhs(i)
            end do
            write(fileunit, '()')

            !test if |rhs - A p| < epsilon. This might fail even though we did everything correct, since we minimize
            !the preconditioned residual norm |D^{-1} (rhs - A p)| and not the residual norm |rhs - A p| in the linear solver
            !nrm = norm2(matmul(traversal%A, traversal%p) - traversal%rhs)
            !assert_le(nrm, cfg%r_epsilon * 1.0e-7 * cfg%r_p0)

            close(unit=fileunit, iostat=i_error); assert_eq(i_error, 0)

			deallocate(traversal%A, stat=i_error); assert_eq(i_error, 0)
			deallocate(traversal%qp, stat=i_error); assert_eq(i_error, 0)
			deallocate(traversal%rhs, stat=i_error); assert_eq(i_error, 0)
		end subroutine

		!******************
		!Geometry operators
		!******************

		subroutine element_op(traversal, section, element)
 			type(t_swe_lse_output_traversal), intent(inout)	:: traversal
			type(t_grid_section), intent(inout)					:: section
			type(t_element_base), intent(inout)					:: element

			!local variables

			real (kind = GRID_SR)       :: A(3, 3)
			real (kind = GRID_SR)       :: qp(3)

			real (kind = GRID_SR)       :: r_indices(3)		!< point data indices
			integer (kind = GRID_SI)    :: indices(3)		!< point data indices
			integer                     :: i, j

			call gm_A%read(element, A)
			call gv_qp%read(element, qp)
			call gv_r%read(element, r_indices)
            !call gv_rhs%read(element, rhs)

			indices(:) = int(r_indices(:), kind=GRID_SI)
			!qp = samoa_basis_p_dofs_to_values(qp)

            do i = 1, 3
                if (indices(i) .ge. 0) then
                    do j = 1, 3
                        if (indices(j) .ge. 0) then
                            traversal%A(indices(j), indices(i)) = traversal%A(indices(j), indices(i)) + A(j, i)
                        else
                            traversal%rhs(indices(i)) = traversal%rhs(indices(i)) - A(j, i) * qp(j)
                        end if
                    end do

                    traversal%qp(indices(i)) = qp(i)
               end if
            end do
		end subroutine

		subroutine node_first_touch_op_array(traversal, section, nodes)
 			type(t_swe_lse_output_traversal), intent(inout)	:: traversal
 			type(t_grid_section), intent(in)				    :: section
			type(t_node_data), intent(inout)				    :: nodes(:)

			integer (kind = GRID_SI)						    :: i, j

            do j = 1, size(nodes)
                do i = 1, 1
                    call pre_dof_op(traversal%i_point_data_index, nodes(j)%data_pers%r(i), nodes(j)%data_pers%is_dirichlet_boundary(i))
                end do
            end do
		end subroutine

		subroutine node_first_touch_op_scalar(traversal, section, node)
 			type(t_swe_lse_output_traversal), intent(inout)	:: traversal
 			type(t_grid_section), intent(in)				    :: section
			type(t_node_data), intent(inout)				    :: node

			integer (kind = GRID_SI)						    :: i

			do i = 1, 1
				call pre_dof_op(traversal%i_point_data_index, node%data_pers%r(i), .false.)
			end do
		end subroutine

		elemental subroutine pre_dof_op(i_point_data_index, r, is_dirichlet)
 			integer(kind = GRID_SI), intent(inout)			:: i_point_data_index
			real(kind = GRID_SR), intent(out)				:: r
			logical, intent(in)                             :: is_dirichlet

            if (is_dirichlet) then
                r = real(-1, GRID_SR)
            else
                r = real(i_point_data_index, GRID_SR)
                i_point_data_index = i_point_data_index + 1
            end if
		end subroutine
	END MODULE
#endif

