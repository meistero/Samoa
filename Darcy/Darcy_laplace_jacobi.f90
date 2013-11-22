! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_DARCY)
	MODULE Darcy_laplace_jacobi_traversal
		use SFC_edge_traversal

		use Samoa_darcy

        implicit none

        type num_traversal_data
            real (kind = GRID_SR)				:: r_sq
            real (kind = GRID_SR)				:: alpha
        end type

		type(darcy_gv_p)						:: gv_p
		type(darcy_gv_r_temp)					:: gv_r_temp
		type(darcy_gv_mat_diagonal)				:: gv_mat_diagonal

#		define _GT_NAME							t_darcy_jacobi_traversal

#		if (_DARCY_P_EDGE_SIZE > 0)
#			define _GT_EDGES
#			define _GT_EDGES_TEMP
#		endif

#		define _GT_NODES
#		define _GT_NODES_TEMP

#		define _GT_PRE_TRAVERSAL_GRID_OP		pre_traversal_grid_op
#		define _GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op
#		define	_GT_PRE_TRAVERSAL_OP			pre_traversal_op

#		define _GT_ELEMENT_OP					element_op

#		define _GT_NODE_FIRST_TOUCH_OP		    node_first_touch_op
#		define _GT_NODE_LAST_TOUCH_OP		    node_last_touch_op
#		define _GT_NODE_REDUCE_OP		        node_reduce_op
#		define _GT_INNER_NODE_LAST_TOUCH_OP		inner_node_last_touch_op
#		define _GT_INNER_NODE_REDUCE_OP		    inner_node_reduce_op

#		define _GT_NODE_MERGE_OP		        node_merge_op

#		include "SFC_generic_traversal_ringbuffer.f90"

		subroutine pre_traversal_grid_op(traversal, grid)
  			type(t_darcy_jacobi_traversal), intent(inout)					:: traversal
 			type(t_grid), intent(inout)							        :: grid

			call scatter(traversal%alpha, traversal%children%alpha)
		end subroutine

 		subroutine post_traversal_grid_op(traversal, grid)
  			type(t_darcy_jacobi_traversal), intent(inout)					:: traversal
 			type(t_grid), intent(inout)							        :: grid

			call reduce(traversal%r_sq, traversal%children%r_sq, MPI_SUM, .true.)
		end subroutine

 		subroutine pre_traversal_op(traversal, section)
 			type(t_darcy_jacobi_traversal), intent(inout)				:: traversal
  			type(t_grid_section), intent(inout)				:: section

			traversal%r_sq = 0.0_GRID_SR
		end subroutine

		!*******************************
		!Geometry operators
		!*******************************

		!element

		subroutine element_op(traversal, section, element)
 			type(t_darcy_jacobi_traversal), intent(inout)				:: traversal
 			type(t_grid_section), intent(inout)							:: section
			type(t_element_base), intent(inout), target		:: element

			!local variables

			real(kind = GRID_SR), dimension(_DARCY_P_SIZE)	:: p
			real(kind = GRID_SR), dimension(_DARCY_P_SIZE)	:: r
			real(kind = GRID_SR), dimension(_DARCY_P_SIZE)	:: mat_diagonal

			call gv_p%read(element, p)

			!call element operator
			call alpha_volume_op(section, element, p, r, mat_diagonal, element%cell%data_pers%permeability)

			call gv_r_temp%add(element, r)
			call gv_mat_diagonal%add(element, mat_diagonal)
		end subroutine

		! first touches

		elemental subroutine node_first_touch_op(traversal, section, node)
 			type(t_darcy_jacobi_traversal), intent(in)				:: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(inout)			:: node

			call pre_dof_op(node%data_temp%r, node%data_temp%mat_diagonal)
		end subroutine

		!last touches

		elemental subroutine node_last_touch_op(traversal, section, node)
 			type(t_darcy_jacobi_traversal), intent(in)				:: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(inout)			:: node


			if (node%position(1) > 0.0_GRID_SR .and. node%position(1) < 1.0_GRID_SR) then
				call post_dof_op(traversal, node%data_pers%p, node%data_temp%r, node%data_temp%mat_diagonal)
			end if
		end subroutine

		elemental subroutine inner_node_last_touch_op(traversal, section, node)
 			type(t_darcy_jacobi_traversal), intent(in)				:: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(inout)			:: node

			call post_dof_op(traversal, node%data_pers%p, node%data_temp%r, node%data_temp%mat_diagonal)
		end subroutine

		pure subroutine node_reduce_op(traversal, section, node)
 			type(t_darcy_jacobi_traversal), intent(inout)				:: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(in)			:: node

			integer (kind = GRID_SI)					:: i

			if (node%position(1) > 0.0_GRID_SR .and. node%position(1) < 1.0_GRID_SR) then
				do i = 1, _DARCY_P_NODE_SIZE
					call reduce_dof_op(traversal, node%data_temp%r(i))
				end do
			end if
		end subroutine

		pure subroutine inner_node_reduce_op(traversal, section, node)
 			type(t_darcy_jacobi_traversal), intent(inout)				:: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(in)			:: node

			integer (kind = GRID_SI)					:: i

			do i = 1, _DARCY_P_NODE_SIZE
				call reduce_dof_op(traversal, node%data_temp%r(i))
			end do
		end subroutine

		elemental subroutine node_merge_op(local_node, neighbor_node)
 			type(t_node_data), intent(inout)			    :: local_node
			type(t_node_data), intent(in)				    :: neighbor_node

            local_node%data_temp%r = local_node%data_temp%r + neighbor_node%data_temp%r
            local_node%data_temp%mat_diagonal = local_node%data_temp%mat_diagonal + neighbor_node%data_temp%mat_diagonal
        end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		elemental subroutine pre_dof_op(r, mat_diagonal)
 			real(kind = GRID_SR), intent(out)			:: r
			real(kind = GRID_SR), intent(out)			:: mat_diagonal

			!init temp variables to 0

			r = 0.0_GRID_SR
			mat_diagonal = epsilon(1.0_GRID_SR)
		end subroutine

		pure subroutine alpha_volume_op(section, element, p, r, mat_diagonal, permeability)
 			type(t_grid_section), intent(inout)							:: section
			type(t_element_base), intent(inout)									:: element
			real(kind = GRID_SR), dimension(_DARCY_P_SIZE), intent(in)			:: p
			real(kind = GRID_SR), dimension(_DARCY_P_SIZE), intent(out)			:: r
			real(kind = GRID_SR), dimension(_DARCY_P_SIZE), intent(out)			:: mat_diagonal
			real (kind = GRID_SR), intent(in)									:: permeability

			integer (kind = GRID_SI)											:: i

			!add up matrix diagonal
			forall (i = 1 : _DARCY_P_SIZE)
				mat_diagonal(i) = permeability * section%stiffness_matrix(i, i)
			end forall

			!add up residual
			r = -permeability * MATMUL(section%stiffness_matrix, p)
		end subroutine

		elemental subroutine post_dof_op(traversal, p, r, mat_diagonal)
 			type(t_darcy_jacobi_traversal), intent(in)	    :: traversal
			real(kind = GRID_SR), intent(inout)			:: p
			real(kind = GRID_SR), intent(inout)			:: r
			real(kind = GRID_SR), intent(in)			:: mat_diagonal

			!Jacobi-step

			!add r / diag(A) to the unknown p
			r = r / mat_diagonal
			p = p + traversal%alpha * r
		end subroutine

		pure subroutine reduce_dof_op(traversal, r)
 			type(t_darcy_jacobi_traversal), intent(inout)				:: traversal
			real(kind = GRID_SR), intent(in)			:: r

			traversal%r_sq = traversal%r_sq + (r * r)
		end subroutine
	END MODULE

    MODULE Darcy_laplace_jacobi
		use linear_solver
		use SFC_edge_traversal

		use Darcy_laplace_jacobi_traversal

        implicit none

        type, extends(t_linear_solver)      :: t_darcy_jacobi_solver
            type(t_darcy_jacobi_traversal)  :: jacobi

            contains

            procedure, pass :: solve
        end type

        contains

		!*******************************
		!Module interface
		!*******************************

		!> Solves a poisson equation using a Jacobi solver
		!> \returns		number of iterations performed
		function solve(solver, grid) result(i_iteration)
 			class(t_darcy_jacobi_solver), intent(inout)					:: solver
 			type(t_grid), intent(inout)							        :: grid

			integer (kind = GRID_SI)									:: i_iteration
			real (kind = GRID_SR)										:: r_t1, r_t2
			real (kind = GRID_SR)										:: r_sq_old

            !$omp master
			_log_write(3, '(A, ES14.7)') "  Jacobi solver, max residual error:", grid%r_epsilon
            !$omp end master

			!set step size to some initial value
			solver%jacobi%alpha = 1.0_GRID_SR
			r_sq_old = 1.0_GRID_SR

			do i_iteration = 1, huge(1_GRID_SI)
                !$omp master
				_log_write(3, '(A, I0, A, F7.4, A, ES17.10)') "   i: ", i_iteration, ", alpha: ", solver%jacobi%alpha, ", res: ", sqrt(solver%jacobi%r_sq)
                !$omp end master

				!do a jacobi step
				call solver%jacobi%traverse(grid)

				!adjust step size
				if (solver%jacobi%r_sq > r_sq_old) then
					solver%jacobi%alpha = 0.95_GRID_SR * solver%jacobi%alpha
				else
					solver%jacobi%alpha = solver%jacobi%alpha + 0.0001_GRID_SR
				end if

				if (sqrt(solver%jacobi%r_sq) < grid%r_epsilon * grid%r_p0) then
					exit
				end if

				r_sq_old = solver%jacobi%r_sq
			end do

            !$omp master
			_log_write(3, '(A, T30, I0)') "  Jacobi iterations:", i_iteration
			solver%stats = solver%jacobi%stats
			!$omp end master
		end function
	END MODULE
#endif
