! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_DARCY)
	MODULE Darcy_laplace_cg1
		use SFC_edge_traversal

		use Samoa_darcy
		use Darcy_data_types

		implicit none

        type num_traversal_data
            real (kind = GRID_SR)				:: beta						!< update ratio
            real (kind = GRID_SR)				:: d_A_d					!< d^T * A * d
        end type

		type(darcy_gv_p)						:: gv_p
		type(darcy_gv_r)						:: gv_r
		type(darcy_gv_d)						:: gv_d
		type(darcy_gv_A_d)						:: gv_A_d
		type(darcy_gv_mat_diagonal)				:: gv_mat_diagonal

#		define _GT_NAME							t_cg1_traversal

#		if (_DARCY_P_EDGE_SIZE > 0)
#			define _GT_EDGES
#		endif

#		define _GT_NODES
#		define _GT_NO_COORDS

#		define _GT_PRE_TRAVERSAL_GRID_OP		pre_traversal_grid_op
#		define _GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op
#		define _GT_PRE_TRAVERSAL_OP				pre_traversal_op
#		define _GT_ELEMENT_OP					element_op

#		define _GT_NODE_FIRST_TOUCH_OP			node_first_touch_op
#		define _GT_NODE_REDUCE_OP			    node_reduce_op
#		define _GT_INNER_NODE_REDUCE_OP		    inner_node_reduce_op

#		define _GT_NODE_MERGE_OP		        node_merge_op

#		include "SFC_generic_traversal_ringbuffer.f90"

		subroutine pre_traversal_grid_op(traversal, grid)
  			type(t_cg1_traversal), intent(inout)					:: traversal
 			type(t_grid), intent(inout)							        :: grid

			call scatter(traversal%beta, traversal%children%beta)
		end subroutine

 		subroutine post_traversal_grid_op(traversal, grid)
  			type(t_cg1_traversal), intent(inout)					:: traversal
 			type(t_grid), intent(inout)							        :: grid

			call reduce(traversal%d_A_d, traversal%children%d_A_d, MPI_SUM, .true.)
		end subroutine

 		pure subroutine pre_traversal_op(traversal, section)
  			type(t_cg1_traversal), intent(inout)						:: traversal
  			type(t_grid_section), intent(inout)							:: section

			traversal%d_A_d = 0.0_GRID_SR
		end subroutine

		!*******************************
		!Geometry operators
		!*******************************

		!element

		pure subroutine element_op(traversal, section, element)
  			type(t_cg1_traversal), intent(inout)							:: traversal
 			type(t_grid_section), intent(inout)					:: section
			type(t_element_base), intent(inout), target			:: element

			real(kind = GRID_SR)		:: d(_DARCY_P_SIZE)
			real(kind = GRID_SR)		:: A_d(_DARCY_P_SIZE)

			call gv_d%read(element, d)

			!call element operator
			call alpha_volume_op(section%stiffness_matrix, d, A_d, element%cell%data_pers%permeability)

			call gv_A_d%add(element, A_d)
		end subroutine

		!first touches

		elemental subroutine node_first_touch_op(traversal, section, node)
  			type(t_cg1_traversal), intent(in)							:: traversal
 			type(t_grid_section), intent(in)						:: section
			type(t_node_data), intent(inout)				:: node

			call pre_dof_op(traversal%beta, node%data_pers%r, node%data_pers%d, node%data_pers%A_d)
		end subroutine

		!last touches

		elemental subroutine node_reduce_op(traversal, section, node)
  			type(t_cg1_traversal), intent(inout)							:: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(in)			:: node

			integer (kind = 1)					:: i

			if (.not. any(node%data_temp%is_dirichlet_boundary)) then
				do i = 1, _DARCY_P_NODE_SIZE
					call reduce_dof_op(traversal%d_A_d, node%data_pers%d(i), node%data_pers%A_d(i))
				end do
			end if
		end subroutine

		elemental subroutine inner_node_reduce_op(traversal, section, node)
 			type(t_cg1_traversal), intent(inout)			:: traversal
 			type(t_grid_section), intent(in)		        :: section
			type(t_node_data), intent(in)				    :: node

			integer (kind = 1)					            :: i

            do i = 1, _DARCY_P_NODE_SIZE
                call reduce_dof_op(traversal%d_A_d, node%data_pers%d(i), node%data_pers%A_d(i))
            end do
        end subroutine

		elemental subroutine node_merge_op(local_node, neighbor_node)
 			type(t_node_data), intent(inout)			    :: local_node
			type(t_node_data), intent(in)				    :: neighbor_node

            local_node%data_pers%A_d = local_node%data_pers%A_d + neighbor_node%data_pers%A_d

            !assert_pure(local_node%data_temp%is_dirichlet_boundary .eqv. neighbor_node%data_temp%is_dirichlet_boundary)
            !assert_pure(real(local_node%data_temp%mat_diagonal(1)) .eq. real(neighbor_node%data_temp%mat_diagonal(1)))

            !assert_pure(real(local_node%data_pers%p(1)) .eq. real(neighbor_node%data_pers%p(1)))
            !assert_pure(real(local_node%data_pers%r(1)) .eq. real(neighbor_node%data_pers%r(1)))
            !assert_pure(real(local_node%data_pers%d(1)) .eq. real(neighbor_node%data_pers%d(1)))
        end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		elemental subroutine pre_dof_op(beta, r, d, A_d)
  			real (kind = GRID_SR), intent(in)		:: beta
 			real (kind = GRID_SR), intent(in)		:: r
			real (kind = GRID_SR), intent(inout)	:: d
			real (kind = GRID_SR), intent(out)	    :: A_d

			d = r + beta * d
			A_d = 0.0_GRID_SR
		end subroutine

		pure subroutine alpha_volume_op(A, d, A_d, permeability)
 			real (kind = GRID_SR), intent(in) 	:: A(_DARCY_P_SIZE, _DARCY_P_SIZE)
			real (kind = GRID_SR), intent(in)	:: d(_DARCY_P_SIZE)
			real (kind = GRID_SR), intent(out)	:: A_d(_DARCY_P_SIZE)
			real (kind = GRID_SR), intent(in)	:: permeability

			A_d = permeability * matmul(A, d)
		end subroutine

		elemental subroutine reduce_dof_op(d_A_d, d, A_d)
  			real (kind = GRID_SR), intent(inout)	:: d_A_d
			real (kind = GRID_SR), intent(in)		:: d
			real (kind = GRID_SR), intent(in)	    :: A_d

			d_A_d = d_A_d + (d * A_d)
		end subroutine
	END MODULE

	MODULE Darcy_laplace_cg2
		use SFC_edge_traversal

		use Samoa_darcy

		implicit none

        type num_traversal_data
            real (kind = GRID_SR)					:: alpha					!< step size
            real (kind = GRID_SR)					:: r_C_r					!< r^T * C * r
            real (kind = GRID_SR)					:: r_sq						!< r^2
        end type

		type(darcy_gv_p)						:: gv_p
		type(darcy_gv_r)						:: gv_r
		type(darcy_gv_d)						:: gv_d
		type(darcy_gv_A_d)						:: gv_A_d
		type(darcy_gv_mat_diagonal)				:: gv_mat_diagonal

#		define _GT_NAME							t_cg2_traversal

#		if (_DARCY_P_EDGE_SIZE > 0)
#			define _GT_EDGES
#		endif

#		define _GT_NODES
#		define _GT_NO_COORDS

#		define _GT_PRE_TRAVERSAL_GRID_OP		pre_traversal_grid_op
#		define _GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op
#		define _GT_PRE_TRAVERSAL_OP				pre_traversal_op

#		define _GT_ELEMENT_OP					element_op

#		define _GT_NODE_FIRST_TOUCH_OP			node_first_touch_op
#		define _GT_NODE_LAST_TOUCH_OP			node_last_touch_op
#		define _GT_NODE_REDUCE_OP			    node_reduce_op
#		define _GT_INNER_NODE_LAST_TOUCH_OP		inner_node_last_touch_op
#		define _GT_INNER_NODE_REDUCE_OP		    inner_node_reduce_op

#		define _GT_NODE_MERGE_OP		        node_merge_op

#		include "SFC_generic_traversal_ringbuffer.f90"

		subroutine pre_traversal_grid_op(traversal, grid)
  			type(t_cg2_traversal), intent(inout)					:: traversal
 			type(t_grid), intent(inout)							        :: grid

			call scatter(traversal%alpha, traversal%children%alpha)
		end subroutine

 		subroutine post_traversal_grid_op(traversal, grid)
  			type(t_cg2_traversal), intent(inout)					:: traversal
 			type(t_grid), intent(inout)							        :: grid

			call reduce(traversal%r_C_r, traversal%children%r_C_r, MPI_SUM, .true.)
			call reduce(traversal%r_sq, traversal%children%r_sq, MPI_SUM, .true.)
		end subroutine

 		pure subroutine pre_traversal_op(traversal, section)
  			type(t_cg2_traversal), intent(inout)						:: traversal
			type(t_grid_section), intent(inout)							:: section

			traversal%r_C_r = 0.0_GRID_SR
			traversal%r_sq = 0.0_GRID_SR
		end subroutine

		!*******************************
		!Geometry operators
		!*******************************

		pure subroutine element_op(traversal, section, element)
  			type(t_cg2_traversal), intent(in)							:: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_element_base), intent(inout), target		:: element

			!local variables

			real(kind = GRID_SR)	:: mat_diagonal(_DARCY_P_SIZE)

			!call element operator
			call alpha_volume_op(section%stiffness_matrix, mat_diagonal, element%cell%data_pers%permeability)

			call gv_mat_diagonal%add(element, mat_diagonal)
		end subroutine

		elemental subroutine node_first_touch_op(traversal, section, node)
  			type(t_cg2_traversal), intent(in)							:: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(inout)			:: node

			call pre_dof_op(node%data_temp%mat_diagonal)
		end subroutine

		elemental subroutine node_last_touch_op(traversal, section, node)
  			type(t_cg2_traversal), intent(in)							:: traversal
 			type(t_grid_section), intent(in)					:: section
			type(t_node_data), intent(inout)			:: node

			integer (kind = GRID_SI)					:: i, j

			if (.not. any(node%data_temp%is_dirichlet_boundary)) then
                call post_dof_op(traversal%alpha, node%data_pers%p, node%data_pers%r, node%data_pers%d, node%data_pers%A_d, node%data_temp%mat_diagonal)
			end if
		end subroutine

		elemental subroutine inner_node_last_touch_op(traversal, section, node)
  			type(t_cg2_traversal), intent(in)							:: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(inout)			                :: node

			call post_dof_op(traversal%alpha, node%data_pers%p, node%data_pers%r, node%data_pers%d, node%data_pers%A_d, node%data_temp%mat_diagonal)
		end subroutine

		pure subroutine node_reduce_op(traversal, section, node)
  			type(t_cg2_traversal), intent(inout)		    :: traversal
 			type(t_grid_section), intent(in)			    :: section
			type(t_node_data), intent(in)				    :: node

			integer											:: i

			if (.not. any(node%data_temp%is_dirichlet_boundary)) then
				do i = 1, _DARCY_P_NODE_SIZE
					call reduce_dof_op(traversal%r_C_r, traversal%r_sq, node%data_pers%r(i), node%data_temp%mat_diagonal(i))
				end do
			end if
		end subroutine

		pure subroutine inner_node_reduce_op(traversal, section, node)
  			type(t_cg2_traversal), intent(inout)		    :: traversal
 			type(t_grid_section), intent(in)			    :: section
			type(t_node_data), intent(in)				    :: node

			integer											:: i

			do i = 1, _DARCY_P_NODE_SIZE
                call reduce_dof_op(traversal%r_C_r, traversal%r_sq, node%data_pers%r(i), node%data_temp%mat_diagonal(i))
            end do
		end subroutine

		elemental subroutine node_merge_op(local_node, neighbor_node)
 			type(t_node_data), intent(inout)			    :: local_node
			type(t_node_data), intent(in)				    :: neighbor_node

            local_node%data_temp%mat_diagonal = local_node%data_temp%mat_diagonal + neighbor_node%data_temp%mat_diagonal

            !assert_pure(local_node%data_temp%is_dirichlet_boundary .eqv. neighbor_node%data_temp%is_dirichlet_boundary)

            !assert_pure(real(local_node%data_pers%p(1)) .eq. real(neighbor_node%data_pers%p(1)))
            !assert_pure(real(local_node%data_pers%d(1)) .eq. real(neighbor_node%data_pers%d(1)))
            !assert_pure(real(local_node%data_pers%A_d(1)) .eq. real(neighbor_node%data_pers%A_d(1)))
        end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		elemental subroutine pre_dof_op(mat_diagonal)
 			real (kind = GRID_SR), intent(out)			:: mat_diagonal

			mat_diagonal = tiny(1.0_GRID_SR)
		end subroutine

		pure subroutine alpha_volume_op(A, mat_diagonal, permeability)
 			real (kind = GRID_SR), intent(in) 	:: A(_DARCY_P_SIZE, _DARCY_P_SIZE)
			real (kind = GRID_SR), intent(out)	:: mat_diagonal(_DARCY_P_SIZE)
			real (kind = GRID_SR), intent(in)	:: permeability

			integer (kind = GRID_SI)			:: i

			!add up matrix diagonal
			forall (i = 1 : _DARCY_P_SIZE)
				mat_diagonal(i) = permeability * A(i, i)
			end forall
		end subroutine

		elemental subroutine post_dof_op(alpha, p, r, d, A_d, mat_diagonal)
  			real (kind = GRID_SR), intent(in)			:: alpha
			real (kind = GRID_SR), intent(inout)		:: p
			real (kind = GRID_SR), intent(inout)		:: r
			real (kind = GRID_SR), intent(in)			:: d
			real (kind = GRID_SR), intent(in)			:: A_d
			real (kind = GRID_SR), intent(in)			:: mat_diagonal

			p = p + alpha * d
			r = r - alpha * A_d / mat_diagonal
		end subroutine

		elemental subroutine reduce_dof_op(r_C_r, r_sq, r, mat_diagonal)
  			real (kind = GRID_SR), intent(inout)		:: r_C_r
  			real (kind = GRID_SR), intent(inout)		:: r_sq
			real (kind = GRID_SR), intent(in)		    :: r
			real (kind = GRID_SR), intent(in)			:: mat_diagonal

			r_C_r = r_C_r + (r * mat_diagonal * r)
			r_sq = r_sq + (r * r)
		end subroutine
	END MODULE

	MODULE Darcy_laplace_cg2_exact
		use SFC_edge_traversal

		use Samoa_darcy
		use Darcy_data_types

		implicit none

        type num_traversal_data
            real (kind = GRID_SR)					:: alpha					!< step size
            real (kind = GRID_SR)					:: r_C_r					!< r^T * C * r
            real (kind = GRID_SR)					:: r_sq						!< r^2
        end type

		type(darcy_gv_p)						:: gv_p
		type(darcy_gv_r)						:: gv_r
		type(darcy_gv_d)						:: gv_d
		type(darcy_gv_A_d)						:: gv_A_d
		type(darcy_gv_mat_diagonal)				:: gv_mat_diagonal

#		define _GT_NAME							t_cg2_exact_traversal

#		if (_DARCY_P_EDGE_SIZE > 0)
#			define _GT_EDGES
#		endif

#		define _GT_NODES
#		define _GT_NO_COORDS

#		define _GT_PRE_TRAVERSAL_GRID_OP		pre_traversal_grid_op
#		define _GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op
#		define _GT_PRE_TRAVERSAL_OP				pre_traversal_op

#		define _GT_ELEMENT_OP					element_op

#		define _GT_NODE_FIRST_TOUCH_OP			node_first_touch_op
#		define _GT_NODE_LAST_TOUCH_OP			node_last_touch_op
#		define _GT_NODE_REDUCE_OP			    node_reduce_op
#		define _GT_INNER_NODE_LAST_TOUCH_OP		inner_node_last_touch_op
#		define _GT_INNER_NODE_REDUCE_OP		    inner_node_reduce_op

#		define _GT_NODE_MERGE_OP		        node_merge_op

#		include "SFC_generic_traversal_ringbuffer.f90"

		subroutine pre_traversal_grid_op(traversal, grid)
  			type(t_cg2_exact_traversal), intent(inout)					:: traversal
 			type(t_grid), intent(inout)							        :: grid

			call scatter(traversal%alpha, traversal%children%alpha)
		end subroutine

 		subroutine post_traversal_grid_op(traversal, grid)
  			type(t_cg2_exact_traversal), intent(inout)					:: traversal
 			type(t_grid), intent(inout)							        :: grid

			call reduce(traversal%r_C_r, traversal%children%r_C_r, MPI_SUM, .true.)
			call reduce(traversal%r_sq, traversal%children%r_sq, MPI_SUM, .true.)
		end subroutine

 		subroutine pre_traversal_op(traversal, section)
  			type(t_cg2_exact_traversal), intent(inout)					:: traversal
 			type(t_grid_section), intent(inout)							:: section

			traversal%r_C_r = 0.0_GRID_SR
			traversal%r_sq = 0.0_GRID_SR
		end subroutine

		!*******************************
		!Geometry operators
		!*******************************

		subroutine element_op(traversal, section, element)
  			type(t_cg2_exact_traversal), intent(inout)		:: traversal
 			type(t_grid_section), intent(inout)				:: section
			type(t_element_base), intent(inout), target		:: element

			!local variables

			real(kind = GRID_SR)	:: p(_DARCY_P_SIZE), r(_DARCY_P_SIZE), mat_diagonal(_DARCY_P_SIZE)

			call gv_p%read(element, p)

			!call alpha volume operator
			call alpha_volume_op(section%stiffness_matrix, p, r, mat_diagonal, element%cell%data_pers%permeability)

			call gv_r%add(element, r)
			call gv_mat_diagonal%add(element, mat_diagonal)
		end subroutine

		elemental subroutine node_first_touch_op(traversal, section, node)
  			type(t_cg2_exact_traversal), intent(in)							:: traversal
 			type(t_grid_section), intent(in)						:: section
			type(t_node_data), intent(inout)				:: node

			call pre_dof_op(node%data_pers%r, node%data_temp%mat_diagonal)
		end subroutine

		elemental subroutine node_last_touch_op(traversal, section, node)
  			type(t_cg2_exact_traversal), intent(in)			:: traversal
 			type(t_grid_section), intent(in)				:: section
			type(t_node_data), intent(inout)				:: node

			integer											:: i

			if (.not. any(node%data_temp%is_dirichlet_boundary)) then
                call post_dof_op(traversal%alpha, node%data_pers%p, node%data_pers%r, node%data_pers%d, node%data_pers%A_d, node%data_temp%mat_diagonal)
			else
				node%data_pers%r = 0.0_GRID_SR
			end if
		end subroutine

		elemental subroutine inner_node_last_touch_op(traversal, section, node)
  			type(t_cg2_exact_traversal), intent(in)							:: traversal
 			type(t_grid_section), intent(in)				:: section
			type(t_node_data), intent(inout)				:: node

			call post_dof_op(traversal%alpha, node%data_pers%p, node%data_pers%r, node%data_pers%d, node%data_pers%A_d, node%data_temp%mat_diagonal)
		end subroutine

		pure subroutine node_reduce_op(traversal, section, node)
  			type(t_cg2_exact_traversal), intent(inout)							:: traversal
 			type(t_grid_section), intent(in)						:: section
			type(t_node_data), intent(in)				:: node

			integer											:: i

			if (node%position(1) > 0.0_GRID_SR .and. node%position(1) < 1.0_GRID_SR) then
				do i = 1, _DARCY_P_NODE_SIZE
					call reduce_dof_op(traversal%r_C_r, traversal%r_sq, node%data_pers%r(i), node%data_temp%mat_diagonal(i))
				end do
			end if
		end subroutine

		pure subroutine inner_node_reduce_op(traversal, section, node)
  			type(t_cg2_exact_traversal), intent(inout)							:: traversal
 			type(t_grid_section), intent(in)						:: section
			type(t_node_data), intent(in)				:: node

			integer											:: i

            do i = 1, _DARCY_P_NODE_SIZE
                call reduce_dof_op(traversal%r_C_r, traversal%r_sq, node%data_pers%r(i), node%data_temp%mat_diagonal(i))
            end do
		end subroutine

		elemental subroutine node_merge_op(local_node, neighbor_node)
 			type(t_node_data), intent(inout)			    :: local_node
			type(t_node_data), intent(in)				    :: neighbor_node

            local_node%data_pers%r = local_node%data_pers%r + neighbor_node%data_pers%r
            local_node%data_temp%mat_diagonal = local_node%data_temp%mat_diagonal + neighbor_node%data_temp%mat_diagonal

            assert_pure(real(local_node%data_pers%p(1)) .eq. real(neighbor_node%data_pers%p(1)))
            assert_pure(real(local_node%data_pers%d(1)) .eq. real(neighbor_node%data_pers%d(1)))
            assert_pure(real(local_node%data_pers%A_d(1)) .eq. real(neighbor_node%data_pers%A_d(1)))
        end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		elemental subroutine pre_dof_op(r, mat_diagonal)
 			real (kind = GRID_SR), intent(out)			:: r
			real (kind = GRID_SR), intent(out)			:: mat_diagonal

			r = 0.0_GRID_SR
			mat_diagonal = tiny(1.0_GRID_SR)
		end subroutine

		pure subroutine alpha_volume_op(A, p, r, mat_diagonal, permeability)
 			real (kind = GRID_SR), intent(in) 	:: A(_DARCY_P_SIZE, _DARCY_P_SIZE)
			real (kind = GRID_SR), intent(in)	:: p(_DARCY_P_SIZE)
			real (kind = GRID_SR), intent(out)	:: r(_DARCY_P_SIZE)
			real (kind = GRID_SR), intent(out)	:: mat_diagonal(_DARCY_P_SIZE)
			real (kind = GRID_SR), intent(in)	:: permeability

			integer (kind = GRID_SI)			:: i

			!add up matrix diagonal
			forall (i = 1 : _DARCY_P_SIZE)
				mat_diagonal(i) = permeability * A(i, i)
			end forall

			r = -permeability * matmul(A, p)
		end subroutine

		elemental subroutine post_dof_op(alpha, p, r, d, A_d, mat_diagonal)
  			real (kind = GRID_SR), intent(in)		    :: alpha
			real (kind = GRID_SR), intent(inout)		:: p
			real (kind = GRID_SR), intent(inout)		:: r
			real (kind = GRID_SR), intent(in)			:: d
			real (kind = GRID_SR), intent(in)			:: A_d
			real (kind = GRID_SR), intent(in)			:: mat_diagonal

			p = p + alpha * d
			r = (r - alpha * A_d) / mat_diagonal
		end subroutine

		elemental subroutine reduce_dof_op(r_C_r, r_sq, r, mat_diagonal)
  			real (kind = GRID_SR), intent(inout)		:: r_C_r
 			real (kind = GRID_SR), intent(inout)		:: r_sq
			real (kind = GRID_SR), intent(in)		    :: r
			real (kind = GRID_SR), intent(in)			:: mat_diagonal

			r_C_r = r_C_r + (r * mat_diagonal * r)
			r_sq = r_sq + (r * r)
		end subroutine
	END MODULE

	MODULE Darcy_laplace_cg
		use linear_solver
		use SFC_edge_traversal

		use Darcy_laplace_cg1
		use Darcy_laplace_cg2
		use Darcy_laplace_cg2_exact

        implicit none

        type, extends(t_linear_solver) :: t_darcy_cg_solver
            type(t_cg1_traversal)           :: cg1
            type(t_cg2_traversal)           :: cg2
            type(t_cg2_exact_traversal)     :: cg2_exact

            contains

            procedure, pass :: solve
        end type

		public

		contains

		!> Solves a poisson equation using a CG solver
		!> \returns		number of iterations performed
		function solve(solver, grid) result(i_iteration)
			class(t_darcy_cg_solver), intent(inout)					    :: solver
			type(t_grid), intent(inout)									:: grid

			integer (kind = GRID_SI)									:: i_iteration
			real (kind = GRID_SR)										:: r_sq, d_A_d, r_C_r, r_C_r_old, alpha, beta

            !$omp master
			_log_write(3, '(2X, A, ES14.7)') "CG solver, max residual error:", grid%r_epsilon * grid%r_p0
            !$omp end master

			!set step sizes to 0
			alpha = 0.0_GRID_SR
			beta = 0.0_GRID_SR

			!compute initial residual

            solver%cg2_exact%alpha = alpha
			call solver%cg2_exact%traverse(grid)
            r_sq = solver%cg2_exact%r_sq
            r_C_r = solver%cg2_exact%r_C_r
            _log_write(5, '(4X, A, ES17.10, A, ES17.10)') "r^T r: ", r_sq, " r^T C r: ", r_C_r

			do i_iteration = 0, huge(1_GRID_SI)
                !$omp master
				_log_write(2, '(3X, A, I0, A, F0.10, A, F0.10, A, ES17.10)')  "i: ", i_iteration, ", alpha: ", alpha, ", beta: ", beta, ", res: ", sqrt(r_sq)
                !$omp end master

				if (sqrt(r_sq) < grid%r_epsilon * grid%r_p0) then
					exit
				end if

				!first step: compute search direction d (d_new = r + beta * d_old), the respective residual update A d and the scalar d^T A d
				solver%cg1%beta = beta
				call solver%cg1%traverse(grid)
				d_A_d = solver%cg1%d_A_d
                _log_write(5, '(4X, A, ES17.10)') "d A d: ", d_A_d

				!compute step size alpha = r^T C r / d^T A d
				alpha = r_C_r / d_A_d
				r_C_r_old = r_C_r

				!second step: apply unknowns update (alpha * d) and residual update (alpha * A d)
				!every once in a while, we compute the residual r = b - A x explicitly to limit the numerical error
				if (iand(i_iteration, z'ff') == z'ff') then
                    solver%cg2_exact%alpha = alpha
					call solver%cg2_exact%traverse(grid)
                    r_sq = solver%cg2_exact%r_sq
                    r_C_r = solver%cg2_exact%r_C_r

                    !$omp master
                    if (iand(i_iteration, z'3ff') == z'3ff') then
                        _log_write(1, '(3X, A, I0, A, F0.10, A, F0.10, A, ES17.10)')  "i: ", i_iteration, ", alpha: ", alpha, ", beta: ", beta, ", res: ", sqrt(r_sq)
                    end if
                    !$omp end master
				else
                    solver%cg2%alpha = alpha
					call solver%cg2%traverse(grid)
                    r_sq = solver%cg2%r_sq
                    r_C_r = solver%cg2%r_C_r
				end if

                _log_write(5, '(4X, A, ES17.10, A, ES17.10)') "r^T r: ", r_sq, " r^T C r: ", r_C_r

				!compute beta = r^T C r (new) / r^T C r (old)
				beta = r_C_r / r_C_r_old
			end do

            !$omp master
			_log_write(2, '(2X, A, T24, I0)') "CG iterations:", i_iteration
            solver%stats = solver%cg1%stats + solver%cg2%stats + solver%cg2_exact%stats
            !$omp end master
		end function
	END MODULE
#endif
