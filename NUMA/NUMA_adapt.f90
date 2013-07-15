! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_NUMA)
	MODULE NUMA_Adapt
		use SFC_edge_traversal
        use Conformity

		use Samoa_NUMA
		use Tools_noise
		use NUMA_initialize
		use NUMA_euler_timestep

		implicit none

		type num_traversal_data
				type(t_state), dimension(_NUMA_CELL_SIZE, 2)							:: Q_in
		end type

		type(t_gv_Q)							:: gv_Q
		type(t_lfs_flux)						:: lfs_flux

#		define _GT_NAME							t_numa_adaption_traversal

#		define _GT_EDGES

#		define _GT_PRE_TRAVERSAL_OP				pre_traversal_op

#		define _GT_TRANSFER_OP					transfer_op
#		define _GT_REFINE_OP					refine_op
#		define _GT_COARSEN_OP					coarsen_op

#		define _GT_CELL_TO_EDGE_OP				cell_to_edge_op
#		define _GT_CELL_LAST_TOUCH_OP			cell_last_touch_op

#		include "SFC_generic_adaptive_traversal.f90"

		subroutine post_traversal_grid_op(traversal, grid)
			type(t_numa_adaption_traversal), intent(inout)				:: traversal
			type(t_grid), intent(inout)							        :: grid

			call reduce(grid%d_max, grid%sections%elements_alloc%d_max, MPI_MAX, .true.)
		end subroutine

		subroutine pre_traversal_op(traversal, section)
			type(t_numa_adaption_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)							:: section

			section%d_max = 0
		end subroutine

		!******************
		!Adaption operators
		!******************

		subroutine transfer_op(traversal, section, src_element, dest_element)
 			type(t_numa_adaption_traversal), intent(inout)	                :: traversal
			type(t_grid_section), intent(inout)											:: section
			type(t_traversal_element), intent(inout)									:: src_element
			type(t_traversal_element), intent(inout)									:: dest_element

			type(t_state), dimension(_SWE_CELL_SIZE)									:: Q

			call gv_Q%read( src_element%t_element_base, Q)
			call gv_Q%write( dest_element%t_element_base, Q)
		end subroutine

		subroutine refine_op(traversal, section, src_element, dest_element, refinement_path)
 			type(t_numa_adaption_traversal), intent(inout)	                :: traversal
			type(t_grid_section), intent(inout)													:: section
			type(t_traversal_element), intent(inout)									:: src_element
			type(t_traversal_element), intent(inout)									:: dest_element
			integer, dimension(:), intent(in)											:: refinement_path


			type(t_state), dimension(_NUMA_CELL_SIZE)									:: Q_in
			type(t_state), dimension(_NUMA_CELL_SIZE, 2)									:: Q_out

			integer																		:: i
			!state vector

			call gv_Q%read( src_element%t_element_base, Q_in)

			do i = 1, size(refinement_path)
				call t_basis_Q_split(Q_in%h, 	Q_out(:, 1)%h, 		Q_out(:, 2)%h)
				call t_basis_Q_split(Q_in%p(1),	Q_out(:, 1)%p(1),	Q_out(:, 2)%p(1))
				call t_basis_Q_split(Q_in%p(2),	Q_out(:, 1)%p(2),	Q_out(:, 2)%p(2))
				call t_basis_Q_split(Q_in%e,	Q_out(:, 1)%e,	Q_out(:, 2)%e)
				call t_basis_Q_split(Q_in%h_ref, 	Q_out(:, 1)%h_ref, 		Q_out(:, 2)%h_ref)
				call t_basis_Q_split(Q_in%p_ref(1),	Q_out(:, 1)%p_ref(1),	Q_out(:, 2)%p_ref(1))
				call t_basis_Q_split(Q_in%p_ref(2),	Q_out(:, 1)%p_ref(2),	Q_out(:, 2)%p_ref(2))
				call t_basis_Q_split(Q_in%e_ref,	Q_out(:, 1)%e_ref,	Q_out(:, 2)%e_ref)
				Q_in = Q_out(:, refinement_path(i))
			end do

			call gv_Q%write( dest_element%t_element_base, Q_in)
		end subroutine

		subroutine coarsen_op(traversal, section, src_element, dest_element, refinement_path)
 			type(t_numa_adaption_traversal), intent(inout)	                :: traversal
			type(t_grid_section), intent(inout)													:: section
			type(t_traversal_element), intent(inout)									:: src_element
			type(t_traversal_element), intent(inout)									:: dest_element
			integer, dimension(:), intent(in)											:: refinement_path

			type(t_state), dimension(_NUMA_CELL_SIZE)									:: Q_out
			integer																		:: i

			!state vector

			i = refinement_path(1)
			call gv_Q%read( src_element%t_element_base, traversal%Q_in(:, i))

			if (i > 1) then
				call t_basis_Q_merge(traversal%Q_in(:, 1)%h,		traversal%Q_in(:, 2)%h,		Q_out%h)
				call t_basis_Q_merge(traversal%Q_in(:, 1)%p(1),	traversal%Q_in(:, 2)%p(1),	Q_out%p(1))
				call t_basis_Q_merge(traversal%Q_in(:, 1)%p(2),	traversal%Q_in(:, 2)%p(2),	Q_out%p(2))
				call t_basis_Q_merge(traversal%Q_in(:, 1)%e,		traversal%Q_in(:, 2)%e,		Q_out%e)
				call t_basis_Q_merge(traversal%Q_in(:, 1)%h_ref,		traversal%Q_in(:, 2)%h_ref,		Q_out%h_ref)
				call t_basis_Q_merge(traversal%Q_in(:, 1)%p_ref(1),	traversal%Q_in(:, 2)%p_ref(1),	Q_out%p_ref(1))
				call t_basis_Q_merge(traversal%Q_in(:, 1)%p_ref(2),	traversal%Q_in(:, 2)%p_ref(2),	Q_out%p_ref(2))
				call t_basis_Q_merge(traversal%Q_in(:, 1)%e_ref,		traversal%Q_in(:, 2)%e_ref,		Q_out%e_ref)

				call gv_Q%write( dest_element%t_element_base, Q_out)
			end if
		end subroutine

			subroutine cell_last_touch_op(traversal, section, cell)
 			type(t_numa_adaption_traversal), intent(inout)	                :: traversal
			type(t_grid_section), intent(inout)							:: section
			type(t_cell_data_ptr), intent(inout)				:: cell
			!set maximum depth
			section%d_max = max(section%d_max, cell%geometry%i_depth)
		end subroutine
	END MODULE
#endif
