! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_DARCY)
	MODULE Darcy_error_estimate
		use SFC_edge_traversal

		use Samoa_darcy

        type num_traversal_data
            integer (kind = GRID_DI)			:: i_refinements_issued
        end type

        type(darcy_gv_saturation)				:: gv_saturation
        type(darcy_gv_p)						:: gv_p

#		define _GT_NAME							t_darcy_error_estimate_traversal

#		if (_DARCY_FLOW_EDGE_SIZE > 0)
#			define _GT_EDGES
#		endif

#		define _GT_NODES

#		define _GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op
#		define _GT_PRE_TRAVERSAL_OP				pre_traversal_op
#		define _GT_ELEMENT_OP					element_op

#		define _GT_NODE_MPI_TYPE
#		define _GT_EDGE_MPI_TYPE

#		include "SFC_generic_traversal_ringbuffer.f90"

        subroutine create_node_mpi_type(mpi_node_type)
            integer, intent(out)            :: mpi_node_type

#           if defined(_MPI)
                type(t_node_data)                       :: node
                integer                                 :: blocklengths(2), types(2), disps(2), type_size, i_error
                integer (kind = MPI_ADDRESS_KIND)       :: lb, ub

                blocklengths(1) = 1
                blocklengths(2) = 1

                disps(1) = 0
                disps(2) = sizeof(node)

                types(1) = MPI_LB
                types(2) = MPI_UB

                call MPI_Type_struct(2, blocklengths, disps, types, mpi_node_type, i_error); assert_eq(i_error, 0)
                call MPI_Type_commit(mpi_node_type, i_error); assert_eq(i_error, 0)

                call MPI_Type_size(mpi_node_type, type_size, i_error); assert_eq(i_error, 0)
                call MPI_Type_get_extent(mpi_node_type, lb, ub, i_error); assert_eq(i_error, 0)

                assert_eq(0, lb)
                assert_eq(0, type_size)
                assert_eq(sizeof(node), ub)
#           endif
        end subroutine

        subroutine create_edge_mpi_type(mpi_edge_type)
            integer, intent(out)            :: mpi_edge_type

#           if defined(_MPI)
                type(t_edge_data)                       :: edge
                integer                                 :: blocklengths(2), types(2), disps(2), type_size, i_error
                integer (kind = MPI_ADDRESS_KIND)       :: lb, ub

                blocklengths(1) = 1
                blocklengths(2) = 1

                disps(1) = 0
                disps(2) = sizeof(edge)

                types(1) = MPI_LB
                types(2) = MPI_UB

                call MPI_Type_struct(2, blocklengths, disps, types, mpi_edge_type, i_error); assert_eq(i_error, 0)
                call MPI_Type_commit(mpi_edge_type, i_error); assert_eq(i_error, 0)

                call MPI_Type_size(mpi_edge_type, type_size, i_error); assert_eq(i_error, 0)
                call MPI_Type_get_extent(mpi_edge_type, lb, ub, i_error); assert_eq(i_error, 0)

                assert_eq(0, lb)
                assert_eq(0, type_size)
                assert_eq(sizeof(edge), ub)
#           endif
        end subroutine

 		subroutine post_traversal_grid_op(traversal, grid)
  			type(t_darcy_error_estimate_traversal), intent(inout)			:: traversal
 			type(t_grid), intent(inout)							        :: grid

			call reduce(traversal%i_refinements_issued, traversal%children%i_refinements_issued, MPI_SUM, .true.)
		end subroutine

		subroutine pre_traversal_op(traversal, section)
            type(t_darcy_error_estimate_traversal), intent(inout)         :: traversal
  			type(t_grid_section), intent(inout)							:: section

			!this variable will be incremented for each cell with a refinement request
			traversal%i_refinements_issued = 0_GRID_DI
		end subroutine

		!*******************************
		!Geometry operators
		!*******************************

		subroutine element_op(traversal, section, element)
 			type(t_darcy_error_estimate_traversal), intent(inout)	:: traversal
 			type(t_grid_section), intent(inout)						:: section
			type(t_element_base), intent(inout)			            :: element

#           if (_DARCY_LAYERS > 0)
                real (kind = GRID_SR)   :: saturation(_DARCY_LAYERS + 1, 3)
                real (kind = GRID_SR)   :: p(_DARCY_LAYERS + 1, 3)
#           else
                real (kind = GRID_SR)   :: saturation(3)
                real (kind = GRID_SR)   :: p(3)
#           endif

			call gv_saturation%read_from_element(element, saturation)
			call gv_p%read_from_element(element, p)

			!call element operator
			call alpha_volume_op(traversal%i_refinements_issued, element%cell%geometry%i_depth, element%cell%geometry%refinement, saturation, p, element%cell%data_pers%base_permeability)
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		subroutine alpha_volume_op(i_refinements_issued, i_depth, i_refinement, saturation, p, base_permeability)
			integer (kind = GRID_DI), intent(inout)		:: i_refinements_issued
			integer (kind = BYTE), intent(in)			:: i_depth
			integer (kind = BYTE), intent(out)			:: i_refinement

#           if (_DARCY_LAYERS > 0)
!               real (kind = GRID_SR), parameter            :: refinement_threshold = 1.0e2_SR
                real (kind = GRID_SR), parameter            :: refinement_threshold = 1.0e8_SR
#           else
                real (kind = GRID_SR), parameter            :: refinement_threshold = 1.0e1_SR
#           endif

            real (kind = GRID_SR)						    :: r_sat_norm, r_p_norm
            logical 									    :: l_coarsen_p, l_coarsen_sat, l_refine_sat, l_relevant

#           if (_DARCY_LAYERS > 0)
                real (kind = GRID_SR), intent(in)		    :: saturation(:,:)
                real (kind = GRID_SR), intent(in)			:: p(:,:)
                real (kind = GRID_SR), intent(in)			:: base_permeability(:,:)

                r_sat_norm = max(maxval(abs(saturation(:, 3) - saturation(:, 2))), maxval(abs(saturation(:, 1) - saturation(:, 2))))
                r_p_norm = max(maxval(abs(p(:, 3) - p(:, 2))), maxval(abs(p(:, 1) - p(:, 2))))
                l_relevant = any(base_permeability > 0.0_GRID_SR)
#           else
                real (kind = GRID_SR), intent(in)		    :: saturation(:)
                real (kind = GRID_SR), intent(in)			:: p(:)
                real (kind = GRID_SR), intent(in)			:: base_permeability

                r_sat_norm = max(abs(saturation(3) - saturation(2)), abs(saturation(1) - saturation(2)))
                r_p_norm = max(abs(p(3) - p(2)), abs(p(1) - p(2)))
                l_relevant = (base_permeability > 0.0_GRID_SR)
#           endif

            !Criteria:
            !refine a variable q if Delta q > threshold * Delta x / (x_max - x_min) * (q_max - q_min)
            !coarsen a variable q if Delta q < threshold / 2 * Delta x / (x_max - x_min) * (q_max - q_min)

			l_refine_sat = r_sat_norm > refinement_threshold * get_edge_size(cfg%i_max_depth)
			l_coarsen_sat = r_sat_norm < refinement_threshold / 2.0_SR * get_edge_size(cfg%i_max_depth)
			l_coarsen_p = r_p_norm < refinement_threshold / 2.0_SR * get_edge_size(cfg%i_max_depth) * cfg%r_p_prod / 4.0_SR

			!* refine the cell if the saturation becomes too steep
			!* coarsen the cell if pressure and saturation are constant within a cell
			!* leave it as is otherwise

			if (i_depth < cfg%i_max_depth .and. l_relevant .and. (l_refine_sat .or. i_depth < cfg%i_min_depth)) then
				_log_write(5, "(A, T30, A, I0, A, L, A, ES9.2)") "  refinement issued:", "depth ", i_depth, ", sat ", l_refine_sat, ", perm ", base_permeability

				i_refinement = 1
				i_refinements_issued = i_refinements_issued + 1_GRID_DI
			else if ((i_depth > cfg%i_min_depth .and. l_coarsen_p .and. l_coarsen_sat) .or. .not. l_relevant) then
				_log_write(5, "(A, T30, A, I0, A, L, A, L, A, ES9.2)") "  coarsening issued:", "depth ", i_depth, ", p ", l_coarsen_p, ", sat ", l_coarsen_sat, ", perm ", base_permeability

				i_refinement = -1
			else
				_log_write(5, "(A, T30, A, I0)") "  keep:", "depth ", i_depth

				i_refinement = 0
			end if
		end subroutine
	END MODULE
#endif
