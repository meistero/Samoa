! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_DARCY)
	MODULE Darcy_Adapt
		use SFC_edge_traversal
		use Conformity

		use Samoa_darcy
		use Tools_noise
		use Darcy_initialize_pressure

        implicit none

        type num_traversal_data
			real (kind = GRID_SR), dimension(_DARCY_P_SIZE, 2)					:: p_in
        end type

		type(darcy_gv_p)							:: gv_p
		type(darcy_gv_saturation)					:: gv_saturation
		type(darcy_gv_volume)						:: gv_volume

#		define	_GT_NAME							t_darcy_adaption_traversal

#		if (_DARCY_U_EDGE_SIZE > 0 || _DARCY_FLOW_EDGE_SIZE > 0)
#			define _GT_EDGES
#		endif

#		define	_GT_NODES

#		define	_GT_TRANSFER_OP						transfer_op
#		define	_GT_REFINE_OP						refine_op
#		define	_GT_COARSEN_OP						coarsen_op

#		define _GT_NODE_FIRST_TOUCH_OP				node_first_touch_op
#		define _GT_NODE_LAST_TOUCH_OP				node_last_touch_op
#		define _GT_INNER_NODE_LAST_TOUCH_OP			inner_node_last_touch_op

#		define _GT_NODE_MERGE_OP				    node_merge_op
#		define _GT_EDGE_MPI_TYPE

#		include "SFC_generic_adaptive_traversal.f90"

        subroutine create_edge_mpi_type(mpi_edge_type)
            integer, intent(out)            :: mpi_edge_type

            type(t_edge_data)               :: edge
            integer                         :: blocklengths(2), types(2), disps(2), i_error, extent

#           if defined(_MPI)
                blocklengths(1) = 1
                blocklengths(2) = 1

                disps(1) = 0
                disps(2) = sizeof(edge)

                types(1) = MPI_LB
                types(2) = MPI_UB

                call MPI_Type_struct(2, blocklengths, disps, types, mpi_edge_type, i_error); assert_eq(i_error, 0)
                call MPI_Type_commit(mpi_edge_type, i_error); assert_eq(i_error, 0)

                call MPI_Type_extent(mpi_edge_type, extent, i_error); assert_eq(i_error, 0)
                assert_eq(sizeof(edge), extent)

                call MPI_Type_size(mpi_edge_type, extent, i_error); assert_eq(i_error, 0)
                assert_eq(0, extent)
#           endif
        end subroutine

		!******************
		!Geometry operators
		!******************

		subroutine transfer_op(traversal, grid, src_element, dest_element)
 			type(t_darcy_adaption_traversal), intent(inout)							:: traversal
 			type(t_grid_section), intent(inout)							            :: grid
			type(t_traversal_element), intent(inout)									:: src_element
			type(t_traversal_element), intent(inout)									:: dest_element

			real (kind = GRID_SR), dimension(_DARCY_P_SIZE)								:: p
			real (kind = GRID_SR), dimension(_DARCY_FLOW_SIZE)							:: saturation
			real (kind = GRID_SR), dimension(_DARCY_FLOW_SIZE)							:: volume

			!pressure

			call gv_p%read( src_element%t_element_base, p)
			call gv_p%write( dest_element%t_element_base, p)

			!saturation

			call gv_saturation%read( src_element%t_element_base, saturation)

			volume = src_element%cell%geometry%get_volume() * [0.25_GRID_SR, 0.5_GRID_SR, 0.25_GRID_SR]

			call gv_saturation%add( dest_element%t_element_base, volume * saturation)
			call gv_volume%add( dest_element%t_element_base, volume)

			!permeability

			dest_element%cell%data_pers%base_permeability = src_element%cell%data_pers%base_permeability
		end subroutine

		subroutine refine_op(traversal, grid, src_element, dest_element, refinement_path)
 			type(t_darcy_adaption_traversal), intent(inout)							:: traversal
 			type(t_grid_section), intent(inout)										:: grid
			type(t_traversal_element), intent(inout)								:: src_element
			type(t_traversal_element), intent(inout)								:: dest_element
			integer, dimension(:), intent(in)										:: refinement_path

			real (kind = GRID_SR), dimension(_DARCY_P_SIZE)							:: p_in
			real (kind = GRID_SR), dimension(_DARCY_P_SIZE, 2)						:: p_out
			real (kind = GRID_SR), dimension(_DARCY_FLOW_SIZE)						:: saturation_in
			real (kind = GRID_SR), dimension(_DARCY_FLOW_SIZE, 2)					:: saturation_out
			real (kind = GRID_SR), dimension(_DARCY_FLOW_SIZE)						:: volume
			integer																	:: i

			call gv_p%read( src_element%t_element_base, p_in)
			call gv_saturation%read( src_element%t_element_base, saturation_in)

			volume = src_element%cell%geometry%get_volume() * [0.25_GRID_SR, 0.5_GRID_SR, 0.25_GRID_SR]

			do i = 1, size(refinement_path)
				!pressure
				call samoa_basis_p_split(p_in, p_out(:, 1), p_out(:, 2))
				p_in = p_out(:, refinement_path(i))

				!saturation
				saturation_out(:, 1) = [ saturation_in(1), 0.50_GRID_SR * (saturation_in(1) + saturation_in(2)), saturation_in(2) ]
				saturation_out(:, 2) = [ saturation_in(2), 0.50_GRID_SR * (saturation_in(2) + saturation_in(3)), saturation_in(3) ]

				saturation_in = saturation_out(:, refinement_path(i))
				volume = 0.5_GRID_SR * volume
			end do

			call gv_p%write( dest_element%t_element_base, p_in)

			call gv_saturation%add( dest_element%t_element_base, volume * saturation_in)
			call gv_volume%add( dest_element%t_element_base, volume)

			!permeability

			dest_element%cell%data_pers%base_permeability = get_base_permeability(grid, samoa_barycentric_to_world_point(dest_element%transform_data, [1.0_GRID_SR / 3.0_GRID_SR, 1.0_GRID_SR / 3.0_GRID_SR]), dest_element%cell%geometry%i_depth / 2_GRID_SI)
		end subroutine

		subroutine coarsen_op(traversal, grid, src_element, dest_element, refinement_path)
  			type(t_darcy_adaption_traversal), intent(inout)							    :: traversal
			type(t_grid_section), intent(inout)										    :: grid
			type(t_traversal_element), intent(inout)									:: src_element
			type(t_traversal_element), intent(inout)									:: dest_element
			integer, dimension(:), intent(in)											:: refinement_path

			real (kind = GRID_SR), dimension(_DARCY_FLOW_SIZE, 2)						:: saturation_in
			real (kind = GRID_SR), dimension(_DARCY_P_SIZE)								:: p_out
			real (kind = GRID_SR), dimension(_DARCY_FLOW_SIZE)							:: saturation_out
			real (kind = GRID_SR), dimension(_DARCY_FLOW_SIZE)							:: volume
			integer																		:: i

			!pressure
			i = refinement_path(1)
			call gv_p%read( src_element%t_element_base, traversal%p_in(:, i))

			if (i > 1) then
				call samoa_basis_p_merge(traversal%p_in(:, 1), traversal%p_in(:, 2), p_out)
				call gv_p%write( dest_element%t_element_base, p_out)
				!permeability (compute the average)

				dest_element%cell%data_pers%base_permeability = dest_element%cell%data_pers%base_permeability + 0.5_GRID_SR * src_element%cell%data_pers%base_permeability
			else
				dest_element%cell%data_pers%base_permeability = 0.5_GRID_SR * src_element%cell%data_pers%base_permeability
			end if

			!saturation

			saturation_in = 0.0_GRID_SR
			call gv_saturation%read( src_element%t_element_base, saturation_in(:, i))

			volume = src_element%cell%geometry%get_volume() * [0.25_GRID_SR, 0.5_GRID_SR, 0.25_GRID_SR]

			saturation_out = [ saturation_in(1, 1) + saturation_in(2, 1), &
				0.5_GRID_SR * (saturation_in(3, 1) + saturation_in(2, 1) + saturation_in(1, 2) + saturation_in(2, 2)), &
				saturation_in(3, 2) + saturation_in(2, 2) ]

			call gv_saturation%add( dest_element%t_element_base, volume * saturation_out)
			call gv_volume%add( dest_element%t_element_base, volume)

			dest_element%cell%data_pers%base_permeability = get_base_permeability(grid, samoa_barycentric_to_world_point(dest_element%transform_data, [ 1.0_GRID_SR / 3.0_GRID_SR, 1.0_GRID_SR / 3.0_GRID_SR ]), dest_element%cell%geometry%i_depth / 2_GRID_SI)
		end subroutine

		!first touch ops

		elemental subroutine node_first_touch_op(traversal, grid, node)
 			type(t_darcy_adaption_traversal), intent(in)							:: traversal
 			type(t_grid_section), intent(in)							:: grid
			type(t_node_data), intent(inout)				:: node

			call pre_dof_op(node%data_pers%saturation, node%data_temp%volume, node%data_pers%d, node%data_pers%A_d)
		end subroutine

		!merge ops

		elemental subroutine node_merge_op(local_node, neighbor_node)
 			type(t_node_data), intent(inout)			    :: local_node
			type(t_node_data), intent(in)				    :: neighbor_node

            local_node%data_pers%saturation = local_node%data_pers%saturation + neighbor_node%data_pers%saturation
            local_node%data_temp%volume = local_node%data_temp%volume + neighbor_node%data_temp%volume
        end subroutine

		!last touch ops

		subroutine node_last_touch_op(traversal, section, nodes)
 			type(t_darcy_adaption_traversal), intent(in)	    :: traversal
 			type(t_grid_section), intent(inout)				    :: section
			type(t_node_data), intent(inout)					:: nodes(:)

            integer :: i

            do i = 1, size(nodes)
                call inner_node_last_touch_op(traversal, section, nodes(i))

                if (norm2(nodes(i)%position - cfg%r_pos_prod) < 1.0e-2) then
                    nodes(i)%data_pers%p = cfg%r_p_prod
                    nodes(i)%data_temp%is_dirichlet_boundary = .true.
                else if(norm2(nodes(i)%position - cfg%r_pos_in) < 1.0e-2) then
                    nodes(i)%data_pers%p = cfg%r_p_in
                    nodes(i)%data_pers%saturation = 1.0_GRID_SR
                    nodes(i)%data_temp%is_dirichlet_boundary = .true.
                else
                    nodes(i)%data_temp%is_dirichlet_boundary = .false.
                end if
            end do
		end subroutine

		subroutine inner_node_last_touch_op(traversal, section, node)
 			type(t_darcy_adaption_traversal), intent(in)    :: traversal
 			type(t_grid_section), intent(inout)			    :: section
			type(t_node_data), intent(inout)			    :: node

			call post_dof_op(node%data_pers%saturation, node%data_temp%volume)

			node%data_pers%phi = get_porosity(section, node%position)
		end subroutine
		!*******************************
		!Volume and DoF operators
		!*******************************

		elemental subroutine pre_dof_op(saturation, volume, d, A_d)
			real (kind = GRID_SR), intent(out)		:: saturation
			real (kind = GRID_SR), intent(out)		:: volume
			real (kind = GRID_SR), intent(out)		:: d
			real (kind = GRID_SR), intent(out)		:: A_d

			saturation = 0.0_GRID_SR
			volume = 0.0_GRID_SR
			d = 0.0_GRID_SR
			A_d = 0.0_GRID_SR
		end subroutine

		elemental subroutine post_dof_op(saturation, volume)
 			real (kind = GRID_SR), intent(inout)	:: saturation
			real (kind = GRID_SR), intent(in)		:: volume

			saturation = saturation / volume
		end subroutine
	END MODULE
#endif
