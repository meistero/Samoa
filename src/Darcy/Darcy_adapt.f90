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
        end type

		type(darcy_gv_p)							    :: gv_p
		type(darcy_gv_saturation)					    :: gv_saturation
		type(darcy_gv_volume)						    :: gv_volume
		type(darcy_gv_mat_diagonal)					    :: gv_p_volume
        type(darcy_gv_boundary_condition)               :: gv_boundary_condition

#		define	_GT_NAME							t_darcy_adaption_traversal

#		define	_GT_NODES

#		define	_GT_TRANSFER_OP						transfer_op
#		define	_GT_REFINE_OP						refine_op
#		define	_GT_COARSEN_OP						coarsen_op

#		define _GT_CELL_FIRST_TOUCH_OP				cell_first_touch_op
#		define _GT_NODE_FIRST_TOUCH_OP				node_first_touch_op
#		define _GT_NODE_LAST_TOUCH_OP				inner_node_last_touch_op
#		define _GT_INNER_NODE_LAST_TOUCH_OP			inner_node_last_touch_op

#		define _GT_NODE_MERGE_OP				    node_merge_op
#		define _GT_EDGE_MPI_TYPE

#		include "SFC_generic_adaptive_traversal.f90"

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

		!******************
		!Geometry operators
		!******************

		subroutine transfer_op(traversal, section, src_element, dest_element)
 			type(t_darcy_adaption_traversal), intent(inout)							:: traversal
 			type(t_grid_section), intent(inout)							            :: section
			type(t_traversal_element), intent(inout)									:: src_element
			type(t_traversal_element), intent(inout)									:: dest_element

			real (kind = GRID_SR)   :: p_in(_DARCY_LAYERS + 1, 3), saturation_in(_DARCY_LAYERS + 1, 3), porosity(_DARCY_LAYERS + 1), volume(3)
			integer					:: i, level

            !make sure, the effective volume never turns out to be 0
            porosity = epsilon(1.0_SR)

#           if (_DARCY_LAYERS > 0)
                porosity(1 : _DARCY_LAYERS) = porosity(1 : _DARCY_LAYERS) + 0.5_SR * src_element%cell%data_pers%porosity
                porosity(2 : _DARCY_LAYERS + 1) = porosity(2 : _DARCY_LAYERS + 1) + 0.5_SR * src_element%cell%data_pers%porosity
#           else
                porosity = porosity + src_element%cell%data_pers%porosity
#           endif

			call gv_p%read_from_element( src_element%t_element_base, p_in)
			call gv_saturation%read_from_element( src_element%t_element_base, saturation_in)

            volume(1) = src_element%cell%geometry%get_volume() * 0.25_SR
            volume(2) = src_element%cell%geometry%get_volume() * 0.5_SR
            volume(3) = src_element%cell%geometry%get_volume() * 0.25_SR

			call gv_p%add_to_element( dest_element%t_element_base, [p_in(:, 1) * volume(1), p_in(:, 2) * volume(2), p_in(:, 3) * volume(3)])
			call gv_p_volume%add_to_element( dest_element%t_element_base, [spread(volume(1), 1, _DARCY_LAYERS + 1), spread(volume(2), 1, _DARCY_LAYERS + 1), spread(volume(3), 1, _DARCY_LAYERS + 1)])
			call gv_saturation%add_to_element( dest_element%t_element_base, [saturation_in(:, 1) * volume(1) * porosity, saturation_in(:, 2) * volume(2) * porosity, saturation_in(:, 3) * volume(3) * porosity])
			call gv_volume%add_to_element( dest_element%t_element_base, [volume(1) * porosity, volume(2) * porosity, volume(3) * porosity])

			!permeability

			dest_element%cell%data_pers%base_permeability = src_element%cell%data_pers%base_permeability
			dest_element%cell%data_pers%porosity = src_element%cell%data_pers%porosity

            call set_boundary_conditions(dest_element%t_element_base)
		end subroutine

		subroutine refine_op(traversal, section, src_element, dest_element, refinement_path)
 			type(t_darcy_adaption_traversal), intent(inout)							:: traversal
 			type(t_grid_section), intent(inout)										:: section
			type(t_traversal_element), intent(inout)								:: src_element
			type(t_traversal_element), intent(inout)								:: dest_element
			integer, intent(in)										                :: refinement_path(:)

			real (kind = GRID_SR)   :: p_in(_DARCY_LAYERS + 1, 3), saturation_in(_DARCY_LAYERS + 1, 3), porosity(_DARCY_LAYERS + 1), volume(3), weights(3, 3), r_cells(4)
            integer					:: level, i

            !make sure, the effective volume never turns out to be 0
            porosity = epsilon(1.0_SR)

#           if (_DARCY_LAYERS > 0)
                porosity(1 : _DARCY_LAYERS) = porosity(1 : _DARCY_LAYERS) + 0.5_SR * src_element%cell%data_pers%porosity
                porosity(2 : _DARCY_LAYERS + 1) = porosity(2 : _DARCY_LAYERS + 1) + 0.5_SR * src_element%cell%data_pers%porosity
#           else
                porosity = porosity + src_element%cell%data_pers%porosity
#           endif

			call gv_p%read_from_element( src_element%t_element_base, p_in)
			call gv_saturation%read_from_element( src_element%t_element_base, saturation_in)

			r_cells = [0.0_SR, 0.25_SR, 0.75_SR, 1.0_SR]

            do level = 1, size(refinement_path)
                if (refinement_path(level) > 1) then
                    r_cells = 0.5_SR + r_cells / 2.0_SR
                else
                    r_cells = r_cells / 2.0_SR
                end if
            end do

            weights = transpose(reshape([max(0.0_SR, min(0.25_SR, r_cells(2:4)) - max(0.0_SR, r_cells(1:3))), &
                        max(0.0_SR, min(0.75_SR, r_cells(2:4)) - max(0.25_SR, r_cells(1:3))), &
                        max(0.0_SR, min(1.0_SR, r_cells(2:4)) - max(0.75_SR, r_cells(1:3)))], [3, 3]))

            weights = src_element%cell%geometry%get_volume() * weights

            forall (level = 1:_DARCY_LAYERS + 1, i = 1:3)
                p_in(level, i) = dot_product(weights(:, i), p_in(level, :))
                saturation_in(level, i) = dot_product(weights(:, i), saturation_in(level, :))
            end forall

            volume(1) = sum(weights(:, 1))
            volume(2) = sum(weights(:, 2))
            volume(3) = sum(weights(:, 3))

 			call gv_p%add_to_element(dest_element%t_element_base, p_in)
			call gv_p_volume%add_to_element( dest_element%t_element_base, [spread(volume(1), 1, _DARCY_LAYERS + 1), spread(volume(2), 1, _DARCY_LAYERS + 1), spread(volume(3), 1, _DARCY_LAYERS + 1)])
			call gv_saturation%add_to_element( dest_element%t_element_base, [saturation_in(:, 1) * porosity, saturation_in(:, 2) * porosity, saturation_in(:, 3) * porosity])
			call gv_volume%add_to_element( dest_element%t_element_base, [volume(1) * porosity, volume(2) * porosity, volume(3) * porosity])

			!permeability & porosity

            call get_permeability_and_porosity_at_element(section, dest_element%t_element_base, dest_element%cell%data_pers%base_permeability, dest_element%cell%data_pers%porosity)

            call set_boundary_conditions(dest_element%t_element_base)
        end subroutine

		subroutine coarsen_op(traversal, section, src_element, dest_element, coarsening_path)
  			type(t_darcy_adaption_traversal), intent(inout)							    :: traversal
			type(t_grid_section), intent(inout)										    :: section
			type(t_traversal_element), intent(inout)									:: src_element
			type(t_traversal_element), intent(inout)									:: dest_element
			integer, dimension(:), intent(in)											:: coarsening_path

			real (kind = GRID_SR)   :: p_in(_DARCY_LAYERS + 1, 3), saturation_in(_DARCY_LAYERS + 1, 3), porosity(_DARCY_LAYERS + 1), volume(3), weights(3, 3), r_cells(4)
			integer					:: i, level

            !make sure, the effective volume never reaches 0
            porosity = epsilon(1.0_SR)

#           if (_DARCY_LAYERS > 0)
                porosity(1 : _DARCY_LAYERS) = porosity(1 : _DARCY_LAYERS) + 0.5_SR * src_element%cell%data_pers%porosity
                porosity(2 : _DARCY_LAYERS + 1) = porosity(2 : _DARCY_LAYERS + 1) + 0.5_SR * src_element%cell%data_pers%porosity
#           else
                porosity = porosity + src_element%cell%data_pers%porosity
#           endif

			call gv_p%read_from_element( src_element%t_element_base, p_in)
			call gv_saturation%read_from_element( src_element%t_element_base, saturation_in)

			r_cells = [0.0_SR, 0.25_SR, 0.75_SR, 1.0_SR]

            do level = 1, size(coarsening_path)
                if (coarsening_path(level) > 1) then
                    r_cells = 0.5_SR + r_cells / 2.0_SR
                else
                    r_cells = r_cells / 2.0_SR
                end if
            end do

            weights = transpose(reshape([max(0.0_SR, min(0.25_SR, r_cells(2:4)) - max(0.0_SR, r_cells(1:3))), &
                        max(0.0_SR, min(0.75_SR, r_cells(2:4)) - max(0.25_SR, r_cells(1:3))), &
                        max(0.0_SR, min(1.0_SR, r_cells(2:4)) - max(0.75_SR, r_cells(1:3)))], [3, 3]))

            weights = dest_element%cell%geometry%get_volume() * weights

            forall (level = 1:_DARCY_LAYERS + 1, i = 1:3)
                p_in(level, i) = dot_product(weights(:, i), p_in(level, :))
                saturation_in(level, i) = dot_product(weights(:, i), saturation_in(level, :))
            end forall

            volume(1) = sum(weights(:, 1))
            volume(2) = sum(weights(:, 2))
            volume(3) = sum(weights(:, 3))

 			call gv_p%add_to_element(dest_element%t_element_base, p_in)
			call gv_p_volume%add_to_element( dest_element%t_element_base, [spread(volume(1), 1, _DARCY_LAYERS + 1), spread(volume(2), 1, _DARCY_LAYERS + 1), spread(volume(3), 1, _DARCY_LAYERS + 1)])
			call gv_saturation%add_to_element( dest_element%t_element_base, [saturation_in(:, 1) * porosity, saturation_in(:, 2) * porosity, saturation_in(:, 3) * porosity])
			call gv_volume%add_to_element( dest_element%t_element_base, [volume(1) * porosity, volume(2) * porosity, volume(3) * porosity])

            !permeability and porosity: compute the average

            dest_element%cell%data_pers%base_permeability = mean_invert(mean_transform(dest_element%cell%data_pers%base_permeability) + (0.5_SR ** size(coarsening_path) * mean_transform(src_element%cell%data_pers%base_permeability)))
            dest_element%cell%data_pers%porosity = dest_element%cell%data_pers%porosity + (0.5 ** size(coarsening_path)) * src_element%cell%data_pers%porosity

            call set_boundary_conditions(dest_element%t_element_base)
		end subroutine

		subroutine set_boundary_conditions(element)
			type(t_element_base), intent(inout)    :: element

			real (kind = GRID_SR)   :: pos_well(2), pos_vertex(2, 3)
			integer (kind = SI)     :: boundary_condition(3)
            integer :: i

#           if defined(_ASAGI)
                do i = 1, size(cfg%r_pos_in, 2)
                    pos_well = samoa_world_to_barycentric_point(element%transform_data, cfg%r_pos_in(:, i))

                    if (norm2(pos_well) < 1.0_SR + epsilon(1.0_SR)) then
                        if (pos_well(1) > -epsilon(1.0_SR) .and. pos_well(2) > -epsilon(1.0_SR) .and. 1.0_SR - (pos_well(1) + pos_well(2)) > -epsilon(1.0_SR)) then
                            !injection well:

                            boundary_condition = 0

                            if (pos_well(1) > 0.5_SR) then
                                boundary_condition(1) = i
                            else if (pos_well(2) > 0.5_SR) then
                                boundary_condition(3) = i
                            else
                                boundary_condition(2) = i
                            end if

                            call gv_boundary_condition%add_to_element(element, boundary_condition)
                        end if
                    end if
                end do

                do i = 1, size(cfg%r_pos_prod, 2)
                    pos_well = samoa_world_to_barycentric_point(element%transform_data, cfg%r_pos_prod(:, i))

                    if (norm2(pos_well) < 1.0_SR + epsilon(1.0_SR)) then
                        if (pos_well(1) > -epsilon(1.0_SR) .and. pos_well(2) > -epsilon(1.0_SR) .and. 1.0_SR - (pos_well(1) + pos_well(2)) > -epsilon(1.0_SR)) then
                            !production well:

                            boundary_condition = 0_SI

                            if (pos_well(1) > 0.5_SR) then
                                boundary_condition(1) = -i
                            else if (pos_well(2) > 0.5_SR) then
                                boundary_condition(3) = -i
                            else
                                boundary_condition(2) = -i
                            end if

                            call gv_boundary_condition%add_to_element(element, boundary_condition)
                        end if
                    end if
                end do
#           else
                pos_vertex(:, 1) = element%nodes(1)%ptr%position
                pos_vertex(:, 2) = element%nodes(2)%ptr%position
                pos_vertex(:, 3) = element%nodes(3)%ptr%position

                if (any(pos_vertex(1, :) == 0.0 .or. pos_vertex(1, :) == 1.0_SR)) then
                    boundary_condition = 0_SI

                    where (pos_vertex(1, :) == 0.0_SR)
                        boundary_condition = 1_SI
                    else where (pos_vertex(1, :) == 1.0_SR)
                        boundary_condition = -1_SI
                    end where

                    call gv_boundary_condition%add_to_element(element, boundary_condition)
                end if
#           endif
        end subroutine

		!first touch ops
		elemental subroutine cell_first_touch_op(traversal, section, cell)
 			type(t_darcy_adaption_traversal), intent(in)							:: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_cell_data_ptr), intent(inout)				:: cell

            cell%data_pers%base_permeability = mean_invert(0.0_SR)
			cell%data_pers%porosity = 0.0_SR
		end subroutine

		subroutine node_first_touch_op(traversal, section, node)
 			type(t_darcy_adaption_traversal), intent(in)	:: traversal
 			type(t_grid_section), intent(in)				:: section
			type(t_node_data), intent(inout)				:: node

			call pre_dof_op(node%data_pers%saturation, node%data_temp%volume, node%data_pers%p, node%data_temp%mat_diagonal, node%data_pers%d, node%data_pers%A_d, node%data_pers%rhs)
            call init_boundary_conditions(node%data_pers%boundary_condition)
        end subroutine

		!merge ops

		elemental subroutine node_merge_op(local_node, neighbor_node)
 			type(t_node_data), intent(inout)			    :: local_node
			type(t_node_data), intent(in)				    :: neighbor_node

            local_node%data_pers%p = local_node%data_pers%p + neighbor_node%data_pers%p
            local_node%data_pers%saturation = local_node%data_pers%saturation + neighbor_node%data_pers%saturation
            local_node%data_temp%volume = local_node%data_temp%volume + neighbor_node%data_temp%volume
            local_node%data_temp%mat_diagonal = local_node%data_temp%mat_diagonal + neighbor_node%data_temp%mat_diagonal
			call gv_boundary_condition%add(local_node, neighbor_node%data_pers%boundary_condition)
        end subroutine

		!last touch ops

		subroutine node_last_touch_op(traversal, section, nodes)
 			type(t_darcy_adaption_traversal), intent(in)	    :: traversal
 			type(t_grid_section), intent(inout)				    :: section
			type(t_node_data), intent(inout)					:: nodes(:)

            integer :: i

            do i = 1, size(nodes)
                call inner_node_last_touch_op(traversal, section, nodes(i))
            end do
		end subroutine

		subroutine inner_node_last_touch_op(traversal, section, node)
 			type(t_darcy_adaption_traversal), intent(in)    :: traversal
 			type(t_grid_section), intent(inout)			    :: section
			type(t_node_data), intent(inout)			    :: node

			call post_dof_op(node%data_pers%saturation, node%data_pers%p, node%data_temp%volume, node%data_temp%mat_diagonal)
		end subroutine
		!*******************************
		!Volume and DoF operators
		!*******************************

		elemental subroutine pre_dof_op(saturation, volume, p, mat_diagonal, d, A_d, rhs)
			real (kind = GRID_SR), intent(out)		:: saturation, volume, p, mat_diagonal, d, A_d, rhs

			saturation = 0.0_GRID_SR
			volume = 0.0_GRID_SR
			p = 0.0_GRID_SR
			mat_diagonal = 0.0_GRID_SR
			d = 0.0_GRID_SR
			A_d = 0.0_GRID_SR
			rhs = 0.0_SR
		end subroutine

		elemental subroutine init_boundary_conditions(boundary_condition)
			integer (kind = SI), intent(out)        :: boundary_condition

            boundary_condition = 0_SI
        end subroutine

		elemental subroutine post_dof_op(saturation, p, volume, p_volume)
 			real (kind = GRID_SR), intent(inout)	:: saturation
 			real (kind = GRID_SR), intent(inout)	:: p
			real (kind = GRID_SR), intent(in)		:: volume
			real (kind = GRID_SR), intent(in)	    :: p_volume

            !assert_pure(volume > 0.0_SR)

            saturation = saturation / volume
            p = p / p_volume

            !assert_pure(saturation <= 1.0_SR)
		end subroutine
	END MODULE
#endif
