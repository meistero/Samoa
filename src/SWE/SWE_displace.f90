! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_SWE)
	MODULE SWE_Displace
		use SFC_edge_traversal
		use SWE_euler_timestep
		use SWE_initialize_bathymetry

		use Samoa_swe

		implicit none

        type num_traversal_data
            integer (kind = GRID_DI)			:: i_refinements_issued
        end type

		type(t_gv_Q)							:: gv_Q
		type(t_lfs_flux)						:: lfs_flux

#		define _GT_NAME							t_swe_displace_traversal

#		define _GT_EDGES

#		define _GT_ELEMENT_OP					element_op

#		define _GT_CELL_TO_EDGE_OP				cell_to_edge_op

#		define _GT_NODE_MPI_TYPE

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

		!******************
		!Geometry operators
		!******************

		subroutine element_op(traversal, section, element)
			type(t_swe_displace_traversal), intent(inout)		:: traversal
			type(t_grid_section), intent(inout)					:: section
			type(t_element_base), intent(inout)					:: element

			type(t_state)			                            :: Q(_SWE_CELL_SIZE)

			call gv_Q%read(element, Q)

			call alpha_volume_op(traversal, section, element, Q)

			call gv_Q%write(element, Q)
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		subroutine alpha_volume_op(traversal, section, element, Q)
			type(t_swe_displace_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)							:: section
			type(t_element_base), intent(inout)						    :: element
			type(t_state), intent(inout)	                            :: Q(:)

			real (kind = GRID_SR)		                                :: db(_SWE_CELL_SIZE)

			!evaluate initial function values at dof positions and compute DoFs
#           if defined(_ASAGI)
                db = -Q%b + get_bathymetry_at_element(section, element, section%r_time)
                Q%h = Q%h + db
                Q%b = Q%b + db
#           endif

            !no coarsening while the earthquake takes place
			element%cell%geometry%refinement = max(0, element%cell%geometry%refinement)
		end subroutine
	END MODULE
#endif
