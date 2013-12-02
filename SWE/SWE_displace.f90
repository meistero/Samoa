! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_SWE)
	MODULE SWE_Displace
		use SFC_edge_traversal
		use SWE_euler_timestep
		use SWE_initialize

		use Samoa_swe

		implicit none

        type num_traversal_data
            integer (kind = GRID_SI)			:: i_refinements_issued
        end type

		type(t_gv_Q)							:: gv_Q
		type(t_lfs_flux)						:: lfs_flux

#		define _GT_NAME							t_swe_displace_traversal

#		define _GT_EDGES

#		define _GT_ELEMENT_OP					element_op

#		define _GT_CELL_TO_EDGE_OP				cell_to_edge_op

#		include "SFC_generic_traversal_ringbuffer.f90"


		!******************
		!Geometry operators
		!******************

		subroutine element_op(traversal, section, element)
			type(t_swe_displace_traversal), intent(inout)				    :: traversal
			type(t_grid_section), intent(inout)							:: section
			type(t_element_base), intent(inout)					:: element

			type(t_state), dimension(_SWE_CELL_SIZE)			:: Q

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
			type(t_state), dimension(_SWE_CELL_SIZE), intent(out)	    :: Q

			integer (kind = GRID_SI)								    :: i
			real (kind = GRID_SR)		                                :: db

			!evaluate initial function values at dof positions and compute DoFs

			do i = 1, _SWE_CELL_SIZE
                db = -Q(i)%b + get_bathymetry(section, samoa_barycentric_to_world_point(element%transform_data, t_basis_Q_get_dof_coords(i)), section%r_time, element%cell%geometry%i_depth / 2_GRID_SI)
				Q(i)%h = Q(i)%h + db
				Q(i)%b = Q(i)%b + db
			end do

            !no coarsening while the earthquake takes place
			if (element%cell%geometry%i_depth < section%i_max_depth .and. any(Q%h .ne. 0.0)) then
                element%cell%geometry%refinement = 1
            else
                element%cell%geometry%refinement = 0
            end if
		end subroutine
	END MODULE
#endif
