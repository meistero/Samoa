! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


!*****************************************************************
!
! AUTHOR: Csaba Vigh (THE BOSS, THE CHAMPION, THE HERO)
!
! function:
!	Defines data structures for FINE LEVEL Triangle cells.
!	A FINE LEVEL trianle contains data that is needed by any
!	traversal routine on the leaf level.
!
!*****************************************************************
#include "Compilation_control.f90"
MODULE SFC_fine_grid
use SFC_grid

implicit none

PUBLIC

#if defined(_ASSERT)
	INTERFACE assignment(=)
		module procedure assign_fem_tri_from_fine
	end INTERFACE
#endif

INTERFACE assignment(=)
	module procedure assign_fine_tri_from_fem
end INTERFACE

!*****************************************************************

CONTAINS

function is_color_boundary(cell)
	type (fine_triangle), intent(in)			:: cell
	logical (kind = GRID_SL)					:: is_color_boundary

	is_color_boundary = (cell%get_color_edge_type() == OLD_BND .or. cell%get_color_edge_type() == NEW_BND)
end function is_color_boundary


! simple takeover of values: FINE = FEM_TRI
subroutine assign_fine_tri_from_fem(fine_tri, fem)
	type(fine_triangle), intent(out)				:: fine_tri
	type(fem_triangle), intent(in)				:: fem
	! local variables
	integer (kind = 1)							:: previous_edge_index, color_edge_index, next_edge_index

	fine_tri%i_plotter_type = fem%i_plotter_type
	fine_tri%i_turtle_type = 4 - abs(fem%i_turtle_type)

	previous_edge_index = get_previous_edge_index(fem%i_turtle_type)
	color_edge_index = get_color_edge_index(fem%i_turtle_type)
	next_edge_index = get_next_edge_index(fem%i_turtle_type)

	if (is_domain_boundary_edge(fem, previous_edge_index)) then
		call fine_tri%set_previous_edge_type(OLD_BND)
	else
		call fine_tri%set_previous_edge_type(OLD)
	endif

	if (is_domain_boundary_edge(fem, color_edge_index) .and. fem%i_turtle_type < 0) then
        call fine_tri%set_color_edge_type(OLD_BND)
    else if (is_domain_boundary_edge(fem, color_edge_index) .and. fem%i_turtle_type > 0) then
        call fine_tri%set_color_edge_type(NEW_BND)
	else if (fem%i_turtle_type < 0) then
		call fine_tri%set_color_edge_type(OLD)
	else
		call fine_tri%set_color_edge_type(NEW)
	endif

	if (is_domain_boundary_edge(fem, next_edge_index)) then
		call fine_tri%set_next_edge_type(NEW_BND)
	else
		call fine_tri%set_next_edge_type(NEW)
	endif

	if (fine_tri%i_turtle_type == V) then
		fine_tri%l_color_edge_color = .not. fem%l_color
	else
		fine_tri%l_color_edge_color = fem%l_color
	endif

	fine_tri%i_depth = fem%i_depth
end subroutine assign_fine_tri_from_fem


# if defined(_ASSERT)
	! rebuilding the FEM_TRIANGLE = FINE
	subroutine assign_fem_tri_from_fine(fem_tri, fine)
		type(fem_triangle), intent(out)				:: fem_tri
		type(fine_triangle), intent(in)				:: fine

		fem_tri%triangle%i_triangle = TRIANGLE_LEAF		! always a LEAF triangle

		!fem_tri%r_leg_size = fine%r_leg_size

		!fem_tri%r_coords(:,:) = fine%r_coords(:,:, l_forward)
		!fem_tri%l_boundary_data(:,:) = fine%l_boundary_data(:,:, l_forward)

		fem_tri%i_plotter_type = fine%i_plotter_type
		fem_tri%i_turtle_type = fine%i_turtle_type

		if (fine%i_turtle_type == V) then
			fem_tri%l_color = .not. fine%l_color_edge_color
		else
			fem_tri%l_color = fine%l_color_edge_color
		endif
		fem_tri%i_depth = fine%i_depth

		!fem_tri%is_n_th_child = ...
	end subroutine assign_fem_tri_from_fine
# endif

end MODULE SFC_fine_grid
