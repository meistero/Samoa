! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


!*****************************************************************
!
! AUTHOR: Csaba Vigh
!
! function:
!	Defines data structures for Triangle cells in form of a binary tree.
!
!	Contains basic routines for generating Adaptive Triangular Grids.
!
! COMMENTS:
!	use different affixes for different data types:
!	  i*_:	integer (kind = GRID_SI) variables
!	  r*_:	real (kind = GRID_SR) variables
!	  p*_:	pointer variables
!	  l*_:	logical variables
!
!*****************************************************************
#include "Compilation_control.f90"

MODULE SFC_grid
	!use ifLPORT
	use SFC_data_types

	implicit none

	PUBLIC

	! a tree representation with temporary stacks provided for traversal
	type triangle_tree
		! root triangle data: coords(:,1) and coords(:,3) are the hypotenuse
		real (kind = GRID_SR), dimension(2,3)	:: r_coords			! 2D coordinates of triangle

#if defined(_ASSERT)
		! size of a leg (LEFT_LEG or RIGHT_LEG): r_leg_size = sqrt(  [r_coords(:, 1) - r_coords(:, 2)]^2  )
		real (kind = GRID_SR)					:: r_leg_size
#endif
		! edge boundary data: boundary_data(:,2) is the hypotenuse
		logical (kind = GRID_SL), dimension(2,3)	:: l_boundary_data	! 00 - internal, 01 - Dirichlet, 10 - Neumann, 11 - Process Boundary

		! type for SFC algorithms, indexed by 'l_forward
		integer (kind = 1)                  		:: i_turtle_type		! k1n =  1, k2n =  2, k3n =  3,
																	! k1o = -1, k2o = -2, k3o = -3

		logical (kind = GRID_SL)					:: l_color			! reference color
		! depth level
		integer (kind = 1)				:: i_depth		! grid level


		integer (kind = 1)						:: i_max_depth
		integer (kind = GRID_SI)				:: i_num_elements		! maintenance needed at creation and adaptation
		integer (kind = GRID_SI)				:: i_num_fine_elements	! maintenance needed at creation and adaptation

		! MPI parallelization indexes for FINE ELEMENTS (inclusive)
		integer (kind = GRID_SI)				:: i_my_fine_elements_start		! cells belonging to this process start here
		integer (kind = GRID_SI)				:: i_my_fine_elements_end		! cells belonging to this process end here
		integer (kind = GRID_SI)				:: i_fine_elements_index		! current cell being processed
		integer (kind = GRID_SI)				:: global_region_index_start	! add this to 'i_fine_elements_index'

		! ESTIMATION: Memory needed for (additional - removed) triangles
		! additional is counted before an adaptation
		! removed was counted after the previous adaptation (that much memory is still available)
		integer (kind=GRID_SI)					:: i_count_additional_triangles	! counted "a priori" before adaptation
		integer (kind=GRID_SI)					:: i_count_removed_triangles	! counted "a posteriori" after adaptation



		!type(triangle), dimension(:), allocatable	:: triangles1, triangles2
		! indices for traversal
		type(triangle), dimension(:), pointer	:: p_input_triangles, p_output_triangles
		integer (kind = GRID_SI)				:: i_input_index, i_output_index
		integer (kind = GRID_SI)				:: i_input_min, i_input_max
		integer (kind = GRID_SI)				:: i_output_min, i_output_max

		!Omer
		integer (kind = GRID_SI)				:: i_matrix_type
	end type triangle_tree

!*****************************************************************

	type(triangle_tree), dimension(:), allocatable, target	:: tree					! 4 initial triangles allocated during execution
	integer (kind = GRID_SI)					:: num_initial_triangles

	!! SFC Triangle cells (integer) ---------- contains also BOUNDARY INFO for the 'color_edge'
	integer (kind = GRID_SI), parameter			:: TRIANGLE_NODE = -1			! this triangle is a BINARY TREE NODE, has children
	integer (kind = GRID_SI), parameter			:: TRIANGLE_DIRICHLET = -20		! Leaf, DIRICHLET BOUNDARY (geometric) finest level triangle
	integer (kind = GRID_SI), parameter			:: TRIANGLE_NEUMANN = -30		! Leaf, NEUMANN BOUNDARY (geometric) finest level triangle
	integer (kind = GRID_SI), parameter			:: TRIANGLE_LEAF = -10			! this triangle is a leaf, finest level triangle

	! ToDo: DEPRECATED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Boundary types
	logical (kind = GRID_SL), dimension(2), parameter		:: internal_edge		 = (/ .false., .false. /)
	logical (kind = GRID_SL), dimension(2), parameter		:: dirichlet_boundary	 = (/ .false.,  .true. /)
	logical (kind = GRID_SL), dimension(2), parameter		:: neumann_boundary		 = (/  .true., .false. /)
	logical (kind = GRID_SL), dimension(2), parameter		:: process_boundary		 = (/  .true.,  .true. /)

	CONTAINS

subroutine initialize_triangle_constants(first_leg_size)
	implicit none
	real (kind = GRID_SR), intent(in)	:: first_leg_size
	! local variables
	integer					:: i

end subroutine initialize_triangle_constants

function get_region_number(global_index)
	integer (kind = GRID_SI), intent(in)		:: global_index
	integer (kind = GRID_SI)					:: get_region_number
	integer (kind = GRID_SI)					:: i

	! search the biggest 'i' for which 'sfc_region_start(i) <= global_index'
	get_region_number = -1
end function get_region_number

! return tree number to which FINE_TRIANGLE(global_index) belongs
function get_tree_number(global_index)
	integer (kind = GRID_SI), intent(in)		:: global_index
	integer (kind = GRID_SI)					:: get_tree_number
	integer (kind = GRID_SI)					:: i

	! search the biggest 'i' for which 'tree(i)%global_region_index_start <= global_index'
	get_tree_number = -1

	do i = 1, num_initial_triangles
		if (tree(i)%global_region_index_start < global_index) then
			get_tree_number = i
		else
			exit
		endif
	end do
end function get_tree_number

function is_in_range(p_tree)
	type (triangle_tree), intent(in)			:: p_tree
	logical										:: is_in_range

	is_in_range =  (p_tree%i_my_fine_elements_start <= p_tree%i_fine_elements_index .and. &
					p_tree%i_fine_elements_index <= p_tree%i_my_fine_elements_end)
end function is_in_range

function is_on_range_boundary(p_tree)
	type (triangle_tree), intent(in)			:: p_tree
	logical										:: is_on_range_boundary
	is_on_range_boundary =  (p_tree%i_my_fine_elements_start == p_tree%i_fine_elements_index .or. &
							p_tree%i_fine_elements_index == p_tree%i_my_fine_elements_end)
end function is_on_range_boundary


function get_num_fine_elements(triangles, num_coarse_triangles)
	type(triangle_tree), dimension(:), intent(in)	:: triangles
	integer (kind = GRID_SI), intent(in)		:: num_coarse_triangles
	integer (kind = GRID_SI)					:: get_num_fine_elements
	! locals
	integer		:: i

	get_num_fine_elements = 0
	do i = 1, num_coarse_triangles
		get_num_fine_elements = get_num_fine_elements + triangles(i)%i_num_fine_elements
	end do
end function get_num_fine_elements

function get_num_elements(triangles, num_coarse_triangles)
	type(triangle_tree), dimension(:), intent(in)	:: triangles
	integer (kind = GRID_SI), intent(in)		:: num_coarse_triangles
	integer (kind = GRID_SI)					:: get_num_elements
	! locals
	integer		:: i

	get_num_elements = 0
	do i = 1, num_coarse_triangles
		get_num_elements = get_num_elements + triangles(i)%i_num_elements
	end do
end function get_num_elements

subroutine print_memory_usage()
	! memory related
	integer				:: m_size, m_resident, m_share, m_text, m_library, m_data_stack

	call read_memory_usage(m_size, m_resident, m_share, m_text, m_library, m_data_stack)
!            unsigned size; //       total program size
!            unsigned resident;//   resident set size
!            unsigned share;//      shared pages
!            unsigned text;//       text (code)
!            unsigned lib;//        library
!            unsigned data;//       data/stack
	PRINT *, "size 			= ", m_size
	PRINT *, "resident 		= ", m_resident
	PRINT *, "share 		= ", m_share
	PRINT *, "text 			= ", m_text
	PRINT *, "library 		= ", m_library
	PRINT *, "data/stack 	= ", m_data_stack
end subroutine print_memory_usage

subroutine read_memory_usage(m_size, m_resident, m_share, m_text, m_library, m_data_stack)
	integer, intent(out)					:: m_size, m_resident, m_share, m_text, m_library, m_data_stack
	! local variables
	character (len = 80)					:: file_name
	integer									:: file_unit = 10101
	integer									:: file_ios
	integer									:: pid

	pid = 0 !getpid()
!	PRINT *, "Process ID = ", pid

	write (file_name, '(a, i0, a)') '/proc/', pid, '/statm'

!	PRINT *, "File name = ", file_name
	open(unit = file_unit, file = file_name, action= 'read', status='old', iostat=file_ios)
	if (file_ios /= 0) then
		PRINT *, "ERROR: '", file_name, "' could not be opened!"
		PRINT *, "IOSTAT = ", file_ios
	else
		read (unit = file_unit, fmt=*) m_size, m_resident, m_share, m_text, m_library, m_data_stack
	end if
	close(unit=file_unit)
end subroutine read_memory_usage

subroutine PRINT_cell(fem_tri)
	type (fem_triangle), intent(in)			:: fem_tri

end subroutine PRINT_cell

function is_leaf_triangle(fem_tri)
	type (fem_triangle), intent(in)			:: fem_tri
	logical									:: is_leaf_triangle

	is_leaf_triangle = (fem_tri%triangle%i_triangle /= TRIANGLE_NODE)
end function is_leaf_triangle


function is_node_triangle(fem_tri)
	type (fem_triangle), intent(in)			:: fem_tri
	logical									:: is_node_triangle

	is_node_triangle = (fem_tri%triangle%i_triangle == TRIANGLE_NODE)
end function is_node_triangle

!
! Reads the first fem_triangle of a tree.
subroutine read_first_fem_triangle(p_triangle_tree, first_fem_triangle)
	type (triangle_tree), intent(inout), target			:: p_triangle_tree 		! intent(inout) (for the contents of the target)
	type (fem_triangle), intent(out)					:: first_fem_triangle

	type (triangle), dimension(:), pointer				:: p_triangle_tree_tmp

	p_triangle_tree_tmp => p_triangle_tree%p_input_triangles
	p_triangle_tree%p_input_triangles => p_triangle_tree%p_output_triangles
	p_triangle_tree%p_output_triangles => p_triangle_tree_tmp

	! Get triangle coordinates and boundary data
    p_triangle_tree%i_fine_elements_index = 0
    p_triangle_tree%i_input_index = p_triangle_tree%i_input_min
    first_fem_triangle%r_coords(:,:) = p_triangle_tree%r_coords(:,:)
    first_fem_triangle%l_boundary_data(:,:) = p_triangle_tree%l_boundary_data(:,:)

	first_fem_triangle%triangle = p_triangle_tree%p_input_triangles(p_triangle_tree%i_input_index)
	p_triangle_tree%i_input_index = p_triangle_tree%i_input_index + 1

	first_fem_triangle%i_plotter_type = p_triangle_tree%i_matrix_type
	first_fem_triangle%i_turtle_type = p_triangle_tree%i_turtle_type
	first_fem_triangle%l_color = p_triangle_tree%l_color
	first_fem_triangle%i_depth = p_triangle_tree%i_depth

	if (is_leaf_triangle(first_fem_triangle)) then
		p_triangle_tree%i_fine_elements_index = p_triangle_tree%i_fine_elements_index + 1
	endif
end subroutine read_first_fem_triangle

!
! Reads the first child fem_triangle
!
subroutine read_first_child(p_triangle_tree, parent_fem_triangle, child_fem_triangle)
	type (triangle_tree), intent(inout), target			:: p_triangle_tree 		! intent(inout) (for the contents of the target)
	type (fem_triangle), intent(inout), target			:: parent_fem_triangle
	type (fem_triangle), intent(out)					:: child_fem_triangle

	integer(kind = 1), dimension(1 : 2, -3 : 3), parameter	:: sfc_child_type = reshape([ 1,-2, -1,-3, -2,-3,	0, 0, 2,-3,	1, 3, 1, 2], [2, 7])
	integer(kind = 1), dimension(-8:8), parameter 	        :: i_plotter_child_type = [ 3, 2, 1, 8, 7, 6, 5, 4,     0,  -6, -7, -8, -1, -2, -3, -4, -5 ]

	! only first child computes parent bisection
	parent_fem_triangle%r_bisection = 0.5_GRID_SR * (parent_fem_triangle%r_coords(:, 1) + parent_fem_triangle%r_coords(:, 3))

	child_fem_triangle%triangle = p_triangle_tree%p_input_triangles(p_triangle_tree%i_input_index)
	p_triangle_tree%i_input_index = p_triangle_tree%i_input_index + 1

	! Get triangle coordinates, boundary data and SFC_type from Parent
	! SFC Nodes (r_coords: 1, 4, 2)
	child_fem_triangle%r_coords(:,1) = parent_fem_triangle%r_coords(:,1)
	child_fem_triangle%r_coords(:,2) = parent_fem_triangle%r_bisection(:)
	child_fem_triangle%r_coords(:,3) = parent_fem_triangle%r_coords(:,2)

	! SFC Edges (l_boundary_data: internal, 3, 2)
	child_fem_triangle%l_boundary_data(:,1) = .false.	! internal edge
	child_fem_triangle%l_boundary_data(:,2) = parent_fem_triangle%l_boundary_data(:,3)
	child_fem_triangle%l_boundary_data(:,3) = parent_fem_triangle%l_boundary_data(:,2)

	! SFC type
	child_fem_triangle%i_turtle_type = sfc_child_type(1, parent_fem_triangle%i_turtle_type)
	! reference color
	child_fem_triangle%l_color = .NOT. parent_fem_triangle%l_color

	child_fem_triangle%i_plotter_type = i_plotter_child_type(parent_fem_triangle%i_plotter_type)

	! grid level
	child_fem_triangle%i_depth = parent_fem_triangle%i_depth + 1

	if (is_leaf_triangle(child_fem_triangle)) then
		p_triangle_tree%i_fine_elements_index = p_triangle_tree%i_fine_elements_index + 1
	endif
end subroutine read_first_child

!
! Reads the second child fem_triangle
!
subroutine read_second_child(p_triangle_tree, parent_fem_triangle, child_fem_triangle)
	type (triangle_tree), intent(inout), target	:: p_triangle_tree 		! intent(inout) (for the contents of the target)
	type (fem_triangle), intent(in), target		:: parent_fem_triangle
	type (fem_triangle), intent(out)			:: child_fem_triangle

	integer(kind = 1), dimension(1 : 2, -3 : 3), parameter	:: sfc_child_type = reshape([ 1,-2, -1,-3, -2,-3,	0, 0, 2,-3,	1, 3, 1, 2], [2, 7])
	integer(kind = 1), dimension(-8:8), parameter 	        :: i_plotter_child_type = [ 3, 2, 1, 8, 7, 6, 5, 4,     0,  -6, -7, -8, -1, -2, -3, -4, -5 ]

	child_fem_triangle%triangle = p_triangle_tree%p_input_triangles(p_triangle_tree%i_input_index)
	p_triangle_tree%i_input_index = p_triangle_tree%i_input_index + 1

	! Get triangle coordinates, boundary data and SFC_type
	! SFC Nodes (r_coords: 2, 4, 3)
	child_fem_triangle%r_coords(:,1) = parent_fem_triangle%r_coords(:,2)
	child_fem_triangle%r_coords(:,2) = parent_fem_triangle%r_bisection(:)
	child_fem_triangle%r_coords(:,3) = parent_fem_triangle%r_coords(:,3)

	! SFC Edges (l_boundary_data: 2, 1, internal)
	child_fem_triangle%l_boundary_data(:,1) = parent_fem_triangle%l_boundary_data(:,2)
	child_fem_triangle%l_boundary_data(:,2) = parent_fem_triangle%l_boundary_data(:,1)
	child_fem_triangle%l_boundary_data(:,3) = .false.	! internal edge

	! SFC type
	child_fem_triangle%i_turtle_type = sfc_child_type(2, parent_fem_triangle%i_turtle_type)
	! reference color
	child_fem_triangle%l_color = .NOT. parent_fem_triangle%l_color

	child_fem_triangle%i_plotter_type = -i_plotter_child_type(-parent_fem_triangle%i_plotter_type)

	! grid level
	child_fem_triangle%i_depth = parent_fem_triangle%i_depth + 1

	if (is_leaf_triangle(child_fem_triangle)) then
		p_triangle_tree%i_fine_elements_index = p_triangle_tree%i_fine_elements_index + 1
	endif
end subroutine read_second_child


function sfc_char2type(sfc_char) result(sfc_type)
	CHARACTER (LEN = 3), intent(in)				:: sfc_char
	integer (kind = GRID_SI)					:: sfc_type

		SELECT CASE (sfc_char)
		CASE ('H')
			sfc_type = H
		CASE ('V')
			sfc_type = V
		CASE ('K')
			sfc_type = K
		end SELECT
end function sfc_char2type


!
! Generates a triangle tree uniformly refined up to 'i_max_depth', with starting coordinates and boundary data
! This structure is ready to be read by 'get_next_triangle()'
!
subroutine generate_triangle_tree(p_triangle_tree, i_depth, i_max_depth, r_coords, l_boundary_data, i_turtle_type, l_color)
	implicit none
	type (triangle_tree), intent(inout), target	:: p_triangle_tree		! intent(out) (for the contents of the target)
	integer (kind = 1), intent(in)		        :: i_depth
	integer (kind = 1), intent(in)		        :: i_max_depth
	real (kind = GRID_SR), dimension(2,3), intent(in)	:: r_coords
	logical (kind = GRID_SL), dimension(2,3), intent(in):: l_boundary_data
	integer (kind = 1), intent(in)				:: i_turtle_type
	integer (kind = 1), intent(in)		        :: l_color
	! local variables
	integer (kind = GRID_SI)					:: i_element_index
	integer (kind = GRID_SI)					:: i_max_elements
	type (triangle)								:: tmp_triangle
	real (kind = GRID_SR), dimension(2)			:: r_diff
    integer             	                    :: i_error

	! Allocate memory for p_triangle_tree members (triangles, temp_fem_triangles)
	i_max_elements = 2**i_depth - 1

    allocate(p_triangle_tree%p_input_triangles(i_max_elements), stat = i_error); assert_eq(i_error, 0)
	allocate(p_triangle_tree%p_output_triangles(i_max_elements), stat = i_error); assert_eq(i_error, 0)

	p_triangle_tree%i_max_depth = i_max_depth
	p_triangle_tree%r_coords(:,:) = r_coords(:,:)
#if defined(_ASSERT)
	r_diff = r_coords(1:2, 1) - r_coords(1:2, 2)
	p_triangle_tree%r_leg_size = sqrt(dot_product( r_diff, r_diff))
!PRINT *, rank_MPI, "triangle with coords", r_coords, "has r_leg_size", p_triangle_tree%r_leg_size
!PRINT *, "sqrt(2)", sqrt_2, "sqrt(2)/2", sqrt_2_half
#endif
	p_triangle_tree%l_boundary_data(:,:) = l_boundary_data(:,:)
	p_triangle_tree%i_turtle_type = i_turtle_type

	p_triangle_tree%l_color = l_color
	p_triangle_tree%i_depth = i_depth

	i_element_index = 1
	call generate_triangle_stream(i_depth, 1, i_element_index, p_triangle_tree)

	! i_input_min, i_input_max, i_output_min, i_output_max
	p_triangle_tree%i_input_min = 1
	p_triangle_tree%i_input_max = i_element_index - 1
	p_triangle_tree%i_output_min = 1
	p_triangle_tree%i_output_max = i_element_index - 1


	p_triangle_tree%i_num_elements = i_element_index - 1
	p_triangle_tree%i_num_fine_elements = i_element_index / 2

	p_triangle_tree%i_my_fine_elements_start = 1
	p_triangle_tree%i_my_fine_elements_end = p_triangle_tree%i_num_fine_elements
	p_triangle_tree%i_fine_elements_index = 0

	p_triangle_tree%i_count_additional_triangles = 0
	p_triangle_tree%i_count_removed_triangles = 0
end subroutine generate_triangle_tree
!
! Generates a triangle tree uniformly refined up to 'i_max_depth', with starting coordinates and boundary data
! This structure is ready to be read by 'get_next_triangle()'
!
subroutine delete_triangle_tree(p_triangle_tree)
	implicit none
	type (triangle_tree), intent(inout), target	:: p_triangle_tree
    integer             	                    :: i_error

    deallocate(p_triangle_tree%p_input_triangles, stat = i_error); assert_eq(i_error, 0)
	deallocate(p_triangle_tree%p_output_triangles, stat = i_error); assert_eq(i_error, 0)
end subroutine delete_triangle_tree


! Generate a triangle stream by recursively refining the triangle with bisection
! Here the first (biggest) triangle will be on the 'p_elements(1)' position
RECURSIVE subroutine generate_triangle_stream(i_max_depth, i_depth, i_element_index, p_triangle_tree)
	implicit none
	integer (kind = 1), intent(in)		        :: i_max_depth
	integer (kind = 1), intent(in)		        :: i_depth
	integer (kind = GRID_SI), intent(inout)		:: i_element_index
	type (triangle_tree), intent(inout), target	:: p_triangle_tree		! intent(out) (for the contents of the target)

	if (i_depth < i_max_depth) then ! this triangle will have children

		p_triangle_tree%p_input_triangles(i_element_index) = triangle(TRIANGLE_NODE)
		p_triangle_tree%p_output_triangles(i_element_index) = triangle(TRIANGLE_NODE)

		i_element_index = i_element_index + 1
		CALL generate_triangle_stream(i_max_depth, i_depth + 1_1, i_element_index, p_triangle_tree)
		CALL generate_triangle_stream(i_max_depth, i_depth + 1_1, i_element_index, p_triangle_tree)
	else ! this triangle will be a Leaf

		p_triangle_tree%p_input_triangles(i_element_index) = triangle(TRIANGLE_LEAF)
		p_triangle_tree%p_output_triangles(i_element_index) = triangle(TRIANGLE_LEAF)

		i_element_index = i_element_index + 1
	endif
end subroutine generate_triangle_stream



subroutine refine_triangle(p_triangle_tree, refine_hypotenuse, refine_left_leg, refine_right_leg)
	type (triangle_tree), intent(inout), target	:: p_triangle_tree
	logical (kind = GRID_SL)					:: refine_hypotenuse	! .true. if it should be refined
	logical (kind = GRID_SL)					:: refine_left_leg		! .true. if it should be refined
	logical (kind = GRID_SL)					:: refine_right_leg		! .true. if it should be refined
	! local variables
	type (triangle), target						:: tri

	if (refine_hypotenuse) then
		if (refine_left_leg) then
			if (refine_right_leg) then	! case 3
!PRINT *, "refine_triangle() case 3"
				tri%i_triangle = TRIANGLE_LEAF
				call write_triangle_out(p_triangle_tree, tri)

				tri%i_triangle = TRIANGLE_LEAF
				call write_triangle_out(p_triangle_tree, tri)

				tri%i_triangle = TRIANGLE_NODE
				call write_triangle_out(p_triangle_tree, tri)

				tri%i_triangle = TRIANGLE_LEAF
				call write_triangle_out(p_triangle_tree, tri)

				tri%i_triangle = TRIANGLE_LEAF
				call write_triangle_out(p_triangle_tree, tri)

				tri%i_triangle = TRIANGLE_NODE
				call write_triangle_out(p_triangle_tree, tri)

				p_triangle_tree%i_num_elements = p_triangle_tree%i_num_elements + 6
				p_triangle_tree%i_num_fine_elements = p_triangle_tree%i_num_fine_elements + 3
			else	! case 2a
!PRINT *, "refine_triangle() case 2a"
				tri%i_triangle = TRIANGLE_LEAF
				call write_triangle_out(p_triangle_tree, tri)

				tri%i_triangle = TRIANGLE_LEAF
				call write_triangle_out(p_triangle_tree, tri)

				tri%i_triangle = TRIANGLE_NODE
				call write_triangle_out(p_triangle_tree, tri)

				tri%i_triangle = TRIANGLE_LEAF

				call write_triangle_out(p_triangle_tree, tri)


				p_triangle_tree%i_num_elements = p_triangle_tree%i_num_elements + 4
				p_triangle_tree%i_num_fine_elements = p_triangle_tree%i_num_fine_elements + 2
			endif
		else if (refine_right_leg) then	! case 2b
!PRINT *, "refine_triangle() case 2b"
				tri%i_triangle = TRIANGLE_LEAF
				call write_triangle_out(p_triangle_tree, tri)

				tri%i_triangle = TRIANGLE_LEAF
				call write_triangle_out(p_triangle_tree, tri)

				tri%i_triangle = TRIANGLE_LEAF
				call write_triangle_out(p_triangle_tree, tri)

				tri%i_triangle = TRIANGLE_NODE
				call write_triangle_out(p_triangle_tree, tri)


				p_triangle_tree%i_num_elements = p_triangle_tree%i_num_elements + 4
				p_triangle_tree%i_num_fine_elements = p_triangle_tree%i_num_fine_elements + 2
		else	! case 1
!PRINT *, "refine_triangle() case 1"
			tri%i_triangle = TRIANGLE_LEAF
			call write_triangle_out(p_triangle_tree, tri)

			tri%i_triangle = TRIANGLE_LEAF
			call write_triangle_out(p_triangle_tree, tri)


			p_triangle_tree%i_num_elements = p_triangle_tree%i_num_elements + 2
			p_triangle_tree%i_num_fine_elements = p_triangle_tree%i_num_fine_elements + 1
		endif
	endif
end subroutine refine_triangle


! Writes fem_triangle to output
subroutine write_triangle_out(p_triangle_tree, p_triangle)
	type (triangle_tree), intent(inout)			:: p_triangle_tree
	type (triangle), intent(in)					:: p_triangle

	p_triangle_tree%i_output_index = p_triangle_tree%i_output_index + 1
	p_triangle_tree%p_output_triangles(p_triangle_tree%i_output_index) = p_triangle
end subroutine write_triangle_out

!=============================================================
!
! Returns coordinates of the "i_index"-th node of a "fem_triangle"
! PREcondition: "i_index" must be in range 1:3, as this will not be checked in the function !!!
!
function get_node_coords(fem_tri, i_index) result (r_coordinates)
	type (fem_triangle), intent(in)				:: fem_tri
	integer, intent(in)							:: i_index
	real (kind = GRID_SR), dimension(2)			:: r_coordinates

	r_coordinates(:) = fem_tri%r_coords(:,i_index)
end function get_node_coords

function get_color_edge_color(sfc_type, base_color) result(color)
	integer(kind = 1), intent(in)			:: sfc_type
	logical (kind = GRID_SL), intent(in)	:: base_color
	logical (kind = GRID_SL)				:: color

	select case(sfc_type)
		case(K, H)
			color = base_color
		case (V)
			color = .not. base_color
	end select
end function get_color_edge_color

function get_current_edge_color(sfc_type, base_color) result(color)
	integer(kind = 1), intent(in)			:: sfc_type
	logical (kind = GRID_SL), intent(in)	:: base_color
	logical (kind = GRID_SL)				:: color

	select case(sfc_type)
		case(V, K)
			color = base_color
		case (H)
			color = .not. base_color
	end select
end function get_current_edge_color

function get_next_edge_color(sfc_type, base_color) result(color)
	integer(kind = 1), intent(in)			:: sfc_type
	logical (kind = GRID_SL), intent(in)	:: base_color
	logical (kind = GRID_SL)				:: color

	select case(sfc_type)
		case(V, H)
			color = base_color
		case (K)
			color = .not. base_color
	end select
end function get_next_edge_color

function get_color_edge_index(sfc_type) result(color_edge_index)
	integer(kind = 1), intent(in)			:: sfc_type
	integer (kind = 1)						:: color_edge_index

	select case(sfc_type)
		case(-3, 3)
			color_edge_index = RIGHT_EDGE
		case(-2, 2)
			color_edge_index = HYPOTENUSE
		case(-1, 1)
			color_edge_index = LEFT_EDGE
		case default
			color_edge_index = -1
	end select
end function get_color_edge_index

function get_previous_edge_index(sfc_type) result(previous_edge_index)
	integer(kind = 1), intent(in)			:: sfc_type
	integer (kind = 1)						:: previous_edge_index

	select case(sfc_type)
		case(-3, 3, -2, 2)
			previous_edge_index = LEFT_EDGE
		case(-1, 1)
			previous_edge_index = HYPOTENUSE
		case default
			previous_edge_index = -1
	end select
end function get_previous_edge_index

function get_next_edge_index(sfc_type) result(next_edge_index)
	integer(kind = 1), intent(in)			:: sfc_type
	integer (kind = 1)						:: next_edge_index

	select case(sfc_type)
		case(-3, 3)
			next_edge_index = HYPOTENUSE
		case(-2, 2, -1, 1)
			next_edge_index = RIGHT_EDGE
		case default
			next_edge_index = -1
	end select
end function get_next_edge_index

function is_internal_edge(fem_tri, i_edge_index) !result(is_dirichlet)
	type (fem_triangle), intent(in)				:: fem_tri
	integer (kind = 1), intent(in)				:: i_edge_index
	logical (kind = GRID_SL)					:: is_internal_edge

	is_internal_edge =	((fem_tri%l_boundary_data(1,i_edge_index) .eqv. internal_edge(1)) .and. &
						 (fem_tri%l_boundary_data(2,i_edge_index) .eqv. internal_edge(2)))
end function is_internal_edge

function is_dirichlet_boundary_edge(fem_tri, i_edge_index) !result(is_dirichlet)
	type (fem_triangle), intent(in)				:: fem_tri
	integer (kind = 1), intent(in)				:: i_edge_index
	logical (kind = GRID_SL)					:: is_dirichlet_boundary_edge

	is_dirichlet_boundary_edge =	((fem_tri%l_boundary_data(1,i_edge_index) .eqv. dirichlet_boundary(1)) .and. &
									 (fem_tri%l_boundary_data(2,i_edge_index) .eqv. dirichlet_boundary(2)))
end function is_dirichlet_boundary_edge

function is_neumann_boundary_edge(fem_tri, i_edge_index) !result(is_neumann)
	type (fem_triangle), intent(in)				:: fem_tri
	integer (kind = 1), intent(in)				:: i_edge_index
	logical (kind = GRID_SL)					:: is_neumann_boundary_edge

	is_neumann_boundary_edge =	((fem_tri%l_boundary_data(1,i_edge_index) .eqv. neumann_boundary(1)) .and. &
								 (fem_tri%l_boundary_data(2,i_edge_index) .eqv. neumann_boundary(2)))
end function is_neumann_boundary_edge

function is_domain_boundary_edge(fem_tri, i_edge_index)
	type (fem_triangle), intent(in)				:: fem_tri
	integer (kind = 1), intent(in)				:: i_edge_index
	logical (kind = GRID_SL)					:: is_domain_boundary_edge

	is_domain_boundary_edge = (is_dirichlet_boundary_edge(fem_tri, i_edge_index) .or. is_neumann_boundary_edge(fem_tri, i_edge_index))
end function is_domain_boundary_edge


end MODULE SFC_grid
