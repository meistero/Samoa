! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE



!> Generic local function space
!> Warning: this a template module body and requires preprocessor commands before inclusion.
!> Usage: Inside your module definition, insert
!>
!> #define _LFS_type_NAME		<value>
!> #define _LFS_CELL_SIZE		<value>
!> #define _LFS_EDGE_SIZE		<value>
!> #define _LFS_NODE_SIZE		<value>
!>
!> #include <this_file>
!>
!> where
!>
!> @item _LFS_type_NAME		Local function space type name
!> @item _LFS_CELL_SIZE		Size of the local function space in the cell, may be 0 if the data is not located in the cell
!> @item _LFS_EDGE_SIZE		Size of the local function space on an edge, may be 0 if the data is not located on the edges
!> @item _LFS_NODE_SIZE		Size of the local function space on a node, may be 0 if the data is not located on the nodes
!>
!> The resulting lfs is defined as <lfs_name>
!> @author Oliver Meister

#define _CONC2(X, Y)							X ## _ ## Y
#define _PREFIX(P, X)							_CONC2(P, X)
#define _LFS_(X)								_PREFIX(_LFS_type_NAME, X)

#define _LFS									_LFS_type_NAME

PRIVATE
PUBLIC :: _LFS

!> dummy type for the local function space used for polymorphic function calls
type _LFS
	integer(kind = BYTE), DIMENSION(0)		:: i_dummy

	contains

    procedure, pass :: transform_element_data
    procedure, pass :: transform_element_data_v

	procedure, pass :: transform_edge_data
	procedure, pass :: transform_edge_data_v

	procedure, pass :: transform_cell_data
	procedure, pass :: transform_cell_data_v

    generic :: transform => transform_element_data, transform_element_data_v, transform_edge_data, transform_edge_data_v, transform_cell_data, transform_cell_data_v
end type

contains

!> Transforms cell, edge and node data into element order and vice versa
pure subroutine transform_element_data(lfs, element, xl)
	class(_LFS), intent(in)													:: lfs
	type(t_element_base), intent(in)										:: element

	_LFS_type, dimension(_LFS_CELL_SIZE), intent(inout)			:: xl

	integer (kind = GRID_SI)												:: i

#	if (_LFS_EDGE_SIZE > 1)
		do i = 1, 3
			call transform_edge_data(lfs, element%transform_data%plotter_data%edges(i), xl(3 * _LFS_NODE_SIZE + i * _LFS_EDGE_SIZE + 1 - _LFS_EDGE_SIZE : 3 * _LFS_NODE_SIZE + i * _LFS_EDGE_SIZE))
		end do
#	endif

#	if (_LFS_CELL_SIZE > 1)
		call transform_cell_data(lfs, element%transform_data%plotter_data, xl(3 * (_LFS_NODE_SIZE + _LFS_EDGE_SIZE) + 1 : 3 * (_LFS_NODE_SIZE + _LFS_EDGE_SIZE) + _LFS_CELL_SIZE))
#	endif
end subroutine

pure subroutine transform_element_data_v(lfs, element, xl)
	class(_LFS), intent(in)													:: lfs
	type(t_element_base), intent(in)										:: element

	_LFS_type, dimension(:, :), intent(inout)								:: xl

#	if (_LFS_EDGE_SIZE > 1)
		integer (kind = GRID_SI)											:: i

		do i = 1, 3
			call transform_edge_data_v(lfs, element%transform_data%plotter_data%edges(i), xl(:, 3 * _LFS_NODE_SIZE + i * _LFS_EDGE_SIZE + 1 - _LFS_EDGE_SIZE : 3 * _LFS_NODE_SIZE + i * _LFS_EDGE_SIZE))
		end do
#	endif

#	if (_LFS_CELL_SIZE > 1)
		call transform_cell_data_v(lfs, element%transform_data%plotter_data, xl(:, 3 * (_LFS_NODE_SIZE + _LFS_EDGE_SIZE) + 1 : 3 * (_LFS_NODE_SIZE + _LFS_EDGE_SIZE) + _LFS_CELL_SIZE))
#	endif
end subroutine

!> Transforms data from edge order into element order and vice versa (order is swapped if necessary)
pure subroutine transform_edge_data(lfs, edge, xl)
	class(_LFS), intent(in)												:: lfs
	type(t_edge_transform_data), intent(in)								:: edge

	_LFS_type, dimension(_LFS_EDGE_SIZE), intent(inout)					:: xl

	xl = xl(: : edge%orientation)
end subroutine

!> Transforms data from edge order into element order and vice versa (order is swapped if necessary)
pure subroutine transform_edge_data_v(lfs, edge, xl)
	class(_LFS), intent(in)												:: lfs
	type(t_edge_transform_data), intent(in)								:: edge

	_LFS_type, dimension(:, :), intent(inout)							:: xl

	xl = xl(:, : : edge%orientation)
end subroutine

!> Transforms data from cell order into element order and vice versa (order is swapped if necessary)
pure subroutine transform_cell_data(lfs, cell, xl)
	class(_LFS), intent(in)												:: lfs
	type(t_cell_transform_data), intent(in)								:: cell

	_LFS_type, dimension(_LFS_CELL_SIZE), intent(inout)					:: xl

	xl = xl(: : cell%orientation)
end subroutine

pure subroutine transform_cell_data_v(lfs, cell, xl)
	class(_LFS), intent(in)												:: lfs
	type(t_cell_transform_data), intent(in)								:: cell

	_LFS_type, dimension(:, :), intent(inout)							:: xl

	xl = xl(:, : : cell%orientation)
end subroutine

#undef _LFS
