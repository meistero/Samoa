! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE



!> Generic grid variable
!> Warning: this a template module body and requires preprocessor commands before inclusion.
!> Usage: Inside your module definition, insert
!>
!> #define _GV_type_NAME		<value>
!> #define _GV_NAME				<value>
!> #define _GV_type				<value>
!> #define _GV_COUNT			<value>
!> #define _GV_PERSISTENT		<value>
!> #define _LFS_type_NAME		<value>
!> #define _LFS_CELL_SIZE		<value>
!> #define _LFS_EDGE_SIZE		<value>
!> #define _LFS_NODE_SIZE		<value>
!>
!> #include <this_file>
!>
!> where
!>
!> @item _GV_type_NAME		Grid variable type name
!> @item _GV_NAME			Access pattern for the grid variable, for example "dof(:)%x"
!> @item _GV_type			Grid variable base data type, for example "real".
!> @item _GV_COUNT			Number of components of the grid variable
!> @item _GV_PERSISTENT		if 'true' the variable is stored persistently on the grid, otherwise only temporary
!> @item _LFS_type_NAME		Local function space type name, will be used as prefix for all operations
!> @item _LFS_CELL_SIZE		Size of the local function space in the cell, may be 0 if the data is not located in the cell
!> @item _LFS_EDGE_SIZE		Size of the local function space on an edge, may be 0 if the data is not located on the edges
!> @item _LFS_NODE_SIZE		Size of the local function space on a node, may be 0 if the data is not located on the nodes
!>
!> The resulting grid variable is defined as _GV_type_NAME
!> @author Oliver Meister

#define _CONC2(X, Y)							X ## _ ## Y
#define _PREFIX(P, X)							_CONC2(P, X)
#define _LFS_(X)								_PREFIX(_LFS_type_NAME, X)
#define _GV_(X)									_PREFIX(_GV_type_NAME, X)

#define _LFS									_LFS_type_NAME
#define _LFS_SIZE								(_LFS_CELL_SIZE + 3 * _LFS_EDGE_SIZE + 3 * _LFS_NODE_SIZE)

#define _GV										_GV_type_NAME

#if (_GV_PERSISTENT)
#	define _GV_SUBSET							data_pers
#else
#	define _GV_SUBSET							data_temp
#endif

#define _GV_ACCESS_E(name, entity, subset)		entity%subset%name
#define _GV_ACCESS(entity)						_GV_ACCESS_E(_GV_NAME, entity, _GV_SUBSET)

#define _GV_GET(entity, value)					value = _GV_ACCESS(entity)
#define _GV_SET(entity, value)					_GV_ACCESS(entity) = value
#define _GV_ADD(entity, value)					_GV_ACCESS(entity) = _GV_ACCESS(entity) + value

PRIVATE
PUBLIC :: _GV

!> dummy type for the grid variable used for polymorphic function calls
type _GV
	integer(kind = 1), DIMENSION(0)		:: i_dummy

	contains

	procedure, pass :: read => read_element
	procedure, pass :: write => write_element
	procedure, pass :: add => add_element
end type

contains

!********************************
!LFS implementation
!********************************

!> Reads the grid variable from the grid into the array xl
pure subroutine read_element(gv, element, xl)
	class(_GV), intent(in)										:: gv
	type(t_element_base), intent(in)							:: element

	integer (kind = GRID_SI)									:: i, j

#	if defined(_GV_COUNT)
		_GV_type, dimension(:, :), intent(out)	:: xl

		forall(i = 1 : 3, j = 1 : _GV_COUNT)
#			if (_LFS_NODE_SIZE == 1)
				xl(j, i) = _GV_ACCESS(element%nodes(i)%ptr)(j, 1)
#			elif (_LFS_NODE_SIZE > 1)
				xl(j, i * _LFS_NODE_SIZE + 1 - _LFS_NODE_SIZE : i * _LFS_NODE_SIZE) = &
					_GV_ACCESS(element%nodes(i)%ptr)(j, :)
#			endif

#			if (_LFS_EDGE_SIZE == 1)
				xl(j, 3 * _LFS_NODE_SIZE + i) = _GV_ACCESS(element%edges(i)%ptr)(j, 1)
#			elif (_LFS_EDGE_SIZE > 1)
				xl(j, 3 * _LFS_NODE_SIZE + i * _LFS_EDGE_SIZE + 1 - _LFS_EDGE_SIZE : 3 * _LFS_NODE_SIZE + i * _LFS_EDGE_SIZE) = &
					transform_edge_data(lfs, element%edges(i)%transform_data, &
					_GV_ACCESS(element%edges(i)%ptr)(j, :))
#			endif
		end forall

		forall(j = 1 : _GV_COUNT)
#			if (_LFS_CELL_SIZE == 1)
				xl(j, 3 * (_LFS_NODE_SIZE + _LFS_EDGE_SIZE) + 1) = _GV_ACCESS(element%cell)(j, 1)
#			elif (_LFS_CELL_SIZE > 1)
				xl(j, 3 * (_LFS_NODE_SIZE + _LFS_EDGE_SIZE) + 1 : 3 * (_LFS_NODE_SIZE + _LFS_EDGE_SIZE) + _LFS_CELL_SIZE) = &
					transform_cell_data(lfs, element%transform_data%plotter_data, &
					_GV_ACCESS(element%cell)(j, :))
#			endif
		end forall
#	else
		_GV_type, dimension(:), intent(out)				:: xl

		forall(i = 1 : 3)
#			if (_LFS_NODE_SIZE == 1)
				xl(i) = _GV_ACCESS(element%nodes(i)%ptr)(1)
#			elif (_LFS_NODE_SIZE > 1)
				xl(i * _LFS_NODE_SIZE + 1 - _LFS_NODE_SIZE : i * _LFS_NODE_SIZE) = &
					_GV_ACCESS(element%nodes(i)%ptr)
#			endif

#			if (_LFS_EDGE_SIZE == 1)
				xl(3 * _LFS_NODE_SIZE + i) = _GV_ACCESS(element%edges(i)%ptr)(1)
#			elif (_LFS_EDGE_SIZE > 1)
				xl(3 * _LFS_NODE_SIZE + i * _LFS_EDGE_SIZE + 1 - _LFS_EDGE_SIZE : 3 * _LFS_NODE_SIZE + i * _LFS_EDGE_SIZE) = &
					transform_edge_data(lfs, element%edges(i)%transform_data, &
					_GV_ACCESS(element%edges(i)%ptr))
#			endif
		end forall

#		if (_LFS_CELL_SIZE == 1)
			xl(3 * (_LFS_NODE_SIZE + _LFS_EDGE_SIZE) + 1) = _GV_ACCESS(element%cell)(1)
#		elif (_LFS_CELL_SIZE > 1)
			xl(3 * (_LFS_NODE_SIZE + _LFS_EDGE_SIZE) + 1 : 3 * (_LFS_NODE_SIZE + _LFS_EDGE_SIZE) + _LFS_CELL_SIZE) = &
				transform_cell_data(lfs, element%transform_data%plotter_data, &
				_GV_ACCESS(element%cell))
#		endif
#	endif
end subroutine

!> Writes data from the array xl into the grid variable
pure subroutine write_element(gv, element, xl)
	class(_GV), intent(in)										:: gv
	type(t_element_base), intent(inout)							:: element

	integer (kind = GRID_SI)									:: i, j

#	if defined(_GV_COUNT)
		_GV_type, dimension(:, :), intent(in)	:: xl

		forall (i = 1 : 3, j = 1 : _GV_COUNT)
#			if (_LFS_NODE_SIZE == 1)
					_GV_ACCESS(element%nodes(i)%ptr)(j, 1) = xl(j, i)
#			elif (_LFS_NODE_SIZE > 1)
					_GV_ACCESS(element%nodes(i)%ptr)(j, :) = &
						xl(j, i * _LFS_NODE_SIZE + 1 - _LFS_NODE_SIZE : i * _LFS_NODE_SIZE)
#			endif

#			if (_LFS_EDGE_SIZE == 1)
					_GV_ACCESS(element%edges(i)%ptr)(j, 1) = xl(j, 3 * _LFS_NODE_SIZE + i)
#			elif (_LFS_EDGE_SIZE > 1)
					_GV_ACCESS(element%edges(i)%ptr)(j, :) = &
						transform_edge_data(lfs, element%edges(i)%transform_data, &
						xl(3 * _LFS_NODE_SIZE + i * _LFS_EDGE_SIZE + 1 - _LFS_EDGE_SIZE : 3 * _LFS_NODE_SIZE + i * _LFS_EDGE_SIZE))
#			endif
		end forall

		forall (j = 1 : _GV_COUNT)
#			if (_LFS_CELL_SIZE == 1)
				_GV_ACCESS(element%cell)(j, 1) = xl(j, 3 * (_LFS_NODE_SIZE + _LFS_EDGE_SIZE) + 1)
#			elif (_LFS_CELL_SIZE > 1)
				_GV_ACCESS(element%cell)(j, :) = &
					transform_cell_data(lfs, element%transform_data%plotter_data, &
					xl(j, 3 * (_LFS_NODE_SIZE + _LFS_EDGE_SIZE) + 1 : 3 * (_LFS_NODE_SIZE + _LFS_EDGE_SIZE) + _LFS_CELL_SIZE))
#			endif
		end forall
#	else
		_GV_type, dimension(:), intent(in)				:: xl

		forall(i = 1 : 3)
#			if (_LFS_NODE_SIZE == 1)
					_GV_ACCESS(element%nodes(i)%ptr)(1) = xl(i)
#			elif (_LFS_NODE_SIZE > 1)
					_GV_ACCESS(element%nodes(i)%ptr) = &
						xl(i * _LFS_NODE_SIZE + 1 - _LFS_NODE_SIZE : i * _LFS_NODE_SIZE)
#			endif

#			if (_LFS_EDGE_SIZE == 1)
					_GV_ACCESS(element%edges(i)%ptr)(1) = xl(3 * _LFS_NODE_SIZE + i)
#			elif (_LFS_EDGE_SIZE > 1)
					_GV_ACCESS(element%edges(i)%ptr) = &
						transform_edge_data(lfs, element%edges(i)%transform_data, &
						xl(3 * _LFS_NODE_SIZE + i * _LFS_EDGE_SIZE + 1 - _LFS_EDGE_SIZE : 3 * _LFS_NODE_SIZE + i * _LFS_EDGE_SIZE))
#			endif
		end forall

#		if (_LFS_CELL_SIZE == 1)
			_GV_ACCESS(element%cell)(1) = xl(3 * (_LFS_NODE_SIZE + _LFS_EDGE_SIZE) + 1)
#		elif (_LFS_CELL_SIZE > 1)
			_GV_ACCESS(element%cell) = &
				transform_cell_data(lfs, element%transform_data%plotter_data, &
				xl(3 * (_LFS_NODE_SIZE + _LFS_EDGE_SIZE) + 1 : 3 * (_LFS_NODE_SIZE + _LFS_EDGE_SIZE) + _LFS_CELL_SIZE))
#		endif
#	endif
end subroutine

!> Adds data from the array xl to the grid variable
pure subroutine add_element(gv, element, xl)
	class(_GV), intent(in)										:: gv
	type(t_element_base), intent(inout)							:: element

	integer (kind = GRID_SI)									:: i, j

#	if defined(_GV_COUNT)
		_GV_type, dimension(:, :), intent(in)					:: xl

		forall (i = 1 : 3, j = 1 : _GV_COUNT)
#			if (_LFS_NODE_SIZE == 1)
				_GV_ACCESS(element%nodes(i)%ptr)(j, 1) = _GV_ACCESS(element%nodes(i)%ptr)(j, 1) + xl(j, i)
#			elif (_LFS_NODE_SIZE > 1)
				_GV_ACCESS(element%nodes(i)%ptr)(j, :) = _GV_ACCESS(element%nodes(i)%ptr)(j, :) + &
					xl(j, i * _LFS_NODE_SIZE + 1 - _LFS_NODE_SIZE : i * _LFS_NODE_SIZE)
#			endif

#			if (_LFS_EDGE_SIZE == 1)
				_GV_ACCESS(element%edges(i)%ptr)(j, 1) = _GV_ACCESS(element%edges(i)%ptr)(j, 1) + xl(j, 3 * _LFS_NODE_SIZE + i)
#			elif (_LFS_EDGE_SIZE > 1)
				_GV_ACCESS(element%edges(i)%ptr)(j, :) = _GV_ACCESS(element%edges(i)%ptr)(j, :) + &
					transform_edge_data(lfs, element%edges(i)%transform_data, &
					xl(3 * _LFS_NODE_SIZE + i * _LFS_EDGE_SIZE + 1 - _LFS_EDGE_SIZE : 3 * _LFS_NODE_SIZE + i * _LFS_EDGE_SIZE))
#			endif
		end forall

		forall (j = 1 : _GV_COUNT)
#			if (_LFS_CELL_SIZE == 1)
				_GV_ACCESS(element%cell)(j, 1) = _GV_ACCESS(element%cell)(j, 1) + xl(j, 3 * (_LFS_NODE_SIZE + _LFS_EDGE_SIZE) + 1)
#			elif (_LFS_CELL_SIZE > 1)
				_GV_ACCESS(element%cell)(j, :) = _GV_ACCESS(element%cell)(j, :) + &
					transform_cell_data(lfs, element%transform_data%plotter_data, &
					xl(j, 3 * (_LFS_NODE_SIZE + _LFS_EDGE_SIZE) + 1 : 3 * (_LFS_NODE_SIZE + _LFS_EDGE_SIZE) + _LFS_CELL_SIZE))
#			endif
		end forall
#	else
		_GV_type, dimension(:), intent(in)						:: xl

		forall(i = 1 : 3)
#			if (_LFS_NODE_SIZE == 1)
				_GV_ACCESS(element%nodes(i)%ptr)(1) = _GV_ACCESS(element%nodes(i)%ptr)(1) + xl(i)
#			elif (_LFS_NODE_SIZE > 1)
				_GV_ACCESS(element%nodes(i)%ptr) = _GV_ACCESS(element%nodes(i)%ptr) + &
					xl(i * _LFS_NODE_SIZE + 1 - _LFS_NODE_SIZE : i * _LFS_NODE_SIZE)
#			endif

#			if (_LFS_EDGE_SIZE == 1)
				_GV_ACCESS(element%edges(i)%ptr)(1) = _GV_ACCESS(element%edges(i)%ptr)(1) + xl(3 * _LFS_NODE_SIZE + i)
#			elif (_LFS_EDGE_SIZE > 1)
				_GV_ACCESS(element%edges(i)%ptr) = _GV_ACCESS(element%edges(i)%ptr) + &
					transform_edge_data(lfs, element%edges(i)%transform_data, &
					xl(3 * _LFS_NODE_SIZE + i * _LFS_EDGE_SIZE + 1 - _LFS_EDGE_SIZE : 3 * _LFS_NODE_SIZE + i * _LFS_EDGE_SIZE))
#			endif
		end forall

#		if (_LFS_CELL_SIZE == 1)
			_GV_ACCESS(element%cell)(1) = _GV_ACCESS(element%cell)(1) + xl(3 * (_LFS_NODE_SIZE + _LFS_EDGE_SIZE) + 1)
#		elif (_LFS_CELL_SIZE > 1)
			_GV_ACCESS(element%cell) = _GV_ACCESS(element%cell) + &
				transform_cell_data(lfs, element%transform_data%plotter_data, &
				xl(3 * (_LFS_NODE_SIZE + _LFS_EDGE_SIZE) + 1 : 3 * (_LFS_NODE_SIZE + _LFS_EDGE_SIZE) + _LFS_CELL_SIZE))
#		endif
#	endif
end subroutine

#if (_LFS_EDGE_SIZE > 1)
	!> Transforms data from edge order into element order and vice versa (order is swapped if necessary)
	pure function transform_edge_data(gv, edge, xe) result(xl)
		type(_GV), intent(in)								:: gv
		type(t_edge_transform_data), intent(in)				:: edge

		_GV_type, dimension(_LFS_EDGE_SIZE), intent(in)		:: xe
		_GV_type, dimension(_LFS_EDGE_SIZE)					:: xl

		xl = xe(: : edge%orientation)
	end function
#endif

#if (_LFS_CELL_SIZE > 1)
	!> Transforms data from cell order into element order and vice versa (order is swapped if necessary)
	pure function transform_cell_data(gv, cell, xc) result(xl)
		type(_GV), intent(in)								:: gv
		type(t_cell_transform_data), intent(in)				:: cell

		_GV_type, dimension(_LFS_CELL_SIZE), intent(in)		:: xc
		_GV_type, dimension(_LFS_CELL_SIZE)					:: xl

		xl = xc(: : cell%orientation)
	end function
#endif


#undef _LFS
#undef _LFS_SIZE
#undef _GV
#undef _GV_SUBSET
#undef _GV_type_NAME
#undef _GV_type
#undef _GV_NAME
#undef _GV_COUNT
#undef _GV_PERSISTENT

