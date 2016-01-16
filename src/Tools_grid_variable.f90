! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE



!> Generic grid variable
!> Warning: this a template module body and requires preprocessor commands before inclusion.
!> Usage: Inside your module definition, insert
!>
!> #define _GV_TYPE_NAME		<value>
!> #define _GV_NAME				<value>
!> #define _GV_TYPE				<value>
!> #define _GV_PERSISTENT		<value>
!> #define _GV_CELL_SIZE		<value>
!> #define _GV_EDGE_SIZE		<value>
!> #define _GV_NODE_SIZE		<value>
!>
!> #include <this_file>
!>
!> where
!>
!> @item _GV_TYPE_NAME		Grid variable type name
!> @item _GV_NAME			Access pattern for the grid variable, for example "dof(:)%x"
!> @item _GV_TYPE			Grid variable base data type, for example "real".
!> @item _GV_PERSISTENT		if 'true' the variable is stored persistently on the grid, otherwise only temporarily
!> @item _GV_CELL_SIZE		Size of the cell data, may be 0 if no data is located in cells
!> @item _GV_EDGE_SIZE		Size of the edge data, may be 0 if no data is located on edges
!> @item _GV_NODE_SIZE		Size of the node data, may be 0 if no data is located on nodes
!>
!> The resulting grid variable is defined as _GV_TYPE_NAME
!> @author Oliver Meister

#define _CONC2(X, Y)							X ## _ ## Y
#define _PREFIX(P, X)							_CONC2(P, X)
#define _GV_(X)									_PREFIX(_GV_TYPE_NAME, X)
#define _GV										_GV_TYPE_NAME

#define _GV_SIZE								(_GV_CELL_SIZE + 3 * _GV_EDGE_SIZE + 3 * _GV_NODE_SIZE)

#if (_GV_PERSISTENT)
#	define _GV_SUBSET							data_pers
#else
#	define _GV_SUBSET							data_temp
#endif

#if !defined(_GV_ADD_OP)
#   define _GV_ADD_OP(x, y)                     x + y
#endif

#define _GV_ACCESS_E(name, entity, subset)		entity%subset%name
#define _GV_ACCESS(entity)						_GV_ACCESS_E(_GV_NAME, entity, _GV_SUBSET)

private
public :: _GV

!> dummy type for the grid variable used for polymorphic function calls
type _GV
	contains

	procedure, pass :: read_from_element
	procedure, pass :: read_from_cell
	procedure, pass :: read_from_edge
	procedure, pass :: read_from_node
	procedure, pass :: write_to_element
	procedure, pass :: write_to_cell
	procedure, pass :: write_to_edge
	procedure, pass :: write_to_node
	procedure, pass :: add_to_element
	procedure, pass :: add_to_cell
	procedure, pass :: add_to_edge
	procedure, pass :: add_to_node

	generic :: read => read_from_element, read_from_cell, read_from_edge, read_from_node
	generic :: write => write_to_element, write_to_cell, write_to_edge, write_to_node
	generic :: add => add_to_element, add_to_cell, add_to_edge, add_to_node
end type

contains

pure subroutine write_node(xn, xl)
    _GV_TYPE, intent(out)	:: xn(*)
    _GV_TYPE, intent(in)	:: xl(*)

    xn(1:_GV_NODE_SIZE) = xl(1:_GV_NODE_SIZE)
end subroutine

pure subroutine write_edge(xe, xl)
    _GV_TYPE, intent(out)	:: xe(*)
    _GV_TYPE, intent(in)	:: xl(*)

    xe(1:_GV_EDGE_SIZE) = xl(1:_GV_EDGE_SIZE)
end subroutine

pure subroutine write_cell(xc, xl)
    _GV_TYPE, intent(out)	:: xc(*)
    _GV_TYPE, intent(in)	:: xl(*)

    xc(1:_GV_CELL_SIZE) = xl(1:_GV_CELL_SIZE)
end subroutine

pure subroutine add_node(xn, xl)
    _GV_TYPE, intent(inout)	:: xn(*)
    _GV_TYPE, intent(in)	:: xl(*)

    integer :: i

    !we use a forall loop here as some gfortran versions do not support non-intrinsic base objects for array operations

    forall (i = 1:_GV_NODE_SIZE)
        xn(i) = _GV_ADD_OP(xn(i),  xl(i))
    end forall
end subroutine

pure subroutine add_edge(xe, xl)
    _GV_TYPE, intent(inout)	:: xe(*)
    _GV_TYPE, intent(in)	:: xl(*)

    integer :: i

    !we use a forall loop here as some gfortran versions do not support non-intrinsic base objects for array operations

    forall (i = 1:_GV_EDGE_SIZE)
        xe(i) = _GV_ADD_OP(xe(i), xl(i))
    end forall
end subroutine

pure subroutine add_cell(xc, xl)
    _GV_TYPE, intent(inout)	:: xc(*)
    _GV_TYPE, intent(in)	:: xl(*)

    integer :: i

    !we use a forall loop here as some gfortran versions do not support non-intrinsic base objects for array operations

    forall (i = 1:_GV_CELL_SIZE)
        xc(i) = _GV_ADD_OP(xc(i), xl(i))
    end forall
end subroutine

!********************************
!LFS implementation
!********************************

!> Reads the grid variable from the grid into the array xl
pure subroutine read_from_element(gv, element, xl)
	class(_GV), intent(in)			    :: gv
	type(t_element_base), intent(in)    :: element
	 _GV_TYPE, intent(out)	            :: xl(*)

	integer (kind = GRID_SI)		    :: i

#   if (_GV_NODE_SIZE > 0)
        do i = 1, 3
            call write_node(xl(_GV_NODE_SIZE * i - _GV_NODE_SIZE + 1), _GV_ACCESS(element%nodes(i)%ptr))
        end do
#   endif

#   if (_GV_EDGE_SIZE > 0)
        do i = 1, 3
            call write_edge(xl(_GV_EDGE_SIZE * i + 3 * _GV_NODE_SIZE - _GV_EDGE_SIZE + 1), _GV_ACCESS(element%edges(i)%ptr))
        end do
#   endif

#   if (_GV_CELL_SIZE > 0)
        call write_cell(xl(3 * _GV_NODE_SIZE + 3 * _GV_EDGE_SIZE + 1), _GV_ACCESS(element%cell))
#   endif
end subroutine

!> Writes data from the array xl into the grid variable
pure subroutine write_to_element(gv, element, xl)
	class(_GV), intent(in)			        :: gv
	type(t_element_base), intent(inout)     :: element
	 _GV_TYPE, intent(in)	                :: xl(*)

	integer (kind = GRID_SI)			    :: i

#   if (_GV_NODE_SIZE > 0)
        do i = 1, 3
            call write_node(_GV_ACCESS(element%nodes(i)%ptr), xl(_GV_NODE_SIZE * i - _GV_NODE_SIZE + 1))
        end do
#   endif

#   if (_GV_EDGE_SIZE > 0)
        do i = 1, 3
            call write_edge(_GV_ACCESS(element%edges(i)%ptr), xl(_GV_EDGE_SIZE * i + 3 * _GV_NODE_SIZE - _GV_EDGE_SIZE + 1))
        end do
#   endif

#   if (_GV_CELL_SIZE > 0)
        call write_cell(_GV_ACCESS(element%cell), xl(3 * _GV_NODE_SIZE + 3 * _GV_EDGE_SIZE + 1))
#   endif
end subroutine

!> Adds data from the array xl to the grid variable
pure subroutine add_to_element(gv, element, xl)
	class(_GV), intent(in)			    :: gv
	type(t_element_base), intent(inout) :: element
    _GV_TYPE, intent(in)                :: xl(*)

	integer (kind = GRID_SI)			:: i

#   if (_GV_NODE_SIZE > 0)
        do i = 1, 3
            call add_node(_GV_ACCESS(element%nodes(i)%ptr), xl(_GV_NODE_SIZE * i - _GV_NODE_SIZE + 1))
        end do
#   endif

#   if (_GV_EDGE_SIZE > 0)
        do i = 1, 3
            call add_edge(_GV_ACCESS(element%edges(i)%ptr), xl(_GV_EDGE_SIZE * i + 3 * _GV_NODE_SIZE - _GV_EDGE_SIZE + 1))
        end do
#   endif

#   if (_GV_CELL_SIZE > 0)
        call add_cell(_GV_ACCESS(element%cell), xl(3 * _GV_NODE_SIZE + 3 * _GV_EDGE_SIZE + 1))
#   endif
end subroutine

!> Reads the grid variable from the grid into the array xl
pure subroutine read_from_cell(gv, cell, xl)
	class(_GV), intent(in)			        :: gv
	type(t_cell_data_ptr), intent(in)       :: cell
    _GV_TYPE, intent(inout)	                :: xl(*)

#   if (_GV_CELL_SIZE > 0)
        call write_cell(xl, _GV_ACCESS(cell))
#   endif
end subroutine

!> Writes data from the array xl into the grid variable
pure subroutine write_to_cell(gv, cell, xl)
	class(_GV), intent(in)		            :: gv
	type(t_cell_data_ptr), intent(inout)    :: cell
    _GV_TYPE, intent(in)	                :: xl(*)

#   if (_GV_CELL_SIZE > 0)
        call write_cell(_GV_ACCESS(cell), xl)
#   endif
end subroutine

!> Adds data from the array xl to the grid variable
pure subroutine add_to_cell(gv, cell, xl)
	class(_GV), intent(in)			        :: gv
	type(t_cell_data_ptr), intent(inout)    :: cell
    _GV_TYPE, intent(in)                    :: xl(*)

#   if (_GV_CELL_SIZE > 0)
        call add_cell(_GV_ACCESS(cell), xl)
#   endif
end subroutine

!> Reads the grid variable from the grid into the array xl
pure subroutine read_from_edge(gv, edge, xl)
	class(_GV), intent(in)			    :: gv
	type(t_edge_data), intent(in)       :: edge
    _GV_TYPE, intent(inout)	            :: xl(*)

#   if (_GV_EDGE_SIZE > 0)
        call write_edge(xl, _GV_ACCESS(edge))
#   endif
end subroutine

!> Writes data from the array xl into the grid variable
pure subroutine write_to_edge(gv, edge, xl)
	class(_GV), intent(in)			        :: gv
	type(t_edge_data), intent(inout)        :: edge
    _GV_TYPE, intent(in)	                :: xl(*)

#   if (_GV_EDGE_SIZE > 0)
        call write_edge(_GV_ACCESS(edge), xl)
#   endif
end subroutine

!> Adds data from the array xl to the grid variable
pure subroutine add_to_edge(gv, edge, xl)
	class(_GV), intent(in)			    :: gv
	type(t_edge_data), intent(inout)    :: edge
    _GV_TYPE, intent(in)                :: xl(*)

#   if (_GV_EDGE_SIZE > 0)
        call add_edge(_GV_ACCESS(edge), xl)
#   endif
end subroutine

!> Reads the grid variable from the grid into the array xl
pure subroutine read_from_node(gv, node, xl)
	class(_GV), intent(in)			    :: gv
	type(t_node_data), intent(in)       :: node
    _GV_TYPE, intent(inout)	            :: xl(*)

#   if (_GV_NODE_SIZE > 0)
        call write_node(xl, _GV_ACCESS(node))
#   endif
end subroutine

!> Writes data from the array xl into the grid variable
pure subroutine write_to_node(gv, node, xl)
	class(_GV), intent(in)			    :: gv
	type(t_node_data), intent(inout)    :: node
    _GV_TYPE, intent(in)	            :: xl(*)

#   if (_GV_NODE_SIZE > 0)
        call write_node(_GV_ACCESS(node), xl)
#   endif
end subroutine

!> Adds data from the array xl to the grid variable
pure subroutine add_to_node(gv, node, xl)
	class(_GV), intent(in)			    :: gv
	type(t_node_data), intent(inout)    :: node
    _GV_TYPE, intent(in)                :: xl(*)

#   if (_GV_NODE_SIZE > 0)
        call add_node(_GV_ACCESS(node), xl)
#   endif
end subroutine

#undef _GV
#undef _GV_SUBSET
#undef _GV_TYPE_NAME
#undef _GV_TYPE
#undef _GV_NAME
#undef _GV_PERSISTENT
#undef _GV_ADD_OP
