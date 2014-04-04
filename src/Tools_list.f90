! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


!> Generic container template
!> Warning: this a template module body and requires preprocessor commands before inclusion.
!> Usage: Inside your module definition, insert
!>
!> #define _CNT_DATA_TYPE			<value>
!> #define _CNT_TYPE_NAME			<value>
!>
!> #include <this_file>
!>
!> where
!>
!> @item _CNT_DATA_TYPE				Any derived type, base data type for the container
!> @item _CNT_TYPE_NAME				Container type name, will be used prefix for all operations
!>
!> The resulting container is defined as _CNT_TYPE_NAME, the chunk as _CNT_TYPE_NAME_chunk, methods as _CNT_TYPE_NAME_<method>
!> @author Oliver Meister

!multiple levels of indirection are necessary to properly resolve the names
#define _CONC2(X, Y)			X ## _ ## Y
#define _PREFIX(P, X)			_CONC2(P, X)
#define _CNT_(X)				_PREFIX(_CNT_TYPE_NAME, X)

#define _CNT					_CNT_TYPE_NAME
#define _T						_CNT_DATA_TYPE

PRIVATE add_element, add_after, resize, remove_at, clear, merge, reverse, to_string, is_forward, get_size

type _CNT
	_T, pointer                 :: elements(:) => null(), elements_alloc(:) => null()
	logical                     :: forward = .true.

	contains

    procedure                   :: resize
	procedure, private, pass    :: add_element
	procedure, private, pass    :: add_after
	procedure, private, pass    :: remove_at
	procedure, pass             :: clear
	procedure, nopass           :: merge
	procedure, pass             :: reverse
	procedure, pass             :: to_string
	procedure, pass             :: is_forward
	procedure, pass             :: get_size

	generic :: add => add_after, add_element
	generic :: remove => remove_at
end type

contains


!list

!> Resizes a list
subroutine resize(list, i_elements, i_src_start_arg, i_dest_start_arg, i_count_arg)
	class(_CNT), intent(inout)		                :: list						!< list object
    integer (kind = GRID_SI), intent(in)            :: i_elements
    integer (kind = GRID_SI), intent(in), optional  :: i_src_start_arg, i_dest_start_arg, i_count_arg

    integer (kind = GRID_SI)                :: i_src_start, i_dest_start, i_count
	integer                                 :: i_error
    _T, pointer                             :: elements_temp(:)

    allocate(elements_temp(i_elements), stat = i_error); assert_eq(i_error, 0)

    if (present(i_src_start_arg)) then
        i_src_start = i_src_start_arg
    else
        i_src_start = 1
    end if

    if (present(i_dest_start_arg)) then
        i_dest_start = i_dest_start_arg
    else
        i_dest_start = 1
    end if

    if (present(i_count_arg)) then
        i_count = i_count_arg
    else
        i_count = min(i_elements, size(list%elements_alloc))
    end if

    if (associated(list%elements)) then
        elements_temp(i_dest_start : i_dest_start + i_count - 1) = list%elements_alloc(i_src_start : i_src_start + i_count - 1)
    end if

    list%elements_alloc => elements_temp

    if (list%is_forward()) then
        list%elements => elements_temp
    else
        list%elements => elements_temp(i_elements : 1 : -1)
    end if
end subroutine

!> Adds an element at the end of a list
subroutine add_element(list, element)
	class(_CNT), intent(inout)		        :: list						!< list object
	_T, intent(in)							:: element

	call list%add_after(size(list%elements), element)
end subroutine

!> Adds an element after the current location
subroutine add_after(list, i_current_element, element)
	class(_CNT), intent(inout)		        :: list						!< list object
    integer (kind = GRID_SI), intent(in)    :: i_current_element
	_T, intent(in)							:: element

	integer                                 :: i_error
    _T, pointer                             :: elements_temp(:)

    if (.not. associated(list%elements_alloc)) then
        allocate(elements_temp(1), stat = i_error); assert_eq(i_error, 0)
        elements_temp(1) = element
    else
        allocate(elements_temp(size(list%elements_alloc) + 1), stat = i_error); assert_eq(i_error, 0)
        elements_temp(1 : i_current_element) = list%elements_alloc(1 : i_current_element)
        elements_temp(i_current_element + 1) = element
        elements_temp(i_current_element + 2 : ) = list%elements_alloc(i_current_element + 1 : )

        if (size(list%elements_alloc) > 0) then
            deallocate(list%elements_alloc, stat = i_error); assert_eq(i_error, 0)
        else
            nullify(list%elements_alloc)
        end if
    end if

    list%elements_alloc => elements_temp

    if (list%forward) then
        list%elements => elements_temp
    else
        list%elements => elements_temp(size(elements_temp) : 1 : -1)
    end if
end subroutine

subroutine remove_at(list, i_current_element)
	class(_CNT), intent(inout)		                :: list						!< list object
	integer, intent(in)                     	    :: i_current_element

	integer                                         :: i_error
    _T, pointer                                     :: elements_temp(:)

    assert(associated(list%elements_alloc))
    assert_ge(i_current_element, 1)
    assert_le(i_current_element, size(list%elements_alloc))

    allocate(elements_temp(size(list%elements_alloc) - 1), stat = i_error); assert_eq(i_error, 0)
    elements_temp(1 : i_current_element - 1) = list%elements_alloc(1 : i_current_element - 1)
    elements_temp(i_current_element : size(list%elements_alloc) - 1) = list%elements_alloc(i_current_element + 1 : size(list%elements_alloc))
    deallocate(list%elements_alloc, stat = i_error); assert_eq(i_error, 0)

    list%elements_alloc => elements_temp

    if (list%forward) then
        list%elements => elements_temp
    else
        list%elements => elements_temp(size(elements_temp) : 1 : -1)
    end if
end subroutine

!> merges the list with a second list
function merge(list1, list2) result(list)
	class(_CNT), intent(in)		            :: list1, list2				!< list object
	type(_CNT)		                        :: list						!< list object

	_T, pointer				                :: elements_temp(:)
	integer                                 :: total_size, i_error

    !merge array parts

    total_size = 0
    if (associated(list1%elements)) then
        total_size = total_size + size(list1%elements)
    end if

    if (associated(list1%elements)) then
        total_size = total_size + size(list2%elements)
    end if

    if (total_size > 0) then
        allocate(elements_temp(total_size), stat = i_error); assert_eq(i_error, 0)

        if (associated(list1%elements)) then
            elements_temp(1 : size(list1%elements)) = list1%elements
        end if

        if (associated(list2%elements)) then
            elements_temp(total_size - size(list2%elements) + 1 : total_size) = list2%elements
        end if

        list%elements_alloc => elements_temp
        list%elements => elements_temp
    end if
end function

subroutine reverse(list)
	class(_CNT), intent(inout)		        :: list						!< list object

	_T, dimension(:), pointer				:: elements_temp

    !reverse array part
	if (associated(list%elements)) then
        elements_temp => list%elements(size(list%elements) : 1 : -1)
        list%elements => elements_temp
    end if

    list%forward = .not. list%forward

    if (size(list%elements) > 1) then
        assert_pure((loc(list%elements(1)) < loc(list%elements(size(list%elements)))) .eqv. list%forward)
    end if
end subroutine

subroutine clear(list)
	class(_CNT), intent(inout)		        :: list						!< list object
	integer                                 :: i_error

	if (associated(list%elements)) then
        nullify(list%elements)

        if (size(list%elements_alloc) > 0) then
            deallocate(list%elements_alloc, stat = i_error); assert_eq(i_error, 0)
        else
            nullify(list%elements_alloc)
        end if
    end if
end subroutine

elemental function to_string(list) result(str)
	class(_CNT), intent(in)						:: list
	character (len = 32)						:: str

	write(str, '(A, I0)') "elements: ", size(list%elements)
end function

function is_forward(list)
	class(_CNT), intent(in)     :: list					!< list object
    logical                     :: is_forward

    is_forward = list%forward
end function

!> Returns the size of the list
pure function get_size(list) result(i_elements)
	class(_CNT), intent(in)     :: list						!< list object
    integer (kind = GRID_SI)    :: i_elements

    if (.not. associated(list%elements)) then
        i_elements = 0
    else
        i_elements = size(list%elements)
    end if
end function

#undef _CNT_DATA_TYPE
#undef _CNT_TYPE_NAME
