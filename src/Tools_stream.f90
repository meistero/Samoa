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

#if !defined(_CHUNK_SIZE)
#   define _CHUNK_SIZE 1_SI
#endif

private
public _CNT

!array implementation

type _CNT
	_T, pointer				    :: elements(:) => null()		!< element array
	_T, pointer				    :: elements_alloc(:) => null()	!< element array in allocation order
	integer(kind = GRID_DI)     :: i_current_element = 0		!< stream current element
    logical                     :: forward = .true.
	contains

	procedure, pass :: resize => resize_by_value
	procedure, private, pass :: resize_auto
	procedure, pass :: clear

	procedure, pass :: attach
	procedure, pass :: unattach
	procedure, pass :: trim => trim_stream
	procedure, nopass :: merge => merge_streams

	procedure, pass :: reset => reset_stream
	procedure, pass :: reverse => reverse_stream
	procedure, pass :: get_c_pointer => get_stream_c_pointer
	procedure, pass :: is_forward
	procedure, pass :: get_size

	procedure, private, pass :: current_element
	procedure, private, pass :: next_element
	procedure, private, pass :: read_element
	procedure, private, pass :: write_element
	procedure, private, pass :: add_element

	procedure, pass :: to_string

	generic :: current => current_element
	generic :: next => next_element
	generic :: read => read_element
	generic :: write => write_element
	generic :: add => add_element
end type

contains

!helper functions

function alloc_wrapper(i_elements) result(data)
    integer(kind = GRID_SI), intent(in)     :: i_elements
    _T, pointer                             :: data(:)

#   if defined(_KMP_ALLOC)
	    integer (kind = KMP_POINTER_KIND)               :: i_data
        type(c_ptr)                                     :: c_ptr_data
        _T                                              :: dummy

        i_data = KMP_MALLOC(i_elements * sizeof(dummy))
        c_ptr_data = transfer(i_data, c_ptr_data)
        call C_F_POINTER(c_ptr_data, data, [i_elements])
#   else
	    integer     :: i_error

        allocate(data(i_elements), stat = i_error); assert_eq(i_error, 0)
#   endif
end function

subroutine free_wrapper(data)
    _T, pointer, intent(inout)  :: data(:)

#   if defined(_KMP_ALLOC)
	    integer (kind = KMP_POINTER_KIND)               :: i_data
        type(c_ptr)                                     :: c_ptr_data
        _T                                              :: dummy

        c_ptr_data = c_loc(data)
        i_data = transfer(c_ptr_data, i_data)
        call KMP_FREE(i_data)
#   else
	    integer     :: i_error

        deallocate(data, stat = i_error); assert_eq(i_error, 0)
#   endif
end subroutine


!construction/destruction

!> resizes a self-managed stream
subroutine resize_by_value(stream, i_elements)
	class(_CNT), intent(inout)						:: stream			!< stream object

	integer(kind = GRID_SI), intent(in)             :: i_elements
    _T, pointer                                     :: elements_temp_alloc(:), elements_temp(:)
	integer(kind = GRID_SI)                         :: i
	integer                                         :: i_error

    assert(i_elements .ge. 0)

#   if defined _ASSERT
        if (stream%is_forward() .or. stream%get_size() < 1) then
            assert(.not. associated(stream%elements) .or. associated(stream%elements, stream%elements_alloc))
        else
            assert(associated(stream%elements, stream%elements_alloc(stream%get_size() : 1 : -1)))
        end if
#   endif

    elements_temp_alloc => alloc_wrapper(i_elements)

    if (stream%is_forward() .or. i_elements < 1) then
        elements_temp => elements_temp_alloc
    else
        elements_temp => elements_temp_alloc(i_elements : 1 : -1)
    end if

	if (associated(stream%elements)) then
	    do i = 1, min(i_elements, size(stream%elements))
            elements_temp(i) = stream%elements(i)
        end do

        call free_wrapper(stream%elements_alloc)
    end if

    stream%elements => elements_temp
    stream%elements_alloc => elements_temp_alloc
end subroutine

!> auto-resizes a self-managed stream
subroutine resize_auto(stream)
	class(_CNT), intent(inout)						:: stream			!< stream object

    _T, pointer                                     :: elements_temp(:)
	integer                                         :: i_error

	call stream%resize(stream%get_size() + _CHUNK_SIZE)
end subroutine

subroutine clear(stream)
	class(_CNT), intent(inout)						:: stream			!< stream object

	integer                                         :: i_error

    if (associated(stream%elements_alloc)) then
        call free_wrapper(stream%elements_alloc)
    end if

    nullify(stream%elements)
end subroutine

!> attaches a stream to an array
subroutine attach(stream, elements, is_forward)
	class(_CNT), intent(inout)						:: stream			!< stream object
	_T, target, intent(inout)			            :: elements(:)		!< target array
	logical                                         :: is_forward       !< true if the array is in allocation order, false otherwise

    if (is_forward .or. size(elements) < 1) then
        stream%elements_alloc => elements
    else
        stream%elements_alloc => elements(size(elements) : 1 : -1)
    end if

    stream%elements => elements
	stream%forward = is_forward
	stream%i_current_element = 0

#   if defined _ASSERT
        if (stream%get_size() > 1) then
            assert(stream%is_forward() .eqv. loc(stream%elements(1)) .lt. loc(stream%elements(size(elements))))
        end if
#   endif
end subroutine

!> unattaches a stream from its target
subroutine unattach(stream)
	class(_CNT), intent(inout)					    :: stream					!< stream object

	nullify(stream%elements_alloc)
	nullify(stream%elements)
end subroutine

!> trims a stream to all elements up to the current element
subroutine trim_stream(stream)
	class(_CNT), intent(inout)						:: stream			!< stream object

	stream%elements => stream%elements(1 : stream%i_current_element)
end subroutine

!access

!> accesses the current element in the stream
!> (only valid after a call to read, write or next)
function current_element(stream) result(p_data)
	class(_CNT), intent(inout)			        :: stream					!< stream object
	_T, pointer									:: p_data					!< data pointer

	!check for overflow
	assert_ge(stream%i_current_element, 1)
	assert_le(stream%i_current_element, stream%get_size())

	p_data => stream%elements(stream%i_current_element)
end function

!> increments the current element index, then accesses the current element in the stream
function next_element(stream) result(p_data)
	class(_CNT), intent(inout)			        :: stream					!< stream object
	_T, pointer									:: p_data					!< data pointer

	stream%i_current_element = stream%i_current_element + 1

	!check for overflow
	assert_ge(stream%i_current_element, 1)
	assert_le(stream%i_current_element, stream%get_size())

	p_data => stream%elements(stream%i_current_element)
end function

!> increments the current element index, then reads the current element from the stream
subroutine read_element(stream, data)
	class(_CNT), intent(inout)					:: stream					!< stream object
	_T, intent(out)								:: data						!< data

	stream%i_current_element = stream%i_current_element + 1

	!check for overflow
	assert_ge(stream%i_current_element, 1)
	assert_le(stream%i_current_element, stream%get_size())

	data = stream%elements(stream%i_current_element)
end subroutine

!> increments the current element index, then writes the current element to the stream
subroutine write_element(stream, data)
	class(_CNT), intent(inout)					:: stream					!< stream object
	_T, intent(in)								:: data						!< data

	stream%i_current_element = stream%i_current_element + 1

	!check for overflow
	assert_ge(stream%i_current_element, 1)
	assert_le(stream%i_current_element, stream%get_size())

	stream%elements(stream%i_current_element) = data
end subroutine

!> adds an element to a self-managed stream. If necessary, the size of the stream is increased to fit the element.
subroutine add_element(stream, data)
	class(_CNT), intent(inout)					:: stream					!< stream object
	_T, intent(in)								:: data						!< data

	stream%i_current_element = stream%i_current_element + 1

	!resize if necessary
	if (stream%i_current_element > stream%get_size()) then
        call stream%resize_auto()
    end if

    !check for overflow
	assert_ge(stream%i_current_element, 1)
	assert_le(stream%i_current_element, stream%get_size())

	stream%elements(stream%i_current_element) = data
end subroutine

!> Merges two self-managed streams
function merge_streams(stream1, stream2) result(stream)
	class(_CNT), intent(in)					    :: stream1, stream2     !< stream objects
	type(_CNT)              					:: stream               !< stream objects

	_T, pointer					                :: elements_temp_alloc(:), elements_temp(:)
	integer (kind = GRID_SI)                    :: stream1_size, stream2_size, i

    assert(.not. associated(stream%elements))
    assert(.not. associated(stream%elements_alloc))
    assert_eqv(stream1%is_forward(), stream2%is_forward())

    !merge array parts
    stream1_size = stream1%get_size()
    stream2_size = stream2%get_size()

    if (stream1_size + stream2_size > 0) then
        elements_temp_alloc => alloc_wrapper(stream1_size + stream2_size)

        if (stream1%is_forward()) then
            elements_temp => elements_temp_alloc
        else
            elements_temp => elements_temp_alloc(stream1_size + stream2_size : 1 : -1)
        end if

        if (associated(stream1%elements)) then
            do i = 1, stream1_size
                elements_temp(i) = stream1%elements(i)
            end do
        end if

        if (associated(stream2%elements)) then
            do i = 1, stream2_size
                elements_temp(stream1_size + i) = stream2%elements(i)
            end do
        end if

        stream%elements_alloc => elements_temp_alloc
        stream%elements => elements_temp
        stream%forward = stream1%is_forward()
    end if
end function

subroutine reverse_stream(stream)
	class(_CNT), intent(inout)					:: stream					!< stream object

	_T, pointer					                :: elements_temp(:)
	integer (kind = GRID_SI)                    :: stream_size

    stream_size = stream%get_size()

    if (stream_size > 0) then
        elements_temp => stream%elements(stream_size : 1 : -1)
        stream%elements => elements_temp
    end if

	stream%i_current_element = 0
	stream%forward = .not. stream%forward

#   if defined(_ASSERT)
        if (stream%get_size() > 1) then
            assert_pure((loc(stream%elements(1)) < loc(stream%elements(size(stream%elements)))) .eqv. stream%is_forward())
        end if
#   endif
end subroutine

subroutine reset_stream(stream)
	class(_CNT), intent(inout)					:: stream					!< stream object

	stream%i_current_element = 0
end subroutine

elemental function to_string(stream) result(str)
	class(_CNT), intent(in)						:: stream					!< stream object
	character (len = 32)						:: str

	if (stream%get_size() > 0) then
		write(str, '(A, I0, A, I0)') "elements: ", stream%get_size(), " current: ", stream%i_current_element
	else
		write(str, '(A, I0)') "elements: ", stream%get_size()
	endif
end function

function is_forward(stream)
	class(_CNT), intent(in)     :: stream					!< stream object
    logical                     :: is_forward

#   if defined _ASSERT
        if (stream%get_size() > 1) then
            assert(loc(stream%elements(1)) .lt. loc(stream%elements(size(stream%elements))) .eqv. stream%forward)
        end if
#   endif

    is_forward = stream%forward
end function

!> Returns the size of the list
pure function get_size(stream) result(i_elements)
	class(_CNT), intent(in)     :: stream						!< list object
    integer (kind = GRID_SI)    :: i_elements

    if (.not. associated(stream%elements)) then
        i_elements = 0
    else
        i_elements = size(stream%elements)
    end if
end function

function get_stream_c_pointer(stream) result(ptr)
	class(_CNT), intent(in)					    :: stream					!< stream object
	_T, pointer					                :: ptr
	_T, target  :: dummy

    if (stream%get_size() .eq. 0) then
        ptr => dummy
    else if (stream%is_forward()) then
        ptr => stream%elements(1)
    else
        ptr => stream%elements(size(stream%elements))
    end if
end function

#undef _CNT_DATA_TYPE
#undef _CNT_TYPE_NAME
