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
#   define _CHUNK_SIZE 1_8
#endif

private
public _CNT

!array implementation

type _CNT
	_T, pointer					:: elements(:) => null()		!< element array
	_T, pointer					:: elements_alloc(:) => null()	!< element array in allocation order
	integer*8                   :: i_current_element = 0		!< stream current element

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

	procedure, private, pass :: current_element
	procedure, private, pass :: next_element
	procedure, private, pass :: read_element
	procedure, private, pass :: read_elements
	procedure, private, pass :: write_element
	procedure, private, pass :: write_elements
	procedure, private, pass :: add_element

	procedure, pass :: to_string

	generic :: current => current_element
	generic :: next => next_element
	generic :: read => read_element, read_elements
	generic :: write => write_element, write_elements
	generic :: add => add_element
end type

contains

!construction/destruction

!> resizes a self-managed stream
subroutine resize_by_value(stream, i_elements)
	class(_CNT), intent(inout)						:: stream			!< stream object
	integer*8, intent(in)                           :: i_elements

    _T, pointer                                     :: elements_temp(:)
	integer                                         :: i_error

    assert(.not. associated(stream%elements) .or. associated(stream%elements, stream%elements_alloc))

	if (associated(stream%elements)) then
        allocate(elements_temp(i_elements), stat = i_error); assert_eq(i_error, 0)

        if (stream%is_forward()) then
            elements_temp(1 : min(i_elements, size(stream%elements))) = stream%elements_alloc
            deallocate(stream%elements_alloc, stat = i_error); assert_eq(i_error, 0)
            stream%elements => elements_temp
        else
            elements_temp(max(1, i_elements - size(stream%elements) + 1) : i_elements) = stream%elements_alloc
            deallocate(stream%elements_alloc, stat = i_error); assert_eq(i_error, 0)
            stream%elements => elements_temp(i_elements : 1 : -1)
        end if

        stream%elements_alloc => elements_temp
    else
        allocate(stream%elements_alloc(i_elements), stat = i_error); assert_eq(i_error, 0)
        stream%elements => stream%elements_alloc
    end if
end subroutine

!> auto-resizes a self-managed stream
subroutine resize_auto(stream)
	class(_CNT), intent(inout)						:: stream			!< stream object

    _T, pointer                                     :: elements_temp(:)
	integer                                         :: i_error

	if (associated(stream%elements)) then
        call stream%resize(size(stream%elements) + _CHUNK_SIZE)
    else
        allocate(stream%elements_alloc(_CHUNK_SIZE), stat = i_error); assert_eq(i_error, 0)
    end if

	stream%elements => stream%elements_alloc
end subroutine

subroutine clear(stream)
	class(_CNT), intent(inout)						:: stream			!< stream object

	integer                                         :: i_error

    if (associated(stream%elements_alloc)) then
        deallocate(stream%elements_alloc, stat = i_error); assert_eq(i_error, 0)
    end if

    nullify(stream%elements)
end subroutine

!> attaches a stream to an array
pure subroutine attach(stream, elements)
	class(_CNT), intent(inout)						:: stream			!< stream object
	_T, target, intent(inout)			            :: elements(:)		!< target array

	stream%elements_alloc => elements
	stream%elements => elements
	stream%i_current_element = 0
end subroutine

!> unattaches a stream from its target
pure subroutine unattach(stream)
	class(_CNT), intent(inout)					    :: stream					!< stream object

	nullify(stream%elements_alloc)
	nullify(stream%elements)
end subroutine

!> trims a stream to all elements up to the current element
elemental subroutine trim_stream(stream)
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
	assert_le(stream%i_current_element, size(stream%elements))

	p_data => stream%elements(stream%i_current_element)
end function

!> increments the current element index, then accesses the current element in the stream
function next_element(stream) result(p_data)
	class(_CNT), intent(inout)			        :: stream					!< stream object
	_T, pointer									:: p_data					!< data pointer

	stream%i_current_element = stream%i_current_element + 1

	!check for overflow
	assert_ge(stream%i_current_element, 1)
	assert_le(stream%i_current_element, size(stream%elements))

	p_data => stream%elements(stream%i_current_element)
end function

!> increments the current element index, then reads the current element from the stream
subroutine read_element(stream, data)
	class(_CNT), intent(inout)					:: stream					!< stream object
	_T, intent(out)								:: data						!< data

	stream%i_current_element = stream%i_current_element + 1

	!check for overflow
	assert_ge(stream%i_current_element, 1)
	assert_le(stream%i_current_element, size(stream%elements))

	data = stream%elements(stream%i_current_element)
end subroutine

!> increments the current element index, then reads a block of elements from the stream
subroutine read_elements(stream, data)
	class(_CNT), intent(inout)					:: stream					!< stream object
	_T, intent(out)								:: data(:)					!< data

	stream%i_current_element = stream%i_current_element + size(data)

	!check for overflow
	assert_ge(stream%i_current_element, size(data))
	assert_le(stream%i_current_element, size(stream%elements))

	data = stream%elements(stream%i_current_element - size(data) + 1 : stream%i_current_element)
end subroutine

!> increments the current element index, then writes the current element to the stream
subroutine write_element(stream, data)
	class(_CNT), intent(inout)					:: stream					!< stream object
	_T, intent(in)								:: data						!< data

	stream%i_current_element = stream%i_current_element + 1

	!check for overflow
	assert_ge(stream%i_current_element, 1)
	assert_le(stream%i_current_element, size(stream%elements))

	stream%elements(stream%i_current_element) = data
end subroutine

!> increments the current element index, then writes a chunk of elements to the stream
subroutine write_elements(stream, data)
	class(_CNT), intent(inout)					:: stream					!< stream object
	_T, intent(in)								:: data(:)					!< data

	stream%i_current_element = stream%i_current_element + size(data)

	!check for overflow
	assert_ge(stream%i_current_element, size(data))
	assert_le(stream%i_current_element, size(stream%elements))

	stream%elements(stream%i_current_element - size(data) + 1 : stream%i_current_element) = data
end subroutine

!> adds an element to a self-managed stream. If necessary, the size of the stream is increased to fit the element.
subroutine add_element(stream, data)
	class(_CNT), intent(inout)					:: stream					!< stream object
	_T, intent(in)								:: data						!< data

	stream%i_current_element = stream%i_current_element + 1

	!resize if necessary
	if (stream%i_current_element > size(stream%elements)) then
        call stream%resize_auto()
    end if

    !check for overflow
	assert_ge(stream%i_current_element, 1)
	assert_le(stream%i_current_element, size(stream%elements))

	stream%elements(stream%i_current_element) = data
end subroutine

!> Merges two a self-managed streams
function merge_streams(stream1, stream2) result(stream)
	class(_CNT), intent(in)					    :: stream1, stream2     !< stream objects
	type(_CNT)              					:: stream               !< stream objects

	_T, pointer					                :: elements_temp(:)
	integer                                     :: total_size, i_error

    !merge array parts
    total_size = size(stream1%elements) + size(stream2%elements)

    if (total_size > 0) then
        allocate(elements_temp(total_size), stat = i_error); assert_eq(i_error, 0)

        if (associated(stream1%elements)) then
            elements_temp(1 : size(stream1%elements)) = stream1%elements
        end if

        if (associated(stream2%elements)) then
            elements_temp(total_size - size(stream2%elements) + 1 : total_size) = stream2%elements
        end if

        stream%elements_alloc => elements_temp
        stream%elements => elements_temp
    end if
end function

elemental subroutine reverse_stream(stream)
	class(_CNT), intent(inout)					:: stream					!< stream object

	_T, pointer					                :: elements_temp(:)

	elements_temp => stream%elements(size(stream%elements) : 1 : -1)

	stream%elements => elements_temp
	stream%i_current_element = 0
end subroutine

elemental subroutine reset_stream(stream)
	class(_CNT), intent(inout)					:: stream					!< stream object

	stream%i_current_element = 0
end subroutine

elemental function to_string(stream) result(str)
	class(_CNT), intent(in)						:: stream					!< stream object
	character (len = 32)						:: str

	if (size(stream%elements) > 0) then
		write(str, '(A, I0, A, I0)') "elements: ", size(stream%elements), " current: ", stream%i_current_element
	else
		write(str, '(A, I0)') "elements: ", size(stream%elements)
	endif
end function

function is_forward(stream)
	class(_CNT), intent(in)     :: stream					!< stream object
    logical                     :: is_forward

    assert(associated(stream%elements))
    assert_ge(size(stream%elements), 1)

    is_forward = loc(stream%elements(1)) .le. loc(stream%elements(size(stream%elements)))
end function

function get_stream_c_pointer(stream) result(ptr)
	class(_CNT), intent(inout)					:: stream					!< stream object
	_T, pointer					                :: ptr

    if (.not. associated(stream%elements) .or. size(stream%elements) .eq. 0) then
        ptr => null()
    else if (stream%is_forward()) then
        ptr => stream%elements(1)
    else
        ptr => stream%elements(size(stream%elements))
    end if
end function

#undef _CNT_DATA_TYPE
#undef _CNT_TYPE_NAME
