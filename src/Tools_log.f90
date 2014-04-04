! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

module Tools_mpi
#	if defined(_MPI)
        include 'mpif.h'
#	else
		enum, bind(c)
		    enumerator :: MPI_REQUEST_NULL = 0, MPI_STATUS_SIZE
		    enumerator :: MPI_MAX, MPI_MIN, MPI_SUM, MPI_PROD, MPI_LAND, MPI_LOR
		    enumerator :: MPI_ADDRESS_KIND = 1
		end enum
#   endif

    public
end module

MODULE Tools_log
    use Tools_mpi

	private
	public log_open_file, log_close_file, g_log_file_unit, prefix_sum, postfix_sum, reduce, scatter, raise_error, term_color, term_reset

	interface prefix_sum
        module procedure prefix_sum_i8_sequential
        module procedure prefix_sum_i16_sequential
        module procedure prefix_sum_r32_sequential
        module procedure prefix_sum_r64_sequential
        module procedure prefix_sum_r128_sequential
    end interface

 	interface postfix_sum
        module procedure postfix_sum_i8
        module procedure postfix_sum_i16
        module procedure postfix_sum_r32
        module procedure postfix_sum_r64
        module procedure postfix_sum_r128
    end interface

	interface reduce
        module procedure reduce_i8_mpi
        module procedure reduce_i8
        module procedure reduce_i16_mpi
        module procedure reduce_i16
        module procedure reduce_i32_mpi
        module procedure reduce_i32
        module procedure reduce_i64_mpi
        module procedure reduce_i64
        module procedure reduce_r32_mpi
        module procedure reduce_r32
        module procedure reduce_r64_mpi
        module procedure reduce_r64
        module procedure reduce_r128_mpi
        module procedure reduce_r128
        module procedure reduce_l_mpi
        module procedure reduce_l
    end interface

	interface scatter
        module procedure scatter_i32
        module procedure scatter_i64
        module procedure scatter_r32
        module procedure scatter_r64
        module procedure scatter_r128
        module procedure scatter_l
        module procedure scatter_s
    end interface

	!> global file unit (there's only one log instance possible due to the lack of variadic functions and macros in Fortran,
	!> which makes it impossible to wrap WRITE and PRINT)
	integer (kind=4)    :: g_log_file_unit = 6
	integer, parameter  :: QP = kind(1.0q0)

	contains

	!> Opens an external log file where the log write commands are written to,
	!> if this is not called, the log is written to stdout instead
	subroutine log_open_file(s_file_name)
		character(*)				:: s_file_name

		integer 					:: i_error
		logical						:: l_is_open

		assert_eq(g_log_file_unit, 6)

		do g_log_file_unit = 10, huge(1)
			inquire(unit=g_log_file_unit, opened=l_is_open, iostat=i_error)

			if (i_error == 0) then
				if (.NOT. l_is_open) then
					open(unit=g_log_file_unit, file=s_file_name, status='replace', access='sequential', iostat=i_error); assert_eq(i_error, 0)
					flush(g_log_file_unit)
                    return
				end if
			end if
		end do
	end subroutine

	!> Closes an external log file,
	!> consecutive write calls are written to stdout instead
	subroutine log_close_file()
		integer 					:: i_error

		assert_ne(g_log_file_unit, 6)

		close(unit=g_log_file_unit, iostat=i_error); assert_eq(i_error, 0)
		g_log_file_unit = 6
	end subroutine

	function linear_search(v, s) result(i)
        integer, intent(in)                :: v(:)
        integer, intent(in)                :: s

        integer                            :: i

        i = minloc(abs(v - s), 1)

        assert_eq(v(i), s)
    end function

	function binary_search(v, s) result(i)
        integer, intent(in)                :: v(:)
        integer, intent(in)                :: s

        integer                            :: i_start, i_end, i
        integer, parameter                 :: vector_threshold = 4 - 1

        i_start = 1
        i_end = size(v)

        if (v(i_end) - v(i_start) .ge. 0) then
            do while (i_end - i_start > vector_threshold)
                i = (i_start + i_end) / 2

                select case (v(i) - s)
                    case(0)
                        return
                    case(1:)
                        i_end = i - 1
                    case(:-1)
                        i_start = i + 1
                end select
            end do
        else
            do while (i_end - i_start > vector_threshold)
                i = (i_start + i_end) / 2

                select case (s - v(i))
                    case(0)
                        return
                    case(1:)
                        i_end = i - 1
                    case(:-1)
                        i_start = i + 1
                end select
            end do
        end if

        !switch to a linear search for small arrays
        i = minloc(abs(v(i_start : i_end) - s), 1)

        assert_eq(v(i), s)
    end function

    function interpolation_search(v, s) result(i)
        integer, intent(in)                :: v(:)
        integer, intent(in)                :: s

        integer                            :: i_start, i_end, i
        integer, parameter                 :: vector_threshold = 4 - 1

        i_start = 1
        i_end = size(v)

        if (v(i_end) - v(i_start) .ge. 0) then
            do while (i_end - i_start > vector_threshold)
                i = i_start + ((i_end - i_start) * (s - v(i_start))) / (v(i_end) - v(i_start))

                select case (v(i) - s)
                    case(0)
                        return
                    case (:-1)
                        i_start = i + 1
                    case (1:)
                        i_end = i - 1
                end select
            end do
        else
            do while (i_end - i_start > vector_threshold)
                i = i_start + ((i_end - i_start) * (s - v(i_start))) / (v(i_end) - v(i_start))

                select case (s - v(i))
                    case(0)
                        return
                    case (:-1)
                        i_start = i + 1
                    case (1:)
                        i_end = i - 1
                end select
            end do
        end if

        !switch to a linear search for small arrays
        i = minloc(abs(v(i_start : i_end) - s), 1)

        assert_eq(v(i), s)
    end function

	pure subroutine postfix_sum_i8(y, x)
        integer(kind = selected_int_kind(8)), intent(inout)  	:: y(:)
        integer(kind = selected_int_kind(8)), intent(in)  	    :: x(:)

		call prefix_sum(y, x(size(x) : 1 : -1))
	end subroutine

	pure subroutine prefix_sum_i8_sequential(y, x)
        integer(kind = selected_int_kind(8)), intent(inout)  	:: y(:)
        integer(kind = selected_int_kind(8)), intent(in)  	    :: x(:)
        integer                  	                            :: i

        if (size(x) == 0) then
            return
        end if

        assert_pure(size(y) .eq. size(x))

        y(1) = x(1)

        do i = 2, size(x)
            y(i) = y(i - 1) + x(i)
        end do
	end subroutine

	pure subroutine prefix_sum_i8_vectorized(y, x)
        integer(kind = selected_int_kind(8)), intent(inout)  	:: y(:)
        integer(kind = selected_int_kind(8)), intent(in)  	    :: x(:)
        integer               	                                :: i, k, n

        if (size(x) == 0) then
            return
        end if

        assert_pure(size(y) .eq. size(x))

        n = size(x)

        !Copy array
        y = x

        !Up-Sweep phase
        i = 1
        do while (i .le. n/2)
            do k = 2 * i, n, 2 * i
                y(k) = y(k) + y(k - i)
            end do

            i = 2 * i
        end do

        !Down-Sweep phase
        do while (i .ge. 2)
            i = i / 2

            do k = 3 * i, n, 2 * i
                y(k) = y(k) + y(k - i)
            end do
        end do
	end subroutine

	pure subroutine postfix_sum_i16(y, x)
        integer(kind = selected_int_kind(16)), intent(inout)  	:: y(:)
        integer(kind = selected_int_kind(16)), intent(in)  	:: x(:)

		call prefix_sum(y, x(size(x) : 1 : -1))
	end subroutine

	pure subroutine prefix_sum_i16_sequential(y, x)
        integer(kind = selected_int_kind(16)), intent(inout)  	:: y(:)
        integer(kind = selected_int_kind(16)), intent(in)  	    :: x(:)
        integer                 	                            :: i

        if (size(x) == 0) then
            return
        end if

        assert_pure(size(y) .eq. size(x))

        y(1) = x(1)

        do i = 2, size(x)
            y(i) = y(i - 1) + x(i)
        end do
	end subroutine

	pure subroutine prefix_sum_i16_vectorized(y, x)
        integer(kind = selected_int_kind(16)), intent(inout)  	:: y(:)
        integer(kind = selected_int_kind(16)), intent(in)  	    :: x(:)
        integer(kind = selected_int_kind(16))		            :: tmp
        integer                	                                :: i, k, n

        if (size(x) == 0) then
            return
        end if

        assert_pure(size(y) .eq. size(x))

        n = size(x)

        !Copy array
        y = x

        !Up-Sweep phase
        i = 1
        do while (i .le. n/2)
            do k = 2 * i, n, 2 * i
                y(k) = y(k) + y(k - i)
            end do

            i = 2 * i
        end do

        !Down-Sweep phase
        do while (i .ge. 2)
            i = i / 2

            do k = 3 * i, n, 2 * i
                y(k) = y(k) + y(k - i)
            end do
        end do
	end subroutine

	pure subroutine postfix_sum_r32(y, x)
        real, intent(inout)	:: y(:)
        real, intent(in)	:: x(:)

		call prefix_sum(y, x(size(x) : 1 : -1))
	end subroutine

	pure subroutine prefix_sum_r32_sequential(y, x)
        real, intent(inout)	:: y(:)
        real, intent(in)	:: x(:)
        integer             :: i

        if (size(x) == 0) then
            return
        end if

        assert_pure(size(y) .eq. size(x))

        y(1) = x(1)

        do i = 2, size(x)
            y(i) = y(i - 1) + x(i)
        end do
	end subroutine

	pure subroutine prefix_sum_r32_vectorized(y, x)
        real, intent(inout)	:: y(:)
        real, intent(in)	:: x(:)
        real		        :: tmp
        integer             :: i, k, n

        if (size(x) == 0) then
            return
        end if

        assert_pure(size(y) .eq. size(x))

        n = size(x)

        !Copy array
        y = x

        !Up-Sweep phase
        i = 1
        do while (i .le. n/2)
            do k = 2 * i, n, 2 * i
                y(k) = y(k) + y(k - i)
            end do

            i = 2 * i
        end do

        !Down-Sweep phase
        do while (i .ge. 2)
            i = i / 2

            do k = 3 * i, n, 2 * i
                y(k) = y(k) + y(k - i)
            end do
        end do
	end subroutine

	pure subroutine postfix_sum_r64(y, x)
        double precision, intent(inout)	:: y(:)
        double precision, intent(in)	:: x(:)

		call prefix_sum(y, x(size(x) : 1 : -1))
	end subroutine

	pure subroutine prefix_sum_r64_sequential(y, x)
        double precision, intent(inout)	:: y(:)
        double precision, intent(in)	:: x(:)
        integer                         :: i

        if (size(x) == 0) then
            return
        end if

        assert_pure(size(y) .eq. size(x))

        y(1) = x(1)

        do i = 2, size(x)
            y(i) = y(i - 1) + x(i)
        end do
	end subroutine

	pure subroutine prefix_sum_r64_vectorized(y, x)
        double precision, intent(inout)	:: y(:)
        double precision, intent(in)	:: x(:)
        double precision		        :: tmp
        integer                         :: i, k, n

        if (size(x) == 0) then
            return
        end if

        assert_pure(size(y) .eq. size(x))

        n = size(x)

        !Copy array
        y = x

        !Up-Sweep phase
        i = 1
        do while (i .le. n/2)
            do k = 2 * i, n, 2 * i
                y(k) = y(k) + y(k - i)
            end do

            i = 2 * i
        end do

        !Down-Sweep phase
        do while (i .ge. 2)
            i = i / 2

            do k = 3 * i, n, 2 * i
                y(k) = y(k) + y(k - i)
            end do
        end do
	end subroutine

    pure subroutine postfix_sum_r128(y, x)
        real (kind = QP), intent(inout)	:: y(:)
        real (kind = QP), intent(in)	:: x(:)

		call prefix_sum(y, x(size(x) : 1 : -1))
	end subroutine

	pure subroutine prefix_sum_r128_sequential(y, x)
        real (kind = QP), intent(inout)	:: y(:)
        real (kind = QP), intent(in)	:: x(:)
        integer                         :: i

        if (size(x) == 0) then
            return
        end if

        assert_pure(size(y) .eq. size(x))

        y(1) = x(1)

        do i = 2, size(x)
            y(i) = y(i - 1) + x(i)
        end do
	end subroutine

	pure subroutine prefix_sum_r128_vectorized(y, x)
        real (kind = QP), intent(inout)	:: y(:)
        real (kind = QP), intent(in)	:: x(:)
        real (kind = QP)		        :: tmp
        integer                         :: i, k, n

        if (size(x) == 0) then
            return
        end if

        assert_pure(size(y) .eq. size(x))

        n = size(x)

        !Copy array
        y = x

        !Up-Sweep phase
        i = 1
        do while (i .le. n/2)
            do k = 2 * i, n, 2 * i
                y(k) = y(k) + y(k - i)
            end do

            i = 2 * i
        end do

        !Down-Sweep phase
        do while (i .ge. 2)
            i = i / 2

            do k = 3 * i, n, 2 * i
                y(k) = y(k) + y(k - i)
            end do
        end do
	end subroutine

	subroutine reduce_r32_mpi(s, mpi_op)
        real, intent(inout)		        :: s
		integer, intent(in)		        :: mpi_op

		integer					        :: i_error

#		if defined(_MPI)
            call mpi_allreduce(MPI_IN_PLACE, s, 1, MPI_REAL, mpi_op, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
#		endif
	end subroutine

	subroutine reduce_r32(s, v, mpi_op, global)
        real, intent(inout)		        :: s
        real, intent(in)		        :: v(:)
		integer, intent(in)		        :: mpi_op
		logical, intent(in)             :: global

        select case (mpi_op)
            case (MPI_MAX)
                s = maxval(v)
            case (MPI_MIN)
                s = minval(v)
            case (MPI_SUM)
                s = sum(v)
            case (MPI_PROD)
                s = product(v)
            case default
                assert(.false.)
        end select

#		if defined(_MPI)
            if (global) then
                call reduce_r32_mpi(s, mpi_op)
            end if
#		endif
	end subroutine

	subroutine reduce_r64_mpi(s, mpi_op)
        double precision, intent(inout)         :: s
		integer, intent(in)			            :: mpi_op

		integer					                :: i_error

#		if defined(_MPI)
            call mpi_allreduce(MPI_IN_PLACE, s, 1, MPI_DOUBLE_PRECISION, mpi_op, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
#		endif
	end subroutine

	subroutine reduce_r64(s, v, mpi_op, global)
        double precision, intent(inout)         :: s
        double precision, intent(in)            :: v(:)
		integer, intent(in)			            :: mpi_op
		logical, intent(in)                     :: global

        select case (mpi_op)
            case (MPI_MAX)
                s = maxval(v)
            case (MPI_MIN)
                s = minval(v)
            case (MPI_SUM)
                s = sum(v)
            case (MPI_PROD)
                s = product(v)
            case default
                assert(.false.)
        end select

#		if defined(_MPI)
            if (global) then
                call reduce_r64_mpi(s, mpi_op)
            end if
#		endif
	end subroutine

	subroutine reduce_r128(s, v, mpi_op, global)
        real (kind = QP), intent(inout)         :: s
        real (kind = QP), intent(in)            :: v(:)
		integer, intent(in)			            :: mpi_op
		logical, intent(in)                     :: global

        select case (mpi_op)
            case (MPI_MAX)
                s = maxval(v)
            case (MPI_MIN)
                s = minval(v)
            case (MPI_SUM)
                s = sum(v)
            case (MPI_PROD)
                s = product(v)
            case default
                assert(.false.)
        end select

#		if defined(_MPI)
            if (global) then
                call reduce_r128_mpi(s, mpi_op)
            end if
#		endif
	end subroutine

	subroutine reduce_r128_mpi(s, mpi_op)
        real (kind = QP), intent(inout)         :: s
		integer, intent(in)			            :: mpi_op

		integer					                :: i_error

#		if defined(_MPI)
            call mpi_allreduce(MPI_IN_PLACE, s, 1, MPI_REAL16, mpi_op, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
#		endif
	end subroutine

	subroutine reduce_i8_mpi(s, mpi_op)
      	integer*1, intent(inout)	        :: s
		integer, intent(in)		            :: mpi_op

		integer					            :: i_error

#		if defined(_MPI)
            call mpi_allreduce(MPI_IN_PLACE, s, 1, MPI_INTEGER1, mpi_op, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
#		endif
	end subroutine

	subroutine reduce_i8(s, v, mpi_op, global)
      	integer*1, intent(inout)	        :: s
        integer*1, intent(in)		        :: v(:)
		integer, intent(in)		            :: mpi_op
		logical, intent(in)                 :: global

        select case (mpi_op)
            case (MPI_MAX)
                s = maxval(v)
            case (MPI_MIN)
                s = minval(v)
            case (MPI_SUM)
                s = sum(v)
            case (MPI_PROD)
                s = product(v)
            case default
                assert(.false.)
        end select

#		if defined(_MPI)
            if (global) then
                call reduce_i8_mpi(s, mpi_op)
            end if
#		endif
	end subroutine

	subroutine reduce_i16_mpi(s, mpi_op)
      	integer*2, intent(inout)	        :: s
		integer, intent(in)			        :: mpi_op

		integer					            :: i_error

#		if defined(_MPI)
            call mpi_allreduce(MPI_IN_PLACE, s, 1, MPI_INTEGER2, mpi_op, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
#		endif
	end subroutine

	subroutine reduce_i16(s, v, mpi_op, global)
      	integer*2, intent(inout)	        :: s
        integer*2, intent(in)		        :: v(:)
		integer, intent(in)			        :: mpi_op
		logical, intent(in)                 :: global

        select case (mpi_op)
            case (MPI_MAX)
                s = maxval(v)
            case (MPI_MIN)
                s = minval(v)
            case (MPI_SUM)
                s = sum(v)
            case (MPI_PROD)
                s = product(v)
            case default
                assert(.false.)
        end select

#		if defined(_MPI)
            if (global) then
                call reduce_i16_mpi(s, mpi_op)
            end if
#		endif
	end subroutine

	subroutine reduce_i32_mpi(s, mpi_op)
      	integer*4, intent(inout)	        :: s
		integer, intent(in)			        :: mpi_op

		integer					            :: i_error

#		if defined(_MPI)
            call mpi_allreduce(MPI_IN_PLACE, s, 1, MPI_INTEGER4, mpi_op, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
#		endif
	end subroutine

	subroutine reduce_i32(s, v, mpi_op, global)
      	integer*4, intent(inout)	        :: s
        integer*4, intent(in)		        :: v(:)
		integer, intent(in)			        :: mpi_op
		logical, intent(in)                 :: global

        select case (mpi_op)
            case (MPI_MAX)
                s = maxval(v)
            case (MPI_MIN)
                s = minval(v)
            case (MPI_SUM)
                s = sum(v)
            case (MPI_PROD)
                s = product(v)
            case default
                assert(.false.)
        end select

#		if defined(_MPI)
            if (global) then
                call reduce_i32_mpi(s, mpi_op)
            end if
#		endif
	end subroutine

	subroutine reduce_i64_mpi(s, mpi_op)
      	integer*8, intent(inout)	        :: s
		integer, intent(in)			        :: mpi_op

		integer					            :: i_error

#		if defined(_MPI)
            call mpi_allreduce(MPI_IN_PLACE, s, 1, MPI_INTEGER8, mpi_op, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
#		endif
	end subroutine

	subroutine reduce_i64(s, v, mpi_op, global)
      	integer*8, intent(inout)	        :: s
        integer*8, intent(in)		        :: v(:)
		integer, intent(in)			        :: mpi_op
		logical, intent(in)                 :: global

        select case (mpi_op)
            case (MPI_MAX)
                s = maxval(v)
            case (MPI_MIN)
                s = minval(v)
            case (MPI_SUM)
                s = sum(v)
            case (MPI_PROD)
                s = product(v)
            case default
                assert(.false.)
        end select

#		if defined(_MPI)
            if (global) then
                call reduce_i64_mpi(s, mpi_op)
            end if
#		endif
	end subroutine

	subroutine reduce_l_mpi(s, mpi_op)
        logical, intent(inout)          :: s
		integer, intent(in)		        :: mpi_op

		integer					        :: i_error

#		if defined(_MPI)
            call mpi_allreduce(MPI_IN_PLACE, s, 1, MPI_LOGICAL, mpi_op, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
#		endif
	end subroutine

	subroutine reduce_l(s, v, mpi_op, global)
        logical, intent(inout)          :: s
        logical, intent(in)             :: v(:)
		integer, intent(in)		        :: mpi_op
		logical, intent(in)             :: global

        select case (mpi_op)
            case (MPI_LAND)
                s = all(v)
            case (MPI_LOR)
                s = any(v)
            case default
                assert(.false.)
        end select

#		if defined(_MPI)
            if (global) then
                call reduce_l_mpi(s, mpi_op)
            end if
#		endif
	end subroutine

	subroutine scatter_r32(s, v)
        real, intent(in)    :: s
        real, intent(out)	:: v(:)

        v = s
	end subroutine

	subroutine scatter_r64(s, v)
        double precision, intent(in)		:: s
        double precision, intent(out)	:: v(:)

        v = s
	end subroutine

	subroutine scatter_r128(s, v)
        real (kind = QP), intent(in)		:: s
        real (kind = QP), intent(out)	:: v(:)

        v = s
	end subroutine

	subroutine scatter_i32(s, v)
        integer, intent(in)		:: s
        integer, intent(out)		:: v(:)

        v = s
	end subroutine

	subroutine scatter_i64(s, v)
        integer*2, intent(in)		:: s
        integer*2, intent(out)		:: v(:)

        v = s
	end subroutine

	subroutine scatter_l(s, v)
        logical, intent(in)		:: s
        logical, intent(out)	:: v(:)

        v = s
	end subroutine

	subroutine scatter_s(s, v)
        character(*), intent(in)		:: s
        character(*), intent(out)	    :: v(:)

        v = s
	end subroutine

	pure function term_color(i_color) result(str)
        integer, intent(in) :: i_color
        character(5) :: str

        character(5), parameter :: sequences(0:7) = &
            [CHAR(27) // "[97m", &
            CHAR(27) // "[91m", &
            CHAR(27) // "[92m", &
            CHAR(27) // "[93m", &
            CHAR(27) // "[94m", &
            CHAR(27) // "[95m", &
            CHAR(27) // "[96m", &
            CHAR(27) // "[90m"]

        str = sequences(iand(i_color, 7))
    end function

    pure function term_reset() result(str)
        character(5) :: str

        character(5), parameter :: seq = CHAR(27) // "[0m"

        str = seq
    end function

    pure subroutine raise_error()
        integer :: i

        i = 0
        i = i / i
    end subroutine
end MODULE


