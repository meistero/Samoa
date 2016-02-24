! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE

!> @author Oliver Meister

#if !defined(_T_MPI_TYPE)
#   define _T_CLASS_ARITHMETIC     0
#   define _T_CLASS_LOGICAL        1
#   define _T_CLASS_STRING         2

    module Tools_parallel_operators
        use Tools_mpi
        use Tools_openmp
        use Tools_log

        implicit none

        interface prefix_sum
            module procedure prefix_sum_sequential_i8
            module procedure prefix_sum_sequential_i16
            module procedure prefix_sum_sequential_i32
            module procedure prefix_sum_sequential_i64
            module procedure prefix_sum_sequential_r32
            module procedure prefix_sum_sequential_r64
            module procedure prefix_sum_sequential_r128
            module procedure prefix_sum_sequential_l
            module procedure prefix_sum_mpi_i8
            module procedure prefix_sum_mpi_i16
            module procedure prefix_sum_mpi_i32
            module procedure prefix_sum_mpi_i64
            module procedure prefix_sum_mpi_r32
            module procedure prefix_sum_mpi_r64
            module procedure prefix_sum_mpi_r128
            module procedure prefix_sum_mpi_l
        end interface

        interface reduce
            module procedure reduce_i8
            module procedure reduce_i16
            module procedure reduce_i32
            module procedure reduce_i64
            module procedure reduce_r32
            module procedure reduce_r64
            module procedure reduce_r128
            module procedure reduce_l
            module procedure reduce_mpi_i8
            module procedure reduce_mpi_i16
            module procedure reduce_mpi_i32
            module procedure reduce_mpi_i64
            module procedure reduce_mpi_r32
            module procedure reduce_mpi_r64
            module procedure reduce_mpi_r128
            module procedure reduce_mpi_vector_i8
            module procedure reduce_mpi_vector_i16
            module procedure reduce_mpi_vector_i32
            module procedure reduce_mpi_vector_i64
            module procedure reduce_mpi_vector_r32
            module procedure reduce_mpi_vector_r64
            module procedure reduce_mpi_vector_r128
            module procedure reduce_mpi_l
            module procedure reduce_mpi_vector_l
        end interface

        interface scatter
            module procedure scatter_i8
            module procedure scatter_i16
            module procedure scatter_i32
            module procedure scatter_i64
            module procedure scatter_r32
            module procedure scatter_r64
            module procedure scatter_r128
            module procedure scatter_l
            module procedure scatter_s
        end interface

        integer, parameter :: QP = kind(1.0q0)

        public prefix_sum, reduce, scatter

        contains

#       define _T               integer (kind = selected_int_kind(1))
#       define _T_MPI_TYPE      MPI_INTEGER1
#       define _T_SUFFIX        i8
#       define _T_SIZE          1
#       define _T_CLASS         _T_CLASS_ARITHMETIC
#       include "Tools_parallel_operators.inc"

#       define _T               integer (kind = selected_int_kind(4))
#       define _T_MPI_TYPE      MPI_INTEGER2
#       define _T_SUFFIX        i16
#       define _T_SIZE          2
#       define _T_CLASS         _T_CLASS_ARITHMETIC
#       include "Tools_parallel_operators.inc"

#       define _T               integer (kind = selected_int_kind(8))
#       define _T_MPI_TYPE      MPI_INTEGER4
#       define _T_SUFFIX        i32
#       define _T_SIZE          4
#       define _T_CLASS         _T_CLASS_ARITHMETIC
#       include "Tools_parallel_operators.inc"

#       define _T               integer (kind = selected_int_kind(16))
#       define _T_MPI_TYPE      MPI_INTEGER8
#       define _T_SUFFIX        i64
#       define _T_SIZE          8
#       define _T_CLASS         _T_CLASS_ARITHMETIC
#       include "Tools_parallel_operators.inc"

#       define _T               real
#       define _T_MPI_TYPE      MPI_REAL
#       define _T_SUFFIX        r32
#       define _T_SIZE          4
#       define _T_CLASS         _T_CLASS_ARITHMETIC
#       include "Tools_parallel_operators.inc"

#       define _T               double precision
#       define _T_MPI_TYPE      MPI_DOUBLE_PRECISION
#       define _T_SUFFIX        r64
#       define _T_SIZE          8
#       define _T_CLASS         _T_CLASS_ARITHMETIC
#       include "Tools_parallel_operators.inc"

#       define _T               real(kind = QP)
#       define _T_MPI_TYPE      MPI_REAL16
#       define _T_SUFFIX        r128
#       define _T_SIZE          16
#       define _T_CLASS         _T_CLASS_ARITHMETIC
#       include "Tools_parallel_operators.inc"

#       define _T               logical
#       define _T_MPI_TYPE      MPI_LOGICAL
#       define _T_SUFFIX        l
#       define _T_SIZE          4
#       define _T_CLASS         _T_CLASS_LOGICAL
#       include "Tools_parallel_operators.inc"

#       define _T               character(*)
#       define _T_MPI_TYPE      MPI_CHARACTER
#       define _T_SUFFIX        s
#       define _T_SIZE          -1
#       define _T_CLASS         _T_CLASS_STRING
#       include "Tools_parallel_operators.inc"
    end module
#else

#define _PREFIX(P, X)   _conc3(P,_,X)
#define _T_(x)          _PREFIX(x,_T_SUFFIX)

pure subroutine _T_(prefix_sum_sequential)(y, x)
    _T, intent(inout)	:: y(:)
    _T, intent(in)	    :: x(:)

#   if _T_CLASS == _T_CLASS_ARITHMETIC
        integer         :: i

        if (size(x) == 0) then
            return
        end if

        assert_pure(size(y) .eq. size(x))

        y(1) = x(1)

        do i = 2, size(x)
            y(i) = y(i - 1) + x(i)
        end do
#   else
        assert_pure(.false.)
#   endif
end subroutine

pure subroutine _T_(prefix_sum_vectorized)(y, x)
    _T, intent(inout)  	:: y(:)
    _T, intent(in)      :: x(:)

#   if _T_CLASS == _T_CLASS_ARITHMETIC
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
#   else
       assert_pure(.false.)
#   endif
end subroutine

subroutine _T_(prefix_sum_mpi)(s, mpi_op)
    _T, intent(inout)	        :: s
    integer, intent(in)		    :: mpi_op

#   if defined(_MPI) &&_T_CLASS != _T_CLASS_STRING
        integer				        :: i_error

#       if defined(__GFORTRAN__)
            _T	        :: r

            !HOTFIX: for some reason, MPI_IN_PLACE does not work with gfortran 4.8 and mpich
            call mpi_scan(s, r, 1, _T_MPI_TYPE, mpi_op, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
            s = r
#       else
            call mpi_scan(MPI_IN_PLACE, s, 1, _T_MPI_TYPE, mpi_op, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
#       endif

        assert_eq(sizeof(s), _T_SIZE)
#   endif
end subroutine

subroutine _T_(reduce_local)(s, v, mpi_op)
    _T, intent(inout)	        :: s
    _T, intent(in)		        :: v(:)
    integer, intent(in)		    :: mpi_op

    select case (mpi_op)
#       if _T_CLASS == _T_CLASS_ARITHMETIC
            case (MPI_MAX)
                s = maxval(v)
            case (MPI_MIN)
                s = minval(v)
            case (MPI_SUM)
                s = sum(v)
            case (MPI_PROD)
                s = product(v)
#       elif _T_CLASS == _T_CLASS_LOGICAL
            case (MPI_LAND)
                s = all(v)
            case (MPI_LOR)
                s = any(v)
#       endif

        case default
            assert(.false.)
    end select
end subroutine

subroutine _T_(reduce_mpi)(s, mpi_op)
    _T, intent(inout)	        :: s
    integer, intent(in)		    :: mpi_op

#   if defined(_MPI) &&_T_CLASS != _T_CLASS_STRING
        integer				    :: i_error

#       if defined(__GFORTRAN__)
            _T	        :: r

            !HOTFIX: for some reason, MPI_IN_PLACE does not work with gfortran 4.8 and mpich
            call mpi_allreduce(s, r, 1, _T_MPI_TYPE, mpi_op, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
            s = r
#       else
            call mpi_allreduce(MPI_IN_PLACE, s, 1, _T_MPI_TYPE, mpi_op, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
#       endif

        assert_eq(sizeof(s), _T_SIZE)
#   endif
end subroutine

subroutine _T_(reduce_mpi_vector)(s, mpi_op)
    _T, intent(inout)	        :: s(:)
    integer, intent(in)		    :: mpi_op

#   if defined(_MPI) &&_T_CLASS != _T_CLASS_STRING
        integer				    :: i_error

#       if defined(__GFORTRAN__)
            _T	        :: r(size(s))

            !HOTFIX: for some reason, MPI_IN_PLACE does not work with gfortran 4.8 and mpich
            call mpi_allreduce(s, r, size(s), _T_MPI_TYPE, mpi_op, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
            s = r
#       else
            call mpi_allreduce(MPI_IN_PLACE, s, size(s), _T_MPI_TYPE, mpi_op, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
#       endif

        assert_ge(size(s), 1)
        assert_eq(sizeof(s(1)), _T_SIZE)
#   endif
end subroutine

subroutine _T_(reduce)(s, v, mpi_op, global)
    _T, intent(inout)	        :: s
    _T, intent(in)		        :: v(:)
    integer, intent(in)		    :: mpi_op
    logical, intent(in)         :: global

    call _T_(reduce_local)(s, v, mpi_op)

#   if defined(_MPI)
        if (global) then
            call _T_(reduce_mpi)(s, mpi_op)
        end if
#   endif
end subroutine

subroutine _T_(scatter)(s, v)
    _T, intent(in)		:: s
    _T, intent(out)	    :: v(:)

    v(:) = s
end subroutine

#   undef _T
#   undef _T_
#   undef _T_MPI_TYPE
#   undef _T_SUFFIX
#   undef _T_SIZE
#   undef _T_CLASS
#endif
