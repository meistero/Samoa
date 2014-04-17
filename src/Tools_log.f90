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

	integer 		    :: rank_MPI = 0
	integer 		    :: size_MPI = 1
	integer             :: ref_count_MPI = 0

	contains

    !> initializes the mpi communicator
    subroutine init_mpi()
        integer                             :: i_error, mpi_prov_thread_support
        integer(kind = MPI_ADDRESS_KIND)    :: mpi_tag_upper_bound
        logical                             :: mpi_flag, mpi_is_initialized

#       if defined(_MPI)
            if (ref_count_MPI == 0) then
                call mpi_initialized(mpi_is_initialized, i_error); assert_eq(i_error, 0)

                if (.not. mpi_is_initialized) then
#                   if defined(_OPENMP)
                        call mpi_init_thread(MPI_THREAD_MULTIPLE, mpi_prov_thread_support, i_error); assert_eq(i_error, 0)

						try(mpi_prov_thread_support >= MPI_THREAD_MULTIPLE, "MPI version does not support MPI_THREAD_MULTIPLE")
#                   else
                        call mpi_init(i_error); assert_eq(i_error, 0)
#                   endif
                else
                    !Set MPI reference counter to 1 if samoa did not initialize MPI
                    !This makes sure that it will not finalize MPI either.
                    ref_count_MPI = 1

#                   if defined(_OPENMP)
                        call mpi_query_thread(mpi_prov_thread_support, i_error); assert_eq(i_error, 0)

						try(mpi_prov_thread_support >= MPI_THREAD_MULTIPLE, "MPI version does not support MPI_THREAD_MULTIPLE")
#                   endif
                end if

                call mpi_comm_size(MPI_COMM_WORLD, size_MPI, i_error); assert_eq(i_error, 0)
                call mpi_comm_rank(MPI_COMM_WORLD, rank_MPI, i_error); assert_eq(i_error, 0)

                call mpi_comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, mpi_tag_upper_bound, mpi_flag, i_error); assert_eq(i_error, 0)
				try(mpi_flag, "MPI tag support could not be determined")
				try(mpi_tag_upper_bound >= ishft(1, 30) - 1, "MPI version does not support a sufficient number of tags")
            end if

            ref_count_MPI = ref_count_MPI + 1
#       else
            size_MPI = 1
            rank_MPI = 0
#       endif
    end subroutine

    !> finalizes the mpi communicator
    subroutine finalize_mpi()
        integer                         :: i_error
        logical                         :: mpi_is_finalized

#	    if defined(_MPI)
            ref_count_MPI = ref_count_MPI - 1

            if (ref_count_MPI == 0) then
                call mpi_finalized(mpi_is_finalized, i_error); assert_eq(i_error, 0)

                if (.not. mpi_is_finalized) then
                    call mpi_barrier(MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
                    call mpi_finalize(i_error); assert_eq(i_error, 0)
                end if
            end if
#	    endif
    end subroutine
end module

module Tools_openmp
#	if defined(_OPENMP)
        use omp_lib
#	endif

    public

	integer                     :: i_thread = 1     !in contrast to OpenMP our thread index starts with 1
	!$omp threadprivate(i_thread)

#	if !defined(_OPENMP)
        contains

        function omp_get_thread_num() result(i_thread)
            integer :: i_thread

            i_thread = 0
        end function

        function omp_get_num_threads() result(i_threads)
            integer :: i_threads

            i_threads = 1
        end function

        function omp_get_max_threads() result(i_threads)
            integer :: i_threads

            i_threads = 1
        end function

        subroutine omp_set_num_threads(i_threads)
            integer :: i_threads

            assert_eq(i_threads, 1)
        end subroutine
#   endif
end module

#include "Tools_parallel_operators.inc"

MODULE Tools_log
    use Tools_mpi
    use Tools_openmp

	private
	public get_wtime, log_open_file, log_close_file, g_log_file_unit, raise_error, term_color, term_reset


	!> global file unit (there's only one log instance possible due to the lack of variadic functions and macros in Fortran,
	!> which makes it impossible to wrap WRITE and PRINT)
	integer (kind=selected_int_kind(8))     :: g_log_file_unit = 6

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

    function get_wtime_internal() result(time)
        double precision :: time

        integer(kind = selected_int_kind(16)) :: counts, count_rate

        call system_clock(counts, count_rate)

        time = dble(counts) / dble(count_rate)
    end function

    function get_wtime() result(wtime)
        double precision :: wtime

#       if defined(_OPENMP)
            wtime = omp_get_wtime()
#       elif defined(_MPI)
            wtime = MPI_wtime()
#       else
            wtime = get_wtime_internal()
#       endif
    end function
end MODULE
