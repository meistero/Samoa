! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

module Tools_mpi
#	if defined(_MPI)
#       if defined(__GFORTRAN__)
            use mpi
#       else
            include 'mpif.h'
#       endif
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

				try(mpi_tag_upper_bound >= ishft(1, 28) - 1, "MPI version does not support a sufficient number of tags")
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

MODULE Tools_log
    use Tools_mpi
    use Tools_openmp

	private
	public get_wtime, log_open_file, log_close_file, g_log_file_unit, get_free_file_unit, raise_error, term_color, term_reset, time_to_hrt


	!> global file unit (there's only one log instance possible due to the lack of variadic functions and macros in Fortran,
	!> which makes it impossible to wrap WRITE and PRINT)
	integer (kind=selected_int_kind(8))     :: g_log_file_unit = 6

	contains

    !> Looks for an available unit that can be used to open a file
	function get_free_file_unit() result(file_unit)
		integer 					:: i_error
		logical						:: l_is_open
		integer                     :: file_unit

		do file_unit = 10, huge(1)
			inquire(unit=file_unit, opened=l_is_open, iostat=i_error)

			if (i_error == 0 .and. .not. l_is_open) then
                exit
			end if
		end do
	end function

	!> Opens an external log file where the log write commands are written to,
	!> if this is not called, the log is written to stdout instead
	subroutine log_open_file(s_file_name)
		character(*)				:: s_file_name

		integer 					:: i_error

		assert_eq(g_log_file_unit, 6)

		g_log_file_unit = get_free_file_unit()

        open(unit=g_log_file_unit, file=s_file_name, status='replace', access='sequential', iostat=i_error); assert_eq(i_error, 0)
        flush(g_log_file_unit)
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
        character(4) :: str

        character(4), parameter :: seq = CHAR(27) // "[0m"

        str = seq
    end function

    pure subroutine raise_error()
        integer :: i

        i = 0
        i = i / i
    end subroutine

    function get_wtime() result(wtime)
        double precision :: wtime

#       if defined(_OPENMP)
            wtime = omp_get_wtime()
#       elif defined(_MPI)
            wtime = MPI_wtime()
#       else
            integer(kind = selected_int_kind(16))   :: counts, count_rate

            call system_clock(counts, count_rate)

            time = dble(counts) / dble(count_rate)
#       endif
    end function

    function time_to_hrt(time) result(str)
        double precision, intent(in)    :: time
        character(256)                  :: str

        character(3), parameter         :: s_unit_names(7) = [character(3) :: "d", "h", "min", "s", "ms", "µs", "ns"]
        integer , parameter             :: i_max_units = 3
        integer (kind = 8)              :: i_unit_values(7)
        integer                         :: i_unit, i_units

        i_unit_values(1) = int(time, kind=8) / (24_8 * 60_8 * 60_8)
        i_unit_values(2) = mod(int(time, kind=8) / (60_8 * 60_8), 24_8)
        i_unit_values(3) = mod(int(time, kind=8) / 60_8, 60_8)
        i_unit_values(4) = mod(int(time, kind=8), 60_8)
        i_unit_values(5) = mod(int(time * 1e3_8, kind=8), 1000_8)
        i_unit_values(6) = mod(int(time * 1e6_8, kind=8), 1000_8)
        i_unit_values(7) = mod(int(time * 1e9_8, kind=8), 1000_8)

        str = ""
        i_units = 0

        do i_unit = 1, 7
            if (i_units < i_max_units .and. i_unit_values(i_unit) > 0) then
                write (str, '(A, " ", I0, A)') trim(str), i_unit_values(i_unit), trim(s_unit_names(i_unit))
                i_units = i_units + 1
            end if
        end do
    end function
end MODULE

#include "Tools_parallel_operators.inc"
