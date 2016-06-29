#include "Compilation_control.f90"

module Tools_statistics
    use Config
    use Tools_mpi
    use Tools_log
    use Tools_parallel_operators

    private
    public t_base_statistics, t_statistics, t_adaptive_statistics

    enum, bind(c)
        enumerator :: TRAVERSALS = 1, TRAVERSED_CELLS, TRAVERSED_EDGES, TRAVERSED_NODES, TRAVERSED_MEMORY, STATS_TYPE_END_INT_COUNTER
    end enum

    public :: TRAVERSALS, TRAVERSED_CELLS, TRAVERSED_EDGES, TRAVERSED_NODES, TRAVERSED_MEMORY

    enum, bind(c)
        enumerator :: TRAVERSAL_TIME = 1, PRE_COMPUTE_TIME, INNER_COMPUTE_TIME, POST_COMPUTE_TIME, ASAGI_TIME, SYNC_TIME, BARRIER_TIME, STATS_TYPE_END_TIME_COUNTER
        enumerator :: ALLOCATION_TIME = STATS_TYPE_END_TIME_COUNTER, UPDATE_DISTANCES_TIME, UPDATE_NEIGHBORS_TIME, INTEGRITY_TIME, LOAD_BALANCING_TIME, STATS_TYPE_END_ADAPTIVE_COUNTER
    end enum

    public :: TRAVERSAL_TIME, PRE_COMPUTE_TIME, INNER_COMPUTE_TIME, POST_COMPUTE_TIME, ASAGI_TIME, SYNC_TIME, BARRIER_TIME
    public :: ALLOCATION_TIME, UPDATE_DISTANCES_TIME, UPDATE_NEIGHBORS_TIME, INTEGRITY_TIME, LOAD_BALANCING_TIME

	type t_base_statistics
        private
        integer (kind = selected_int_kind(16)) :: counter_stats(STATS_TYPE_END_INT_COUNTER - 1) = 0_8

		contains

		procedure, private, pass :: assign_base_stats => t_base_statistics_assign

		procedure, public, pass :: add_counter => t_base_statistics_add_counter
		procedure, public, pass :: get_counter => t_base_statistics_get_counter
		procedure, public, pass :: clear_counter => t_base_statistics_clear_counter
	end type

	type, extends(t_base_statistics) :: t_statistics
        private
        double precision :: time_stats(STATS_TYPE_END_TIME_COUNTER - 1) = 0.d0

		contains

		procedure, private, pass :: reduce_local => t_statistics_reduce_local
		procedure, private, pass :: reduce_global => t_statistics_reduce_global
		procedure, private, pass :: add => t_statistics_add
		procedure, private, pass :: negate => t_statistics_negate
		procedure, private, pass :: sub => t_statistics_sub
		procedure, private, pass :: assign => t_statistics_assign

		procedure, public, pass :: start_time => t_statistics_start_time
		procedure, public, pass :: stop_time => t_statistics_stop_time
		procedure, public, pass :: get_time => t_statistics_get_time
		procedure, public, pass :: clear_time => t_statistics_clear_time

        procedure, public, pass :: scale => t_statistics_scale
        procedure, public, pass :: to_string => t_statistics_to_string
		procedure, public, pass :: clear => t_statistics_clear

		generic, public :: reduce => reduce_local, reduce_global

		generic, public :: assignment(=) => assign
		generic, public :: operator(+) => add
		generic, public :: operator(-) => sub, negate
	end type

	type, extends(t_base_statistics) ::  t_adaptive_statistics
        private
        double precision :: adaptive_stats(STATS_TYPE_END_ADAPTIVE_COUNTER - 1) = 0.d0

		contains

		procedure, private, pass :: reduce_local => t_adaptive_statistics_reduce_local
		procedure, private, pass :: reduce_global => t_adaptive_statistics_reduce_global
		procedure, private, pass :: add_adaptive_stats => t_adaptive_statistics_add_adaptive_stats
		procedure, private, pass :: add_stats => t_adaptive_statistics_add_stats
		procedure, private, pass :: assign => t_adaptive_statistics_assign
		procedure, private, pass :: negate => t_adaptive_statistics_negate
		procedure, private, pass :: sub => t_adaptive_statistics_sub

		procedure, public, pass :: start_time => t_adaptive_statistics_start_time
		procedure, public, pass :: stop_time => t_adaptive_statistics_stop_time
		procedure, public, pass :: get_time => t_adaptive_statistics_get_time
		procedure, public, pass :: clear_time => t_adaptive_statistics_clear_time
        procedure, public, pass :: scale => t_adaptive_statistics_scale
        procedure, public, pass :: to_string => t_adaptive_statistics_to_string
		procedure, public, pass :: clear => t_adaptive_statistics_clear

		generic, private :: add => add_stats, add_adaptive_stats
		generic, public :: reduce => reduce_local, reduce_global

        generic, public :: assignment(=) => assign
		generic, public :: operator(+) => add_stats, add_adaptive_stats
		generic, public :: operator(-) => sub, negate
	end type

	contains

    elemental subroutine t_base_statistics_add_counter(s, counter_type, counter_value)
        class(t_base_statistics), intent(inout)			    :: s
        integer, intent(in)			                        :: counter_type
        integer (kind = selected_int_kind(16)), intent(in)  :: counter_value

        assert_pure(counter_type > 0)
        assert_pure(counter_type < STATS_TYPE_END_INT_COUNTER)

		s%counter_stats(counter_type) = s%counter_stats(counter_type) + counter_value
    end subroutine

    elemental function t_base_statistics_get_counter(s, counter_type) result(counter_value)
        class(t_base_statistics), intent(in)	:: s
        integer, intent(in)			            :: counter_type
        integer (kind = selected_int_kind(16))  :: counter_value

        assert_pure(counter_type > 0)
        assert_pure(counter_type < STATS_TYPE_END_INT_COUNTER)

		counter_value = s%counter_stats(counter_type)
    end function

    elemental subroutine t_base_statistics_clear_counter(s, counter_type)
        class(t_base_statistics), intent(inout)			    :: s
        integer, intent(in)			                        :: counter_type

        assert_pure(counter_type > 0)
        assert_pure(counter_type < STATS_TYPE_END_INT_COUNTER)

		s%counter_stats(counter_type) = 0
    end subroutine

    elemental subroutine t_base_statistics_assign(s, s2)
        class(t_base_statistics), intent(inout)			:: s
        type(t_base_statistics), intent(in)			    :: s2

		s%counter_stats = s2%counter_stats
    end subroutine

    subroutine t_statistics_reduce_local(s, v, mpi_op)
        class(t_statistics), intent(inout)	:: s
        type(t_statistics), intent(in)		:: v(:)
        integer, intent(in)                 :: mpi_op

        integer :: i
        integer (kind = selected_int_kind(16))          :: counters(size(v))
        double precision                                :: times(size(v))

        do i = 1, size(s%counter_stats)
            !HACK: Intel Fortran Compiler 13.1 produces wrong results for the following code:
            !call reduce(s%counter_stats(i), v%counter_stats(i), mpi_op, .false.)
            !hence, a workaround with a temporary array:

            counters = v%counter_stats(i)
            call reduce(s%counter_stats(i), counters, mpi_op, .false.)
		end do

        do i = 1, size(s%time_stats)
            !HACK: same here

            times = v%time_stats(i)
            call reduce(s%time_stats(i), times, mpi_op, .false.)
		end do
     end subroutine

    subroutine t_statistics_reduce_global(s, mpi_op)
        class(t_statistics), intent(inout)          :: s
        integer, intent(in)                         :: mpi_op

		call reduce(s%counter_stats, mpi_op)
		call reduce(s%time_stats, mpi_op)
     end subroutine

    elemental subroutine t_statistics_clear(s)
        class(t_statistics), intent(inout)	    :: s

		s = t_statistics()
    end subroutine

    pure subroutine t_statistics_assign(sr, s2)
        class(t_statistics), intent(inout)	    :: sr
        type(t_statistics), intent(in)			:: s2

		sr%counter_stats = s2%counter_stats
		sr%time_stats = s2%time_stats
    end subroutine

    pure function t_statistics_add(s1, s2) result(sr)
        class(t_statistics), intent(in)			:: s1
        type(t_statistics), intent(in)			:: s2
        type(t_statistics)						:: sr

		sr%counter_stats = s1%counter_stats + s2%counter_stats
		sr%time_stats = s1%time_stats + s2%time_stats
    end function

    pure function t_statistics_negate(s1) result(sr)
        class(t_statistics), intent(in)			:: s1
        type(t_statistics)						:: sr

		sr%counter_stats = -s1%counter_stats
		sr%time_stats = -s1%time_stats
    end function

    pure function t_statistics_sub(s1, s2) result(sr)
        class(t_statistics), intent(in)			:: s1
        type(t_statistics), intent(in)			:: s2
        type(t_statistics)						:: sr

		sr = s1%add(s2%negate())
    end function

    pure function t_statistics_scale(s, scaling) result(sr)
        class(t_statistics), intent(in)		            :: s
        double precision, intent(in)                    :: scaling
        type(t_statistics)		                        :: sr

		sr%counter_stats = int(s%counter_stats * scaling, selected_int_kind(16))
		sr%time_stats = s%time_stats * scaling
    end function

    subroutine t_statistics_start_time(s, time_type)
        class(t_statistics), intent(inout)			        :: s
        integer, intent(in)			                        :: time_type

        assert_pure(time_type > 0)
        assert_pure(time_type < STATS_TYPE_END_TIME_COUNTER)

		s%time_stats(time_type) = s%time_stats(time_type) - get_wtime()
    end subroutine

    subroutine t_statistics_stop_time(s, time_type)
        class(t_statistics), intent(inout)			        :: s
        integer, intent(in)			                        :: time_type

        assert_pure(time_type > 0)
        assert_pure(time_type < STATS_TYPE_END_TIME_COUNTER)

		s%time_stats(time_type) = s%time_stats(time_type) + get_wtime()
    end subroutine

    elemental function t_statistics_get_time(s, time_type) result(time_value)
        class(t_statistics), intent(in)			:: s
        integer, intent(in)			            :: time_type
        double precision                        :: time_value

        assert_pure(time_type > 0)
        assert_pure(time_type < STATS_TYPE_END_TIME_COUNTER)

		time_value = s%time_stats(time_type)
    end function

    elemental subroutine t_statistics_clear_time(s, time_type)
        class(t_statistics), intent(inout)			        :: s
        integer, intent(in)			                        :: time_type

        assert_pure(time_type > 0)
        assert_pure(time_type < STATS_TYPE_END_TIME_COUNTER)

		s%time_stats(time_type) = 0.d0
    end subroutine

    pure function t_statistics_to_string(s) result(str)
        class(t_statistics), intent(in)			:: s
		character (len = 512)					:: str

        write(str, '("#travs: ", I0, " time: ", F0.4, " s (comp: ", F0.4, " s (pre: ", F0.4, " s inner: ", F0.4, " s post: ", F0.4, " s) asagi: ", F0.4, " s sync: ", F0.4, " s barr: ", F0.4, " s)")') &
            s%get_counter(traversals), s%get_time(traversal_time),  &
            s%get_time(pre_compute_time) + s%get_time(inner_compute_time) + s%get_time(post_compute_time), &
            s%get_time(pre_compute_time), s%get_time(inner_compute_time), s%get_time(post_compute_time), &
            s%get_time(asagi_time), s%get_time(sync_time), s%get_time(barrier_time)

        if (s%get_time(traversal_time) > 0.0d0) then
            write(str, '(A, " ET: ", F0.4, " M/s  MT: ", F0.4, " GB/s")') trim(str), dble(s%get_counter(traversed_cells)) / (1.0d6 * s%get_time(traversal_time)), &
            dble(s%get_counter(traversed_memory)) / (1024.0d0 * 1024.0d0 * 1024.0d0 * s%get_time(traversal_time))
        else
            write(str, '(A, " ET: ", F0.4, " M/s  MT: ", F0.4, " GB/s")') trim(str), -1.0d0, -1.0d0
        end if
	end function

    subroutine t_adaptive_statistics_reduce_local(s, v, mpi_op)
        class(t_adaptive_statistics), intent(inout)	    :: s
        type(t_adaptive_statistics), intent(in)		    :: v(:)
        integer, intent(in)                             :: mpi_op

        integer :: i
        integer (kind = selected_int_kind(16))          :: counters(size(v))
        double precision                                :: times(size(v))

        do i = 1, size(s%counter_stats)
            !HACK: Intel Fortran Compiler 13.1 produces wrong results for the following code:
            !call reduce(s%counter_stats(i), v%counter_stats(i), mpi_op, .false.)
            !hence, a workaround with a temporary array:

            counters = v%counter_stats(i)
            call reduce(s%counter_stats(i), counters, mpi_op, .false.)
		end do

        do i = 1, size(s%adaptive_stats)
            !HACK: same here

            times = v%adaptive_stats(i)
            call reduce(s%adaptive_stats(i), times, mpi_op, .false.)
		end do
    end subroutine

    subroutine t_adaptive_statistics_reduce_global(s, mpi_op)
        class(t_adaptive_statistics), intent(inout)     :: s
        integer, intent(in)                             :: mpi_op

		call reduce(s%counter_stats, mpi_op)
		call reduce(s%adaptive_stats, mpi_op)
    end subroutine

    elemental subroutine t_adaptive_statistics_clear(s)
        class(t_adaptive_statistics), intent(inout)	    :: s

		s = t_adaptive_statistics()
    end subroutine

    pure subroutine t_adaptive_statistics_assign(sr, s2)
        class(t_adaptive_statistics), intent(inout)	    :: sr
        type(t_adaptive_statistics), intent(in)			:: s2

        sr%counter_stats = s2%counter_stats
		sr%adaptive_stats = s2%adaptive_stats
    end subroutine

    pure function t_adaptive_statistics_add_adaptive_stats(s1, s2) result(sr)
        class(t_adaptive_statistics), intent(in)		:: s1
        type(t_adaptive_statistics), intent(in)		    :: s2
        type(t_adaptive_statistics)						:: sr

		sr%counter_stats = s1%counter_stats + s2%counter_stats
		sr%adaptive_stats = s1%adaptive_stats + s2%adaptive_stats
    end function

    pure function t_adaptive_statistics_add_stats(s1, s2) result(sr)
        class(t_adaptive_statistics), intent(in)		:: s1
        type(t_statistics), intent(in)		    :: s2
        type(t_adaptive_statistics)						:: sr

		sr%counter_stats = s1%counter_stats + s2%counter_stats
		sr%adaptive_stats(:STATS_TYPE_END_TIME_COUNTER - 1) = s1%adaptive_stats(:STATS_TYPE_END_TIME_COUNTER - 1) + s2%time_stats
    end function

    pure function t_adaptive_statistics_negate(s1) result(sr)
        class(t_adaptive_statistics), intent(in)		:: s1
        type(t_adaptive_statistics)						:: sr

		sr%adaptive_stats = -s1%adaptive_stats
		sr%counter_stats = -s1%counter_stats
    end function

    pure function t_adaptive_statistics_sub(s1, s2) result(sr)
        class(t_adaptive_statistics), intent(in)		:: s1
        type(t_adaptive_statistics), intent(in)	        :: s2
        type(t_adaptive_statistics)						:: sr

		sr = s1%add(s2%negate())
    end function

    pure function t_adaptive_statistics_sub_stats(s1, s2) result(sr)
        class(t_adaptive_statistics), intent(in)		:: s1
        type(t_statistics), intent(in)		            :: s2
        type(t_adaptive_statistics)						:: sr

		sr = s1%add(s2%negate())
    end function

    pure function t_adaptive_statistics_scale(s, scaling) result(sr)
        class(t_adaptive_statistics), intent(in)		:: s
        double precision, intent(in)                    :: scaling
        type(t_adaptive_statistics)		                :: sr

		sr%adaptive_stats = s%adaptive_stats * scaling
		sr%counter_stats = int(s%counter_stats * scaling, selected_int_kind(16))
    end function

    subroutine t_adaptive_statistics_start_time(s, time_type)
        class(t_adaptive_statistics), intent(inout)			:: s
        integer, intent(in)			                        :: time_type

        assert_pure(time_type > 0)
        assert_pure(time_type < STATS_TYPE_END_ADAPTIVE_COUNTER)

		s%adaptive_stats(time_type) = s%adaptive_stats(time_type) - get_wtime()
    end subroutine

    subroutine t_adaptive_statistics_stop_time(s, time_type)
        class(t_adaptive_statistics), intent(inout)			:: s
        integer, intent(in)			                        :: time_type

        assert_pure(time_type > 0)
        assert_pure(time_type < STATS_TYPE_END_ADAPTIVE_COUNTER)

		s%adaptive_stats(time_type) = s%adaptive_stats(time_type) + get_wtime()
    end subroutine

    elemental function t_adaptive_statistics_get_time(s, time_type) result(time_value)
        class(t_adaptive_statistics), intent(in)    :: s
        integer, intent(in)			                :: time_type
        double precision                            :: time_value

        assert_pure(time_type > 0)
        assert_pure(time_type < STATS_TYPE_END_ADAPTIVE_COUNTER)

		time_value = s%adaptive_stats(time_type)
    end function

    elemental subroutine t_adaptive_statistics_clear_time(s, time_type)
        class(t_adaptive_statistics), intent(inout)			:: s
        integer, intent(in)			                        :: time_type

        assert_pure(time_type > 0)
        assert_pure(time_type < STATS_TYPE_END_ADAPTIVE_COUNTER)

		s%adaptive_stats(time_type) = 0.d0
    end subroutine

    pure function t_adaptive_statistics_to_string(s) result(str)
        class(t_adaptive_statistics), intent(in)	:: s
		character (len = 512)					    :: str

        write(str, '("#travs: ", I0, " time: ", F0.4, " s (comp: ", F0.4, " s (pre: ", F0.4, " s inner: ", F0.4, " s post: ", F0.4, " s) asagi: ", F0.4, " s sync: ", F0.4, " s barr: ", F0.4, " s update distances: ", F0.4, " s update neighbors: ", F0.4, " s) integrity: ", F0.4, " s load balancing: ", F0.4, " s (de)allocation: ", F0.4, " s")') &
            s%get_counter(traversals), s%get_time(traversal_time), &
            s%get_time(pre_compute_time) + s%get_time(inner_compute_time) + s%get_time(post_compute_time), &
            s%get_time(pre_compute_time), s%get_time(inner_compute_time), s%get_time(post_compute_time), &
            s%get_time(asagi_time), s%get_time(sync_time), s%get_time(barrier_time), &
            s%get_time(update_distances_time), s%get_time(update_neighbors_time), &
            s%get_time(integrity_time), s%get_time(load_balancing_time), s%get_time(allocation_time)

        if (s%get_time(traversal_time) > 0.0d0) then
            write(str, '(A, " ET: ", F0.4, " M/s  MT: ", F0.4, " GB/s")') trim(str), dble(s%get_counter(traversed_cells)) / (1.0d6 * s%get_time(traversal_time)), &
            dble(s%get_counter(traversed_memory)) / (1024.0d0 * 1024.0d0 * 1024.0d0 * s%get_time(traversal_time))
        else
            write(str, '(A, " ET: ", F0.4, " M/s  MT: ", F0.4, " GB/s")') trim(str), -1.0d0, -1.0d0
        end if
	end function
end module
