#include "Compilation_control.f90"

module Tools_statistics
    use Config
    use Tools_mpi
    use Tools_log
    use Tools_parallel_operators

    private
    public t_section_statistics, t_statistics, t_adaptive_statistics

	type t_section_statistics
        double precision					        	    :: r_computation_time = 0.0d0
        double precision                                    :: r_asagi_time = 0.0d0
        integer (kind = selected_int_kind(8))			    :: i_traversals = 0
        integer (kind = selected_int_kind(16))			    :: i_traversed_cells = 0
        integer (kind = selected_int_kind(16))			    :: i_traversed_edges = 0
        integer (kind = selected_int_kind(16))			    :: i_traversed_nodes = 0
        integer (kind = selected_int_kind(16))			    :: i_traversed_memory = 0
	end type

	type, extends(t_section_statistics) :: t_statistics
        double precision					        	    :: r_traversal_time = 0.0d0
        double precision					        	    :: r_pre_compute_time = 0.0d0
        double precision					        	    :: r_post_compute_time = 0.0d0
        double precision					        	    :: r_sync_time = 0.0d0
        double precision					        	    :: r_barrier_time = 0.0d0

		contains

		procedure, private, pass :: reduce_local => t_statistics_reduce_local
		procedure, private, pass :: reduce_global => t_statistics_reduce_global
		procedure, private, pass :: add => t_statistics_add
		procedure, private, pass :: inv => t_statistics_inv
		procedure, private, pass :: sub => t_statistics_sub
		procedure, private, pass :: assign_stats => t_statistics_assign

        procedure, pass :: scale => t_statistics_scale
        procedure, pass :: to_string => t_statistics_to_string
		procedure, pass :: clear => t_statistics_clear

		generic :: assignment(=) => assign_stats
		generic :: reduce => reduce_local, reduce_global
		generic :: operator(+) => add
		generic :: operator(-) => sub
	end type

	type, extends(t_statistics) ::  t_adaptive_statistics
        double precision					        	:: r_allocation_time = 0.0d0
        double precision					        	:: r_update_distances_time = 0.0d0
        double precision					        	:: r_update_neighbors_time = 0.0d0
        double precision					        	:: r_integrity_time = 0.0d0
        double precision					        	:: r_load_balancing_time = 0.0d0

		contains

		procedure, private, pass :: reduce_local_adaptive => t_adaptive_statistics_reduce_local
		procedure, private, pass :: reduce_global => t_adaptive_statistics_reduce_global
		procedure, private, pass :: add_adaptive => t_adaptive_statistics_add
		procedure, private, pass :: inv_adaptive => t_adaptive_statistics_inv
		procedure, private, pass :: sub_adaptive => t_adaptive_statistics_sub
		procedure, private, pass :: assign_adaptive => t_adaptive_statistics_assign

        procedure, pass :: scale_adaptive => t_adaptive_statistics_scale
        procedure, pass :: to_string => t_adaptive_statistics_to_string
		procedure, pass :: clear => t_adaptive_statistics_clear

		generic :: assignment(=) => assign_adaptive
		generic :: reduce => reduce_local_adaptive
		generic :: operator(+) => add_adaptive
		generic :: operator(-) => sub_adaptive
	end type

	contains

    subroutine t_statistics_reduce_local(s, v, mpi_op)
        class(t_statistics), intent(inout)	:: s
        type(t_statistics), intent(in)		:: v(:)
        integer, intent(in)                 :: mpi_op

		call reduce(s%r_traversal_time, v%r_traversal_time, mpi_op, .false.)
        call reduce(s%r_pre_compute_time, v%r_pre_compute_time, mpi_op, .false.)
        call reduce(s%r_post_compute_time, v%r_post_compute_time, mpi_op, .false.)
		call reduce(s%r_computation_time, v%r_computation_time, mpi_op, .false.)
		call reduce(s%r_sync_time, v%r_sync_time, mpi_op, .false.)
		call reduce(s%r_barrier_time, v%r_barrier_time, mpi_op, .false.)
		call reduce(s%r_asagi_time, v%r_asagi_time, mpi_op, .false.)

		call reduce(s%i_traversals, v%i_traversals, mpi_op, .false.)
		call reduce(s%i_traversed_cells, v%i_traversed_cells, mpi_op, .false.)
		call reduce(s%i_traversed_edges, v%i_traversed_edges, mpi_op, .false.)
		call reduce(s%i_traversed_nodes, v%i_traversed_nodes, mpi_op, .false.)
		call reduce(s%i_traversed_memory, v%i_traversed_memory, mpi_op, .false.)
     end subroutine

    subroutine t_adaptive_statistics_reduce_local(s, v, mpi_op)
        class(t_adaptive_statistics), intent(inout)	    :: s
        type(t_adaptive_statistics), intent(in)		    :: v(:)
        integer, intent(in)                             :: mpi_op

		call t_statistics_reduce_local(s, v%t_statistics, mpi_op)

		call reduce(s%r_allocation_time, v%r_allocation_time, mpi_op, .false.)
		call reduce(s%r_integrity_time, v%r_integrity_time, mpi_op, .false.)
		call reduce(s%r_load_balancing_time, v%r_load_balancing_time, mpi_op, .false.)
		call reduce(s%r_update_distances_time, v%r_update_distances_time, mpi_op, .false.)
		call reduce(s%r_update_neighbors_time, v%r_update_neighbors_time, mpi_op, .false.)
    end subroutine

    subroutine t_statistics_reduce_global(s, mpi_op)
        class(t_statistics), intent(inout)          :: s
        integer, intent(in)                         :: mpi_op

		call reduce(s%r_traversal_time, mpi_op)
		call reduce(s%r_computation_time, mpi_op)
        call reduce(s%r_pre_compute_time, mpi_op)
        call reduce(s%r_post_compute_time, mpi_op)
		call reduce(s%r_sync_time, mpi_op)
		call reduce(s%r_barrier_time, mpi_op)
		call reduce(s%r_asagi_time, mpi_op)

		call reduce(s%i_traversals, mpi_op)
		call reduce(s%i_traversed_cells, mpi_op)
		call reduce(s%i_traversed_edges, mpi_op)
		call reduce(s%i_traversed_nodes, mpi_op)
		call reduce(s%i_traversed_memory, mpi_op)
     end subroutine

    subroutine t_adaptive_statistics_reduce_global(s, mpi_op)
        class(t_adaptive_statistics), intent(inout)     :: s
        integer, intent(in)                             :: mpi_op

        call t_statistics_reduce_global(s%t_statistics, mpi_op)

		call reduce(s%r_allocation_time, mpi_op)
		call reduce(s%r_integrity_time, mpi_op)
		call reduce(s%r_load_balancing_time, mpi_op)
		call reduce(s%r_update_distances_time, mpi_op)
		call reduce(s%r_update_neighbors_time, mpi_op)
    end subroutine

    subroutine t_adaptive_statistics_clear(s)
        class(t_adaptive_statistics), intent(inout)	    :: s

		s = t_adaptive_statistics()
    end subroutine


    subroutine t_statistics_clear(s)
        class(t_statistics), intent(inout)	    :: s

		s = t_statistics()
    end subroutine

    pure subroutine t_statistics_assign(sr, s2)
        class(t_statistics), intent(inout)	    :: sr
        type(t_statistics), intent(in)			:: s2

		sr%r_traversal_time = s2%r_traversal_time
		sr%r_computation_time = s2%r_computation_time
        sr%r_pre_compute_time = s2%r_pre_compute_time
        sr%r_post_compute_time = s2%r_post_compute_time
		sr%r_sync_time = s2%r_sync_time
		sr%r_barrier_time = s2%r_barrier_time
		sr%r_asagi_time = s2%r_asagi_time

		sr%i_traversals = s2%i_traversals
		sr%i_traversed_cells = s2%i_traversed_cells
		sr%i_traversed_edges = s2%i_traversed_edges
		sr%i_traversed_nodes = s2%i_traversed_nodes
		sr%i_traversed_memory = s2%i_traversed_memory
    end subroutine

    pure function t_statistics_add(s1, s2) result(sr)
        class(t_statistics), intent(in)			:: s1
        type(t_statistics), intent(in)			:: s2
        type(t_statistics)						:: sr

		sr%r_traversal_time = s1%r_traversal_time + s2%r_traversal_time
		sr%r_computation_time = s1%r_computation_time + s2%r_computation_time
        sr%r_pre_compute_time = s1%r_pre_compute_time + s2%r_pre_compute_time
        sr%r_post_compute_time = s1%r_post_compute_time + s2%r_post_compute_time
		sr%r_sync_time = s1%r_sync_time + s2%r_sync_time
		sr%r_barrier_time = s1%r_barrier_time + s2%r_barrier_time
		sr%r_asagi_time = s1%r_asagi_time + s2%r_asagi_time

		sr%i_traversals = s1%i_traversals + s2%i_traversals
		sr%i_traversed_cells = s1%i_traversed_cells + s2%i_traversed_cells
		sr%i_traversed_edges = s1%i_traversed_edges + s2%i_traversed_edges
		sr%i_traversed_nodes = s1%i_traversed_nodes + s2%i_traversed_nodes
		sr%i_traversed_memory = s1%i_traversed_memory + s2%i_traversed_memory
    end function

    pure subroutine t_adaptive_statistics_assign(sr, s2)
        class(t_adaptive_statistics), intent(inout)	    :: sr
        type(t_adaptive_statistics), intent(in)			:: s2

        sr%t_statistics = s2%t_statistics

		sr%r_allocation_time = s2%r_allocation_time
		sr%r_integrity_time = s2%r_integrity_time
		sr%r_load_balancing_time = s2%r_load_balancing_time
		sr%r_update_distances_time = s2%r_update_distances_time
		sr%r_update_neighbors_time = s2%r_update_neighbors_time
    end subroutine

    pure function t_adaptive_statistics_add(s1, s2) result(sr)
        class(t_adaptive_statistics), intent(in)		:: s1
        type(t_adaptive_statistics), intent(in)		    :: s2
        type(t_adaptive_statistics)						:: sr

        sr%t_statistics = t_statistics_add(s1%t_statistics, s2%t_statistics)

		sr%r_allocation_time = s1%r_allocation_time + s2%r_allocation_time
		sr%r_integrity_time = s1%r_integrity_time + s2%r_integrity_time
		sr%r_load_balancing_time = s1%r_load_balancing_time + s2%r_load_balancing_time
		sr%r_update_distances_time = s1%r_update_distances_time + s2%r_update_distances_time
		sr%r_update_neighbors_time = s1%r_update_neighbors_time + s2%r_update_neighbors_time
    end function

    pure function t_statistics_inv(s1) result(sr)
        class(t_statistics), intent(in)			:: s1
        type(t_statistics)						:: sr

		sr%r_traversal_time = -s1%r_traversal_time
		sr%r_computation_time = -s1%r_computation_time
        sr%r_pre_compute_time = -s1%r_pre_compute_time
        sr%r_post_compute_time = -s1%r_post_compute_time
		sr%r_sync_time = -s1%r_sync_time
		sr%r_barrier_time = -s1%r_barrier_time
		sr%r_asagi_time = -s1%r_asagi_time

		sr%i_traversals = -s1%i_traversals
		sr%i_traversed_cells = -s1%i_traversed_cells
		sr%i_traversed_edges = -s1%i_traversed_edges
		sr%i_traversed_nodes = -s1%i_traversed_nodes
		sr%i_traversed_memory = -s1%i_traversed_memory
    end function

    pure function t_adaptive_statistics_inv(s1) result(sr)
        class(t_adaptive_statistics), intent(in)		:: s1
        type(t_adaptive_statistics)						:: sr

        sr%t_statistics = s1%t_statistics%inv()

		sr%r_allocation_time = -s1%r_allocation_time
		sr%r_integrity_time = -s1%r_integrity_time
		sr%r_load_balancing_time = -s1%r_load_balancing_time
		sr%r_update_distances_time = -s1%r_update_distances_time
		sr%r_update_neighbors_time = -s1%r_update_neighbors_time
    end function

    pure function t_statistics_scale(s, scaling) result(sr)
        class(t_statistics), intent(in)		            :: s
        double precision, intent(in)                    :: scaling
        type(t_statistics)		                        :: sr

		sr%r_traversal_time = s%r_traversal_time * scaling
		sr%r_computation_time = s%r_computation_time * scaling
        sr%r_pre_compute_time = s%r_pre_compute_time * scaling
        sr%r_post_compute_time = s%r_post_compute_time * scaling
		sr%r_sync_time = s%r_sync_time * scaling
		sr%r_barrier_time = s%r_barrier_time * scaling
		sr%r_asagi_time = s%r_asagi_time * scaling

		sr%i_traversals = int(s%i_traversals * scaling, selected_int_kind(8))
		sr%i_traversed_cells = int(s%i_traversed_cells * scaling, selected_int_kind(16))
		sr%i_traversed_edges = int(s%i_traversed_edges * scaling, selected_int_kind(16))
		sr%i_traversed_nodes = int(s%i_traversed_nodes * scaling, selected_int_kind(16))
		sr%i_traversed_memory = int(s%i_traversed_memory * scaling, selected_int_kind(16))
    end function

    pure function t_adaptive_statistics_scale(s, scaling) result(sr)
        class(t_adaptive_statistics), intent(in)		:: s
        double precision, intent(in)                    :: scaling
        type(t_adaptive_statistics)		                :: sr

        sr%t_statistics = t_statistics_scale(s%t_statistics, scaling)

		sr%r_allocation_time = max(0.0d0, s%r_allocation_time * scaling)
		sr%r_integrity_time = max(0.0d0, s%r_integrity_time * scaling)
		sr%r_load_balancing_time = max(0.0d0, s%r_load_balancing_time * scaling)
		sr%r_update_distances_time = max(0.0d0, s%r_update_distances_time * scaling)
		sr%r_update_neighbors_time = max(0.0d0, s%r_update_neighbors_time * scaling)
    end function

    pure function t_statistics_sub(s1, s2) result(sr)
        class(t_statistics), intent(in)			:: s1
        type(t_statistics), intent(in)			:: s2
        type(t_statistics)						:: sr

		sr = t_statistics_add(s1, s2%inv())
    end function

    pure function t_adaptive_statistics_sub(s1, s2) result(sr)
        class(t_adaptive_statistics), intent(in)		:: s1
        type(t_adaptive_statistics), intent(in)	        :: s2
        type(t_adaptive_statistics)						:: sr

		sr = t_adaptive_statistics_add(s1, s2%inv_adaptive())
    end function

    pure function t_statistics_to_string(s) result(str)
        class(t_statistics), intent(in)			:: s
		character (len = 512)					:: str

        write(str, '("#travs: ", I0, " time: ", F0.4, " s (comp: ", F0.4, " s (pre: ", F0.4, " s inner: ", F0.4, " s post: ", F0.4, " s) asagi: ", F0.4, " s sync: ", F0.4, " s barr: ", F0.4, " s)")') &
            s%i_traversals, s%r_traversal_time, s%r_computation_time, s%r_pre_compute_time, s%r_computation_time - s%r_pre_compute_time - s%r_post_compute_time, s%r_post_compute_time, s%r_asagi_time, s%r_sync_time, s%r_barrier_time

        if (s%r_traversal_time > 0.0d0) then
            write(str, '(A, " ET: ", F0.4, " M/s  MT: ", F0.4, " GB/s")') trim(str), dble(s%i_traversed_cells) / (1.0d6 * s%r_traversal_time), dble(s%i_traversed_memory) / (1024.0d0 * 1024.0d0 * 1024.0d0 * s%r_traversal_time)
        else
            write(str, '(A, " ET: ", F0.4, " M/s  MT: ", F0.4, " GB/s")') trim(str), -1.0d0, -1.0d0
        end if
	end function

    pure function t_adaptive_statistics_to_string(s) result(str)
        class(t_adaptive_statistics), intent(in)	:: s
		character (len = 512)					    :: str

        write(str, '("#travs: ", I0, " time: ", F0.4, " s (comp: ", F0.4, " s (pre: ", F0.4, " s inner: ", F0.4, " s post: ", F0.4, " s) asagi: ", F0.4, " s sync: ", F0.4, " s barr: ", F0.4, " s update distances: ", F0.4, " s update neighbors: ", F0.4, " s) integrity: ", F0.4, " s load balancing: ", F0.4, " s (de)allocation: ", F0.4, " s")') &
            s%i_traversals, s%r_traversal_time, s%r_computation_time, s%r_pre_compute_time, s%r_computation_time - s%r_pre_compute_time - s%r_post_compute_time, s%r_post_compute_time, s%r_asagi_time, s%r_sync_time, s%r_barrier_time, s%r_update_distances_time, s%r_update_neighbors_time, &
            s%r_integrity_time, s%r_load_balancing_time, s%r_allocation_time

        if (s%r_traversal_time > 0.0d0) then
            write(str, '(A, " ET: ", F0.4, " M/s  MT: ", F0.4, " GB/s")') trim(str), dble(s%i_traversed_cells) / (1.0d6 * s%r_traversal_time), dble(s%i_traversed_memory) / (1024.0d0 * 1024.0d0 * 1024.0d0 * s%r_traversal_time)
        else
            write(str, '(A, " ET: ", F0.4, " M/s  MT: ", F0.4, " GB/s")') trim(str), -1.0d0, -1.0d0
        end if
	end function
end module
