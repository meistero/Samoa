! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_ASAGI)
#	include 'asagi.f90'
#endif

MODULE SFC_data_types
#	if defined(_TESTS)
		use Tests_data_types
#	elif defined(_HEAT_EQ)
		use Heat_Eq_data_types
#	elif defined(_DARCY)
		use Darcy_data_types
#	elif defined(_SWE)
		use SWE_data_types
#	elif defined(_NUMA)
		use NUMA_data_types
#	elif defined(_PYOP2)
		use pyop2_data_types
#	endif

	use Tools_log
    use Tools_mpi
    use omp_lib

#	if defined(_ASAGI)
		use asagi, asagi_create => grid_create, asagi_open => grid_open, asagi_close => grid_close
#	endif

    implicit none


	!constants

    enum, bind(c)
        enumerator :: K = 1, V = 2, H = 3
    end enum

    enum, bind(c)
        enumerator :: LEFT_NODE = 1, TOP_NODE = 2, RIGHT_NODE = 3
    end enum

    enum, bind(c)
        enumerator :: RIGHT_EDGE = 1, HYPOTENUSE = 2, LEFT_EDGE = 3
    end enum

    enum, bind(c)
        enumerator ::   OLD = 0, NEW = 1, OLD_BND = 2, NEW_BND = 3
        enumerator ::   INNER_OLD       = 16 * OLD + 4 * NEW + OLD
        enumerator ::   INNER_NEW       = 16 * OLD + 4 * NEW + NEW
        enumerator ::   INNER_OLD_BND   = 16 * OLD + 4 * NEW + OLD_BND
        enumerator ::   INNER_NEW_BND   = 16 * OLD + 4 * NEW + NEW_BND
        enumerator ::   FIRST_NEW       = 16 * OLD_BND + 4 * NEW + NEW
        enumerator ::   FIRST_OLD_BND   = 16 * OLD_BND + 4 * NEW + OLD_BND
        enumerator ::   FIRST_NEW_BND   = 16 * OLD_BND + 4 * NEW + NEW_BND
        enumerator ::   LAST_OLD        = 16 * OLD + 4 * NEW_BND + OLD
        enumerator ::   LAST_OLD_BND    = 16 * OLD + 4 * NEW_BND + OLD_BND
        enumerator ::   LAST_NEW_BND    = 16 * OLD + 4 * NEW_BND + NEW_BND
        enumerator ::   SINGLE_OLD_BND  = 16 * OLD_BND + 4 * NEW_BND + OLD_BND
        enumerator ::   SINGLE_NEW_BND  = 16 * OLD_BND + 4 * NEW_BND + NEW_BND
    end enum

    enum, bind(c)
        enumerator ::   RED = -1, GREEN = 0
    end enum

	character (LEN = 5), dimension(RED : GREEN), parameter 				:: color_to_char = [ 'RED', 'GREEN']
	character (LEN = 1), dimension(K : H), parameter 					:: turtle_type_to_char = [ 'K', 'V', 'H' ]
	character (LEN = 11), dimension(OLD : NEW_BND), parameter		    :: edge_type_to_char = [ 'OLD', 'NEW', 'OLD_BND', 'NEW_BND']

    integer (kind = 1), parameter                                 		:: MAX_DEPTH = 8_1 * sizeof(1.0_GRID_SR) - 4_1
	real (kind = GRID_SR), parameter									:: PI = 3.14159265358979323846_GRID_SR 		!< PI. Apparently, "_GRID_SR" is necessary to avoid digit truncation

	integer 				:: rank_MPI = 0
	integer 				:: size_MPI = 1

	real (kind = GRID_SR)   :: r_asagi_time = 0.0_GRID_SR, r_asagi_time_initial = 0.0_GRID_SR

	!********************************
	!Generic scenario data structures
	!********************************

	!---------- this defines a triangle element
	type triangle
		integer (kind = 1)									:: i_triangle			! tree node: -1, tree leaf: -10 (inner cell), -20 (dirichlet boundary), -30 (neumann boundary)
	end type triangle

	type, extends(triangle)									:: fem_triangle
		! SFC nodes: coords(:,1) and coords(:,3) are the nodes of the hypotenuse
		real (kind = GRID_SR), dimension(2,3)				:: r_coords				! 2D coordinates of triangle

		! edge boundary data: boundary_data(:,2) is the hypotenuse
		logical (kind = GRID_SL), dimension(2,3)			:: l_boundary_data
		! 00 - internal, 01 - Dirichlet, 10 - Neumann, 11 - Process Boundary

		! temporary storage
		real (kind = GRID_SR), dimension(2)					:: r_bisection			! = (r_coords(:,1) + r_coords(:3)) / 2	(SFC Nodes)

		integer (kind = 1)							        :: i_plotter_type
		integer (kind = 1)									:: i_turtle_type

		logical (kind = GRID_SL)							:: l_color				! reference color
		integer (kind = 1)									:: i_depth				! grid depth
	end type fem_triangle

	type fine_triangle
		integer (kind = 1)									:: i_edge_types			! encodes the types of previous, color and next edge (OLD, NEW, DOMAIN_BND, PROCESS_BND)	max. 6 bit
		integer (kind = 1)									:: i_depth				! grid depth (0 to MAX_DEPTH)																max. 6? bit
		integer (kind = 1)									:: refinement			! refinement info (-1: coarsen, 0: keep, 1-4: refine once or multiple times)				max. 3 bit
		integer (kind = 1)									:: i_plotter_type		! plotter grammar type for cell orientation (1 to 8)										max. 3 bit
		integer (kind = 1)									:: i_turtle_type		! turtle grammar type for edge/node indexing) (K = 1, V = 2, H = 3)							max. 2 bit
		logical (kind = GRID_SL)							:: l_color_edge_color	! color of the color_edge																	max. 1 bit

		contains

		procedure, pass	:: get_edge_types => cell_get_edge_types
		procedure, pass	:: get_previous_edge_type => cell_get_previous_edge_type
		procedure, pass	:: get_color_edge_type => cell_get_color_edge_type
		procedure, pass	:: get_next_edge_type => cell_get_next_edge_type

		procedure, pass	:: set_edge_types => cell_set_edge_types
		procedure, pass	:: set_previous_edge_type => cell_set_previous_edge_type
		procedure, pass	:: set_color_edge_type => cell_set_color_edge_type
		procedure, pass	:: set_next_edge_type => cell_set_next_edge_type

		procedure, pass	:: get_edge_indices => cell_get_edge_indices
		procedure, pass	:: get_node_indices => cell_get_node_indices

		procedure, pass	:: get_scaling => cell_get_scaling
		procedure, pass	:: get_volume => cell_get_volume
		procedure, pass	:: get_leg_size => cell_get_leg_size
		procedure, pass	:: get_hypo_size => cell_get_hypo_size
		procedure, pass	:: get_edge_sizes => cell_get_edge_sizes

		procedure, pass	:: reverse => cell_reverse
		procedure, pass	:: reverse_inner => cell_reverse_inner
		procedure, pass :: reverse_refinement => cell_reverse_refinement
		procedure, pass :: to_string => cell_to_string
	end type fine_triangle

	!stack & stream data types

	type t_statistics
        real (kind = GRID_SR)					        	:: r_traversal_time = 0.0_GRID_SR
        real (kind = GRID_SR)					        	:: r_computation_time = 0.0_GRID_SR
        real (kind = GRID_SR)					        	:: r_sync_time = 0.0_GRID_SR
        real (kind = GRID_SR)					        	:: r_barrier_time = 0.0_GRID_SR
        integer (kind = GRID_SI)					        :: i_traversals = 0
        integer (kind = GRID_DI)					        :: i_traversed_cells = 0
        integer (kind = GRID_DI)					        :: i_traversed_edges = 0
        integer (kind = GRID_DI)					        :: i_traversed_nodes = 0
        integer (kind = GRID_DI)					        :: i_traversed_memory = 0

		contains

		procedure, private, pass :: reduce_local => t_statistics_reduce_local
		procedure, private, pass :: reduce_global => t_statistics_reduce_global
		procedure, private, pass :: add => t_statistics_add
		procedure, private, pass :: inv => t_statistics_inv
		procedure, private, pass :: sub => t_statistics_sub

        procedure, pass :: to_string => t_statistics_to_string

		generic :: reduce => reduce_local, reduce_global
		generic :: operator(+) => add
		generic :: operator(-) => inv, sub
	end type


	type, extends(t_statistics) ::  t_adaptive_statistics
        real (kind = GRID_SR)					        	:: r_allocation_time = 0.0_GRID_SR
        real (kind = GRID_SR)					        	:: r_update_distances_time = 0.0_GRID_SR
        real (kind = GRID_SR)					        	:: r_update_neighbors_time = 0.0_GRID_SR
        real (kind = GRID_SR)					        	:: r_integrity_time = 0.0_GRID_SR
        real (kind = GRID_SR)					        	:: r_load_balancing_time = 0.0_GRID_SR

		contains

		procedure, private, pass :: reduce_local_adaptive => t_adaptive_statistics_reduce_local
		procedure, private, pass :: reduce_global_adaptive => t_adaptive_statistics_reduce_global
		procedure, private, pass :: add_adaptive => t_adaptive_statistics_add
		procedure, private, pass :: inv_adaptive => t_adaptive_statistics_inv
		procedure, private, pass :: sub_adaptive => t_adaptive_statistics_sub

        procedure, pass :: to_string => t_adaptive_statistics_to_string

		generic :: reduce => reduce_local_adaptive, reduce_global_adaptive
		generic :: operator(+) => add_adaptive
		generic :: operator(-) => inv_adaptive, sub_adaptive
	end type

	type, extends(num_global_data) :: t_global_data
        integer (kind = GRID_SI)					        :: i_sections_per_thread
        integer (kind = 1)					                :: i_min_depth, i_max_depth
        integer (kind = GRID_DI), dimension(RED : GREEN)	:: start_distance, min_distance, end_distance

        integer (kind = GRID_SI), dimension(RED : GREEN)	:: start_dest_stack, end_dest_stack, min_dest_stack, max_dest_stack
        integer (kind = GRID_DI)	                        :: dest_cells, last_dest_cell
        real (kind = GRID_SR)                               :: load, partial_load
        logical	                                            :: l_conform

        contains

        procedure, pass :: to_string => t_global_data_to_string
    end type

	!cell storage

	!> Cell stream data
	type, extends(fine_triangle) 		:: t_cell_stream_data
		type(num_cell_data_pers)		:: data_pers
	end type

	!> Cell local data
	type t_cell_local_data
		type(num_edge_data_temp)		:: data_temp
	end type

	! Edge storage

	!> Edge geometry data structure
	type t_edge_geometry
		logical (kind = GRID_SL)							:: refine				!< if true, the edge will be split for refinement (initially false, if set to true the edge will be refined)
		logical (kind = GRID_SL)							:: coarsen				!< if true, the edge will either be merged or removed for coarsening (initially true, if set to false the edge will not be coarsened)
		logical (kind = GRID_SL)							:: remove				!< if true, the edge will be removed for coarsening (initially true, if remove is set to false the edge will not be removed)

		contains

		procedure, pass :: to_string => edge_to_string
	end type

	!> Crossed edge stream data
	type, extends(t_edge_geometry)							:: t_crossed_edge_stream_data
		type(num_edge_data_pers)							:: data_pers
	end type

	!> Color edge stream data
	type, extends(t_crossed_edge_stream_data)				:: t_color_edge_stream_data
#		if defined(_USE_SKELETON_OP)
			type(num_cell_rep)								:: rep
#		endif
	end type

	!> Color edge stack data
	type, extends(t_color_edge_stream_data)					:: t_edge_data
		type(num_edge_data_temp)							:: data_temp

#		if defined(_USE_SKELETON_OP)
			type(num_cell_update)							:: update
#		endif

        integer (kind = GRID_DI)                            :: min_distance         !< edge minimum distance (only defined for boundary edges)
        logical (kind = GRID_SL)                            :: owned_locally        !< if true, the current section owns the edge, defined only for boundary nodes!
        logical (kind = GRID_SL)                            :: owned_globally       !< if true, the current rank owns the edge, defined only for boundary nodes!
        integer (kind = 1)                                  :: depth                !< edge depth

		type(t_edge_transform_data), pointer				:: transform_data		!< local edge transform data
	end type

	! Node storage

	!> Node geometry data structure
	type t_node_geometry
#		if defined(_STORE_NODE_COORDS)
			real (kind = GRID_SR), dimension(2)				:: position
#		endif
	end type

	!> Node stream data
	type, extends(t_node_geometry)							:: t_node_stream_data
		type(num_node_data_pers)							:: data_pers
	end type

	!> Node stack data
	type, extends(t_node_stream_data)						:: t_node_data
		type(num_node_data_temp)							:: data_temp

#		if .not. defined(_STORE_NODE_COORDS)
			real (kind = GRID_SR), dimension(2)				:: position
#		endif

        integer (kind = GRID_DI)                            :: distance             !< node distance, defined only for boundary nodes!
        logical (kind = GRID_SL)                            :: owned_locally        !< if true, the current section owns the node, defined only for boundary nodes!
        logical (kind = GRID_SL)                            :: owned_globally       !< if true, the current rank owns the node, defined only for boundary nodes!
	end type

	!traversal data types

	!cell geom* + refinement* + pers* + temp* data
	type t_cell_data_ptr
		type(fine_triangle), pointer						:: geometry
		type(num_cell_data_pers), pointer					:: data_pers
		type(num_cell_data_temp), pointer					:: data_temp
	end type

	!> Edge pointer
	type t_edge_ptr
		type(t_edge_data), pointer						    :: ptr					!< pointer to edge data
	end type

	!> Node pointer
	type t_node_ptr
		type(t_node_data), pointer							:: ptr					!< pointer to node data
	end type

	!***********************************************************
	!Generic triangle <-> Reference triangle transformation data
	!***********************************************************

	!> Reference edge transformation data
	type t_edge_transform_data
#		if defined(_USE_SKELETON_OP)
			real (kind = GRID_SR), dimension(2)				:: normal					!< local edge normal in cartesian coordinates
#		endif

		integer (kind = 1)									:: index					!< local edge index
		integer (kind = 1)									:: orientation				!< local edge orientation: -1: backward 1: forward
	end type

	!> Reference cell transformation data for all
	type t_cell_transform_data
		!integer (kind = 1)									:: plotter_type				!< Sierpinski plotter grammar triangle type (appears not to be used anywhere?)
		integer (kind = 1)									:: forward					!< true if forward traversal
		integer (kind = 1)									:: orientation				!< local orientation: -1: backward 1: forward
		real (kind = GRID_SR), DIMENSION(2, 2)				:: jacobian					!< Jacobian of the reference element transformation from barycentric to cartesian coordinates
		real (kind = GRID_SR), DIMENSION(2, 2)				:: jacobian_inv				!< Inverse of the Jacobian
		real (kind = GRID_SR)								:: det_jacobian				!< Determinant of the Jacobian

		type(t_edge_transform_data), DIMENSION(3)			:: edges					!< Reference edge data
	end type

	!> Element-specific data for the generic triangle <-> Reference triangle transformation
	type t_custom_transform_data
		!TODO: real (kind = GRID_SR), DIMENSION(2:2)		:: matrix					!< Element matrix
		!TODO: real (kind = GRID_SR), DIMENSION(2:2)		:: matrix_inv				!< Inverse element matrix
		real (kind = GRID_SR), DIMENSION(:), pointer		:: offset					!< Element offset
		real (kind = GRID_SR)								:: scaling					!< Element scaling
	end type

	!> All data for the generic triangle <-> Reference triangle transformation
	type t_transform_data
		type(t_cell_transform_data), pointer				:: plotter_data				!< Plotter grammar data
		type(t_custom_transform_data)		 				:: custom_data				!< Element-specific custom data
	end type

	type(t_cell_transform_data), DIMENSION(-8 : 8), target	:: ref_plotter_data			!< Reference plotter grammar data for the 16 possible triangle orientations

	interface get_c_pointer
        module procedure t_node_data_get_c_pointer
        module procedure t_edge_data_get_c_pointer
    end interface

	!***********************************************************
	!Element data types
	!***********************************************************

	!> Generic element interface type
	type t_element_base
		type(t_cell_data_ptr)								:: cell												!< cell

		type(t_edge_ptr)									:: previous_edge, next_edge, color_edge				!< edges, by type
		type(t_node_ptr)									:: color_node_in, color_node_out, transfer_node		!< nodes, by type

		type(t_edge_ptr), dimension(3)						:: edges
		type(t_node_ptr), dimension(3)						:: nodes

		type(t_transform_data)								:: transform_data
	end type

	contains

	function t_edge_data_get_c_pointer(array) result(ptr)
        type(t_edge_data), pointer, intent(inout)	:: array(:)
        type(t_edge_data), pointer					:: ptr

        if (.not. associated(array) .or. size(array) .eq. 0) then
            ptr => null()
        else if (loc(array(1)) < loc(array(size(array)))) then
            ptr => array(1)
        else
            ptr => array(size(array))
        end if
	end function

	function t_node_data_get_c_pointer(array) result(ptr)
        type(t_node_data), pointer, intent(inout)	:: array(:)
        type(t_node_data), pointer					:: ptr

        if (.not. associated(array) .or. size(array) .eq. 0) then
            ptr => null()
        else if (loc(array(1)) < loc(array(size(array)))) then
            ptr => array(1)
        else
            ptr => array(size(array))
        end if
	end function

    subroutine t_statistics_reduce_local(s, v)
        class(t_statistics), intent(inout)	:: s
        type(t_statistics), intent(in)		:: v(:)
        real (kind = GRID_SR)               :: scaling

        scaling = 1.0 / real(max(1, size(v)), GRID_SR)

		call reduce(s%r_traversal_time, v%r_traversal_time, MPI_SUM, .false.)
		call reduce(s%r_computation_time, v%r_computation_time, MPI_SUM, .false.)
		call reduce(s%r_sync_time, v%r_sync_time, MPI_SUM, .false.)
		call reduce(s%r_barrier_time, v%r_barrier_time, MPI_SUM, .false.)
		call reduce(s%i_traversals, v%i_traversals, MPI_MAX, .false.)

		call reduce(s%i_traversed_cells, v%i_traversed_cells, MPI_SUM, .false.)
		call reduce(s%i_traversed_edges, v%i_traversed_edges, MPI_SUM, .false.)
		call reduce(s%i_traversed_nodes, v%i_traversed_nodes, MPI_SUM, .false.)
		call reduce(s%i_traversed_memory, v%i_traversed_memory, MPI_SUM, .false.)

		s%r_traversal_time = max(0.0, s%r_traversal_time * scaling)
		s%r_computation_time = max(0.0, s%r_computation_time * scaling)
		s%r_sync_time = max(0.0, s%r_sync_time * scaling)
		s%r_barrier_time = max(0.0, s%r_barrier_time * scaling)
		s%i_traversals = max(1, s%i_traversals)
     end subroutine

    subroutine t_adaptive_statistics_reduce_local(s, v)
        class(t_adaptive_statistics), intent(inout)	    :: s
        type(t_adaptive_statistics), intent(in)		    :: v(:)
        real (kind = GRID_SR)               			:: scaling

        scaling = 1.0 / real(max(1, size(v)), GRID_SR)

		call t_statistics_reduce_local(s, v%t_statistics)

		call reduce(s%r_update_distances_time, v%r_update_distances_time, MPI_SUM, .false.)
		call reduce(s%r_update_neighbors_time, v%r_update_neighbors_time, MPI_SUM, .false.)

		s%r_update_distances_time = max(0.0, s%r_update_distances_time * scaling)
		s%r_update_neighbors_time = max(0.0, s%r_update_neighbors_time * scaling)

		!allocation, integrity and load balancing time is measured globally, not per section
		!so we do not reduce these locally
    end subroutine

    subroutine t_statistics_reduce_global(s)
        class(t_statistics), intent(inout)		        :: s
        real (kind = GRID_SR)               :: scaling

        scaling = 1.0 / real(size_MPI, GRID_SR)

		call reduce(s%r_traversal_time, mpi_op=MPI_SUM)
		call reduce(s%r_computation_time, mpi_op=MPI_SUM)
		call reduce(s%r_sync_time, mpi_op=MPI_SUM)
		call reduce(s%r_barrier_time, mpi_op=MPI_SUM)
		call reduce(s%i_traversals, mpi_op=MPI_MAX)

		call reduce(s%i_traversed_cells, mpi_op=MPI_SUM)
		call reduce(s%i_traversed_edges, mpi_op=MPI_SUM)
		call reduce(s%i_traversed_nodes, mpi_op=MPI_SUM)
		call reduce(s%i_traversed_memory, mpi_op=MPI_SUM)

		s%r_traversal_time = max(0.0, s%r_traversal_time * scaling)
		s%r_computation_time = max(0.0, s%r_computation_time * scaling)
		s%r_sync_time = max(0.0, s%r_sync_time * scaling)
		s%r_barrier_time = max(0.0, s%r_barrier_time * scaling)
		s%i_traversals = max(1, s%i_traversals)
     end subroutine

    subroutine t_adaptive_statistics_reduce_global(s)
        class(t_adaptive_statistics), intent(inout)		        :: s
        real (kind = GRID_SR)               :: scaling

        scaling = 1.0 / real(size_MPI, GRID_SR)

        call t_statistics_reduce_global(s%t_statistics)

		call reduce(s%r_allocation_time, mpi_op=MPI_SUM)
		call reduce(s%r_integrity_time, mpi_op=MPI_SUM)
		call reduce(s%r_load_balancing_time, mpi_op=MPI_SUM)
		call reduce(s%r_update_distances_time, mpi_op=MPI_SUM)
		call reduce(s%r_update_neighbors_time, mpi_op=MPI_SUM)

		s%r_allocation_time = max(0.0, s%r_allocation_time * scaling)
		s%r_integrity_time = max(0.0, s%r_integrity_time * scaling)
		s%r_load_balancing_time = max(0.0, s%r_load_balancing_time * scaling)
		s%r_update_distances_time = max(0.0, s%r_update_distances_time * scaling)
		s%r_update_neighbors_time = max(0.0, s%r_update_neighbors_time * scaling)
     end subroutine

    elemental function t_statistics_add(s1, s2) result(sr)
        class(t_statistics), intent(in)			:: s1, s2
        type(t_statistics)						:: sr

		sr%r_traversal_time = s1%r_traversal_time + s2%r_traversal_time
		sr%r_computation_time = s1%r_computation_time + s2%r_computation_time
		sr%r_sync_time = s1%r_sync_time + s2%r_sync_time
		sr%r_barrier_time = s1%r_barrier_time + s2%r_barrier_time
		sr%i_traversals = s1%i_traversals + s2%i_traversals
		sr%i_traversed_cells = s1%i_traversed_cells + s2%i_traversed_cells
		sr%i_traversed_edges = s1%i_traversed_edges + s2%i_traversed_edges
		sr%i_traversed_nodes = s1%i_traversed_nodes + s2%i_traversed_nodes
		sr%i_traversed_memory = s1%i_traversed_memory + s2%i_traversed_memory
    end function

    elemental function t_adaptive_statistics_add(s1, s2) result(sr)
        class(t_adaptive_statistics), intent(in)		:: s1, s2
        type(t_adaptive_statistics)						:: sr

        sr%t_statistics = t_statistics_add(s1%t_statistics, s2%t_statistics)

		sr%r_allocation_time = s1%r_allocation_time + s2%r_allocation_time
		sr%r_integrity_time = s1%r_integrity_time + s2%r_integrity_time
		sr%r_load_balancing_time = s1%r_load_balancing_time + s2%r_load_balancing_time
		sr%r_update_distances_time = s1%r_update_distances_time + s2%r_update_distances_time
		sr%r_update_neighbors_time = s1%r_update_neighbors_time + s2%r_update_neighbors_time
    end function

    elemental function t_statistics_inv(s1) result(sr)
        class(t_statistics), intent(in)			:: s1
        type(t_statistics)						:: sr

		sr%r_traversal_time = -s1%r_traversal_time
		sr%r_computation_time = -s1%r_computation_time
		sr%r_sync_time = -s1%r_sync_time
		sr%r_barrier_time = -s1%r_barrier_time
		sr%i_traversals = -s1%i_traversals
		sr%i_traversed_cells = -s1%i_traversed_cells
		sr%i_traversed_edges = -s1%i_traversed_edges
		sr%i_traversed_nodes = -s1%i_traversed_nodes
		sr%i_traversed_memory = -s1%i_traversed_memory
    end function

    elemental function t_adaptive_statistics_inv(s1) result(sr)
        class(t_adaptive_statistics), intent(in)		:: s1
        type(t_adaptive_statistics)						:: sr

        sr%t_statistics = -s1%t_statistics

		sr%r_allocation_time = -s1%r_allocation_time
		sr%r_integrity_time = -s1%r_integrity_time
		sr%r_load_balancing_time = -s1%r_load_balancing_time
		sr%r_update_distances_time = -s1%r_update_distances_time
		sr%r_update_neighbors_time = -s1%r_update_neighbors_time
    end function

    elemental function t_statistics_sub(s1, s2) result(sr)
        class(t_statistics), intent(in)			:: s1, s2
        type(t_statistics)						:: sr

		sr = s1 + (-s2)
    end function

    elemental function t_adaptive_statistics_sub(s1, s2) result(sr)
        class(t_adaptive_statistics), intent(in)		:: s1
        type(t_adaptive_statistics), intent(in)		    :: s2
        type(t_adaptive_statistics)						:: sr

		sr = s1 + (-s2)
    end function

    elemental function t_statistics_to_string(s) result(str)
        class(t_statistics), intent(in)			:: s
		character (len = 256)					:: str

        write(str, '("#travs: ", I0, " time: ", F0.4, " s (comp: ", F0.4, " s sync: ", F0.4, " s barr: ", F0.4, " s)")') s%i_traversals, s%r_traversal_time, s%r_computation_time, s%r_sync_time, s%r_barrier_time

        if (s%r_traversal_time > 0.0_GRID_SR) then
            write(str, '(A, " ET: ", F0.4, " M/s  MT: ", F0.4, " GB/s")') trim(str), real(s%i_traversed_cells, GRID_SR) / (1.0e6_GRID_SR * s%r_traversal_time), real(s%i_traversed_memory, GRID_SR) / (1024.0_GRID_SR * 1024.0_GRID_SR * 1024.0_GRID_SR * s%r_traversal_time)
        else
            write(str, '(A, " ET: ", F0.4, " M/s  MT: ", F0.4, " GB/s")') trim(str), -1.0_GRID_SR, -1.0_GRID_SR
        end if
	end function

    elemental function t_adaptive_statistics_to_string(s) result(str)
        class(t_adaptive_statistics), intent(in)	:: s
		character (len = 256)					    :: str

        write(str, '(A, " update distances: ", F0.4, " update neighbors: ", F0.4, " integrity: ", F0.4, "s load balancing: ", F0.4, " s (de)allocation: ", F0.4, " s")'), trim(t_statistics_to_string(s)), s%r_update_distances_time, s%r_update_neighbors_time, s%r_integrity_time, s%r_load_balancing_time, s%r_allocation_time
	end function

    elemental function t_global_data_to_string(gd) result(str)
        class(t_global_data), intent(in)		:: gd
		character (len = 256)					:: str

		write(str, '(A, I0, X, I0, X, I0)') "depth (min, max): ", gd%i_min_depth, gd%i_max_depth
		write(str, '(A, A, F0.4, X, F0.4, X, F0.4)') trim(str), " distance RED (start, min, end): ", decode_distance(gd%start_distance(RED)), decode_distance(gd%min_distance(RED)), decode_distance(gd%end_distance(RED))
		write(str, '(A, A, F0.4, X, F0.4, X, F0.4)') trim(str), " distance GREEN (start, min, end): ", decode_distance(gd%start_distance(GREEN)), decode_distance(gd%min_distance(GREEN)), decode_distance(gd%end_distance(GREEN))
    end function

	!edge geometry type-bound procedures

	elemental function edge_to_string(edge) result(str)
		class(t_edge_geometry), intent(in)		:: edge
		character (len = 64)					:: str

		write(str, '(3(A, L))') "refine:", edge%refine, " coarsen:", edge%coarsen,  " remove:", edge%remove
	end function

	!cell geometry type-bound procedures

	elemental subroutine cell_get_edge_indices(cell, i_previous_edge, i_color_edge, i_next_edge)
		class(fine_triangle), intent(in)		:: cell
		integer (KIND = 1), intent(out)			:: i_previous_edge
		integer (KIND = 1), intent(out)			:: i_color_edge
		integer (KIND = 1), intent(out)			:: i_next_edge

		integer (KIND = 1), dimension(3, 3), parameter :: indices = reshape( \
			[LEFT_EDGE, RIGHT_EDGE, HYPOTENUSE, \
			LEFT_EDGE, HYPOTENUSE, RIGHT_EDGE, \
			HYPOTENUSE, LEFT_EDGE, RIGHT_EDGE], [3, 3])

		i_previous_edge = indices(1, cell%i_turtle_type)
		i_color_edge = indices(2, cell%i_turtle_type)
		i_next_edge = indices(3, cell%i_turtle_type)
	end subroutine

	elemental subroutine cell_get_node_indices(cell, i_color_node_out, i_transfer_node, i_color_node_in)
		class(fine_triangle), intent(in)		:: cell
		integer (KIND = 1), intent(out)			:: i_color_node_out
		integer (KIND = 1), intent(out)			:: i_transfer_node
		integer (KIND = 1), intent(out)			:: i_color_node_in

		integer (KIND = 1), dimension(3, 3), parameter :: indices = reshape( \
			[TOP_NODE, LEFT_NODE, RIGHT_NODE, \
			LEFT_NODE, TOP_NODE, RIGHT_NODE, \
			LEFT_NODE, RIGHT_NODE, TOP_NODE], [3, 3])

		i_color_node_out = indices(1, cell%i_turtle_type)
		i_transfer_node = indices(2, cell%i_turtle_type)
		i_color_node_in = indices(3, cell%i_turtle_type)
	end subroutine

	elemental subroutine cell_get_edge_types(cell, i_previous_edge_type, i_color_edge_type, i_next_edge_type)
		class(fine_triangle), intent(in)		:: cell
		integer (KIND = 1), intent(out)			:: i_previous_edge_type
		integer (KIND = 1), intent(out)			:: i_color_edge_type
		integer (KIND = 1), intent(out)			:: i_next_edge_type

		i_previous_edge_type = 0
		i_color_edge_type = 0
		i_next_edge_type = 0
		call mvbits(cell%i_edge_types, 4, 2, i_previous_edge_type, 0)
		call mvbits(cell%i_edge_types, 2, 2, i_next_edge_type, 0)
		call mvbits(cell%i_edge_types, 0, 2, i_color_edge_type, 0)
	end subroutine

	elemental function cell_get_previous_edge_type(cell) result(i_previous_edge_type)
		class(fine_triangle), intent(in)		:: cell
		integer (KIND = 1)						:: i_previous_edge_type

		i_previous_edge_type = 0
		call mvbits(cell%i_edge_types, 4, 2, i_previous_edge_type, 0)
	end function

	elemental function cell_get_color_edge_type(cell) result(i_color_edge_type)
		class(fine_triangle), intent(in)		:: cell
		integer (KIND = 1)						:: i_color_edge_type

		i_color_edge_type = 0
		call mvbits(cell%i_edge_types, 0, 2, i_color_edge_type, 0)
	end function

	elemental function cell_get_next_edge_type(cell) result(i_next_edge_type)
		class(fine_triangle), intent(in)		:: cell
		integer (KIND = 1)						:: i_next_edge_type

		i_next_edge_type = 0
		call mvbits(cell%i_edge_types, 2, 2, i_next_edge_type, 0)
	end function

	elemental subroutine cell_set_edge_types(cell, i_previous_edge_type, i_color_edge_type, i_next_edge_type)
		class(fine_triangle), intent(inout)		:: cell
		integer (KIND = 1), intent(in)			:: i_previous_edge_type
		integer (KIND = 1), intent(in)			:: i_color_edge_type
		integer (KIND = 1), intent(in)			:: i_next_edge_type

		cell%i_edge_types = 0
		call mvbits(i_previous_edge_type, 0, 2, cell%i_edge_types, 4)
		call mvbits(i_next_edge_type, 0, 2, cell%i_edge_types, 2)
		call mvbits(i_color_edge_type, 0, 2, cell%i_edge_types, 0)
	end subroutine

	elemental subroutine cell_set_previous_edge_type(cell, i_previous_edge_type)
		class(fine_triangle), intent(inout)		:: cell
		integer (KIND = 1), intent(in)			:: i_previous_edge_type

		call mvbits(i_previous_edge_type, 0, 2, cell%i_edge_types, 4)
	end subroutine

	elemental subroutine cell_set_color_edge_type(cell, i_color_edge_type)
		class(fine_triangle), intent(inout)		:: cell
		integer (KIND = 1), intent(in)			:: i_color_edge_type

		call mvbits(i_color_edge_type, 0, 2, cell%i_edge_types, 0)
	end subroutine

	elemental subroutine cell_set_next_edge_type(cell, i_next_edge_type)
		class(fine_triangle), intent(inout)		:: cell
		integer (KIND = 1), intent(in)			:: i_next_edge_type

		call mvbits(i_next_edge_type, 0, 2, cell%i_edge_types, 2)
	end subroutine

	elemental subroutine cell_reverse(cell)
		class(fine_triangle), intent(inout)		:: cell

        !invert edge types: swap previous and next edge and invert old/new bits
        cell%i_edge_types = ishftc(cell%i_edge_types, 2, 6)
        cell%i_edge_types = ishftc(cell%i_edge_types, 2, 4)
        cell%i_edge_types = ieor(cell%i_edge_types, B"010101")

        !invert turtle type: (K, V, H) -> (H, V, K)
        cell%i_turtle_type = 4 - cell%i_turtle_type

        !invert plotter type
        cell%i_plotter_type = -cell%i_plotter_type
	end subroutine

	elemental subroutine cell_reverse_inner(cell)
		class(fine_triangle), intent(inout)		:: cell

        !invert edge types: invert old/new color edge bit
        cell%i_edge_types = ieor(cell%i_edge_types, 1)

        !invert turtle type: (K, V, H) -> (H, V, K)
        cell%i_turtle_type = 4 - cell%i_turtle_type

        !invert plotter type
        cell%i_plotter_type = -cell%i_plotter_type
	end subroutine

	elemental subroutine cell_reverse_refinement(cell)
		class(fine_triangle), intent(inout)		:: cell
		integer (kind = 1), parameter           :: flip_2_3(-1:4) = [-1, 0, 1, 3, 2, 4]

        !invert refinement flag: flip 2 and 3
        cell%refinement = flip_2_3(cell%refinement)
	end subroutine

	!> returns the volume of an element
	elemental function cell_get_scaling(cell) result(scaling)
		class(fine_triangle), intent(in)		:: cell
        integer (kind = 1)                      :: i
		real (kind = GRID_SR)					:: scaling
        real (kind = GRID_SR), parameter, dimension(MAX_DEPTH)		:: r_scalings = [ (0.5_GRID_SR ** i, 0.5_GRID_SR ** i, i = 1, MAX_DEPTH/2) ]

		scaling = r_scalings(cell%i_depth)
	end function

	!> returns the volume of an element
	elemental function cell_get_volume(cell) result(volume)
		class(fine_triangle), intent(in)		:: cell
        integer (kind = 1)                      :: i
		real (kind = GRID_SR)					:: volume
        real (kind = GRID_SR), parameter, dimension(MAX_DEPTH)		:: r_volumes = [ (0.5_GRID_SR ** (i + 1), i = 1, MAX_DEPTH) ]

		volume = r_volumes(cell%i_depth)
	end function

	!> returns the leg size of an element
	elemental function cell_get_leg_size(cell) result(edge_size)
		class(fine_triangle), intent(in)		:: cell
		real (kind = GRID_SR)					:: edge_size

        integer (kind = 1)                      :: i
        real (kind = GRID_SR), parameter, dimension(0 : MAX_DEPTH)	:: r_leg_sizes = [ (sqrt(0.5_GRID_SR) ** i, i = 0, MAX_DEPTH) ]

		edge_size = r_leg_sizes(cell%i_depth)
	end function

	!> returns the leg size of an element
	elemental function cell_get_hypo_size(cell) result(edge_size)
		class(fine_triangle), intent(in)		:: cell
		real (kind = GRID_SR)					:: edge_size

        integer (kind = 1)                      :: i
        real (kind = GRID_SR), parameter, dimension(0 : MAX_DEPTH)	:: r_leg_sizes = [ (sqrt(0.5_GRID_SR) ** i, i = 0, MAX_DEPTH) ]

		edge_size = r_leg_sizes(cell%i_depth - 1)
	end function

	!> returns the edge size of an element
	pure function cell_get_edge_sizes(cell) result(edge_sizes)
		class(fine_triangle), intent(in)		:: cell
		real (kind = GRID_SR)           		:: edge_sizes(3)

        integer (kind = 1)                      :: i
        real (kind = GRID_SR), parameter, dimension(0 : MAX_DEPTH)	:: r_leg_sizes = [ (sqrt(0.5_GRID_SR) ** i, i = 0, MAX_DEPTH) ]

		edge_sizes = [r_leg_sizes(cell%i_depth), r_leg_sizes(cell%i_depth - 1), r_leg_sizes(cell%i_depth)]
	end function

	elemental function get_cell_volume(depth) result(volume)
        integer (kind = 1), intent(in)          :: depth
        integer (kind = 1)                      :: i
		real (kind = GRID_SR)					:: volume

        real (kind = GRID_SR), parameter, dimension(MAX_DEPTH)		:: r_volumes = [ (0.5_GRID_SR ** (i + 1), i = 1, MAX_DEPTH) ]

		volume = r_volumes(depth)
	end function

	elemental function cell_to_string(cell) result(str)
		class(fine_triangle), intent(in)		:: cell
		character (len = 64)					:: str

		write(str, '(A, 1X, A, 1X, A, 1X, A, 1X, A, 1X, I0)') turtle_type_to_char(cell%i_turtle_type), color_to_char(cell%l_color_edge_color), edge_type_to_char(cell%get_previous_edge_type()), edge_type_to_char(cell%get_color_edge_type()), edge_type_to_char(cell%get_next_edge_type()), cell%i_depth
	end function

    elemental function get_edge_size(depth) result(edge_size)
        integer (kind = 1), intent(in)          :: depth
		real (kind = GRID_SR)					:: edge_size

        integer (kind = 1)                      :: i
        real (kind = GRID_SR), parameter, dimension(0 : MAX_DEPTH)  :: r_leg_sizes = [ (sqrt(0.5_GRID_SR) ** i, i = 0, MAX_DEPTH) ]

		edge_size = r_leg_sizes(depth)
    end function

    pure function get_edge_sizes(depth) result(edge_sizes)
        integer (kind = 1), intent(in)          :: depth
		real (kind = GRID_SR)					:: edge_sizes(3)

		edge_sizes = get_edge_size([depth, depth - 1_1, depth])
    end function

    elemental function encode_edge_size(depth) result(code)
        integer (kind = 1), intent(in)          :: depth
		integer (kind = GRID_DI)				:: code

        integer (kind = 1)                      :: i
        integer (kind = GRID_DI), parameter, dimension(0 : MAX_DEPTH + 1_1)  :: codes = [(ishft(1_GRID_DI, MAX_DEPTH / 2_1 - i), ishft(1_GRID_DI, MAX_DEPTH / 2_1 - i), i = 0_1, MAX_DEPTH / 2_1)]

        code = codes(depth)
    end function

    elemental function decode_distance(code) result(distance)
        integer (kind = GRID_DI), intent(in)    :: code
        real (kind = GRID_SR)                   :: distance

        real (kind = GRID_SR), parameter		:: scaling = 1.0_GRID_SR / (dble(ishft(1_GRID_DI, MAX_DEPTH / 2_GRID_DI))) !< square root of 2

        distance = real(code, GRID_SR) * scaling
    end function

    !> initializes the mpi communicator
    subroutine init_mpi()
        integer                         :: i_error, mpi_prov_thread_support
        integer(kind=MPI_ADDRESS_KIND)  :: mpi_tag_upper_bound
        logical                         :: mpi_flag, mpi_is_initialized

#       if defined(_MPI)
            call mpi_initialized(mpi_is_initialized, i_error); assert_eq(i_error, 0)

            if (.not. mpi_is_initialized) then
#               if defined(_OMP)
                    call mpi_init_thread(MPI_THREAD_MULTIPLE, mpi_prov_thread_support, i_error); assert_eq(i_error, 0)
                    assert_ge(mpi_prov_thread_support, MPI_THREAD_MULTIPLE)
#               else
                    call mpi_init(i_error); assert_eq(i_error, 0)
#               endif
            else
#               if defined(_OMP)
                    call mpi_query_thread(mpi_prov_thread_support, i_error); assert_eq(i_error, 0)
                    assert_ge(mpi_prov_thread_support, MPI_THREAD_MULTIPLE)
#               endif
            end if

            call mpi_comm_size(MPI_COMM_WORLD, size_MPI, i_error); assert_eq(i_error, 0)
            call mpi_comm_rank(MPI_COMM_WORLD, rank_MPI, i_error); assert_eq(i_error, 0)

            call mpi_comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, mpi_tag_upper_bound, mpi_flag, i_error); assert_eq(i_error, 0)
            assert(mpi_flag)
            assert_ge(mpi_tag_upper_bound, ishft(1, 30) - 1)
#       else
            size_MPI = 1
            rank_MPI = 0
#       endif
    end subroutine

    !> finalizes the mpi communicator
    subroutine finalize_mpi()
        integer                         :: i_error
        logical                         :: mpi_is_finalized

        call mpi_finalized(mpi_is_finalized, i_error); assert_eq(i_error, 0)

        if (.not. mpi_is_finalized) then
#	        if defined(_MPI)
                call mpi_barrier(MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
                call mpi_finalize(i_error); assert_eq(i_error, 0)
#	        endif
        end if
    end subroutine

	!**********************************
	!Transformation data initialization
	!**********************************

	subroutine init_transform_data()
		type(t_cell_transform_data), pointer				:: p_cell_data
		type(t_edge_transform_data), pointer				:: p_edge_data
		integer (kind = 1)									:: i_plotter_type, i
		integer                                             :: i_error
		real (kind = GRID_SR)								:: r_angle
		real (kind = GRID_SR), dimension(2, 3)				:: edge_vectors, edge_normals
		type(t_global_data)                                 :: global_data

		!set transformation matrices for the 8 different plotter grammar patterns

        do i_plotter_type = -8, 8
            if (i_plotter_type .ne. 0) then
                p_cell_data => ref_plotter_data(i_plotter_type)
                !p_cell_data%plotter_type = i_plotter_type

                !set rotation angle
                r_angle = PI / 4.0_GRID_SR * dble(abs(i_plotter_type))

                if (i_plotter_type > 0) then
                    p_cell_data%orientation = 1
                    p_cell_data%jacobian = reshape([cos(r_angle), sin(r_angle), -sin(r_angle), cos(r_angle)], [2, 2])
                else if (i_plotter_type < 0) then
                    p_cell_data%orientation = -1
                    p_cell_data%jacobian = reshape([-sin(r_angle), cos(r_angle), cos(r_angle), sin(r_angle)], [2, 2])
                end if

                !round values to nearest whole numbers
                !(the matrices should contain only whole numbers, these can be represented as reals without numerical error)
                if (iand(abs(i_plotter_type), 1) == 0) then
                    p_cell_data%jacobian = anint(p_cell_data%jacobian)
                else
                    p_cell_data%jacobian = anint(sqrt(2.0_GRID_SR) * p_cell_data%jacobian)
                end if

                p_cell_data%det_jacobian = p_cell_data%jacobian(1, 1) * p_cell_data%jacobian(2, 2) - p_cell_data%jacobian(1, 2) * p_cell_data%jacobian(2, 1)
                p_cell_data%jacobian_inv = 1.0_GRID_SR / p_cell_data%det_jacobian * reshape([ p_cell_data%jacobian(2, 2), -p_cell_data%jacobian(2, 1), -p_cell_data%jacobian(1, 2), p_cell_data%jacobian(1, 1) ], [ 2, 2 ])

                _log_write(7, '(X, A)') "jacobian: "
                _log_write(7, '(2X, 2(F0.4, X), /, 2X, 2(F0.4, X))') p_cell_data%jacobian
                _log_write(7, '(X, A)') "inverse: "
                _log_write(7, '(2X, 2(F0.4, X), /, 2X, 2(F0.4, X))') p_cell_data%jacobian_inv
                _log_write(7, '(X, A)') "determinant: "
                _log_write(7, '(2X, F0.4)') p_cell_data%det_jacobian

                edge_vectors(:, 1) = matmul(p_cell_data%jacobian, [0.0_GRID_SR, 1.0_GRID_SR])
                edge_vectors(:, 2) = matmul(p_cell_data%jacobian, [1.0_GRID_SR, -1.0_GRID_SR])
                edge_vectors(:, 3) = matmul(p_cell_data%jacobian, [-1.0_GRID_SR, 0.0_GRID_SR])

                edge_normals(:, 1) = matmul([-1.0_GRID_SR, 0.0_GRID_SR], p_cell_data%jacobian_inv)
                edge_normals(:, 2) = matmul([1.0_GRID_SR, 1.0_GRID_SR], p_cell_data%jacobian_inv)
                edge_normals(:, 3) = matmul([0.0_GRID_SR, -1.0_GRID_SR], p_cell_data%jacobian_inv)

                do i = 1, 3
                    p_edge_data => p_cell_data%edges(i)
                    p_edge_data%index = i

                    if (edge_vectors(2, i) == 0.0_GRID_SR) then
                        p_edge_data%orientation = sign(1.0_GRID_SR, dot_product(edge_vectors(:, i), [1.0_GRID_SR, 0.0_GRID_SR]))
                    else if (edge_vectors(1, i) == 0.0_GRID_SR) then
                        p_edge_data%orientation = sign(1.0_GRID_SR, dot_product(edge_vectors(:, i), [0.0_GRID_SR, 1.0_GRID_SR]))
                    else if (edge_vectors(1, i) - edge_vectors(2, i) == 0.0_GRID_SR) then
                        p_edge_data%orientation = sign(1.0_GRID_SR, dot_product(edge_vectors(:, i), [1.0_GRID_SR, 1.0_GRID_SR]))
                    else if (edge_vectors(1, i) + edge_vectors(2, i) == 0.0_GRID_SR) then
                        p_edge_data%orientation = sign(1.0_GRID_SR, dot_product(edge_vectors(:, i), [1.0_GRID_SR, -1.0_GRID_SR]))
                    else
                        !something's wrong
                        assert(.false.)
                    end if

                    edge_vectors(:, i) = dble(p_edge_data%orientation) * edge_vectors(:, i)
                    edge_normals(:, i) = 1.0_GRID_SR / sqrt(dot_product(edge_normals(:, i), edge_normals(:, i))) * edge_normals(:, i)

#	    			if defined(_USE_SKELETON_OP)
                        p_edge_data%normal = edge_normals(:, i)
#	    			endif

                    _log_write(7, '(X, A, I0)') "edge: ", p_edge_data%index
                    _log_write(7, '(2X, A, I0)') "orientation: ", p_edge_data%orientation
                    _log_write(7, '(2X, A, 2(F0.4, 1X))') "vector: ", edge_vectors(:, i)
                    _log_write(7, '(2X, A, 2(F0.4, 1X))') "normal: ", edge_normals(:, i)
                end do
            end if
        end do
	end subroutine
end MODULE
