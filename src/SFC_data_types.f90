! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE

#include "Compilation_control.f90"

MODULE SFC_data_types
#	if defined(_TESTS)
		use Tests_data_types
#	elif defined(_HEAT_EQ)
		use Heat_Eq_data_types
#	elif defined(_DARCY)
		use Darcy_data_types
#	elif defined(_SWE)
		use SWE_data_types
#	elif defined(_FLASH)
		use FLASH_data_types
#	elif defined(_NUMA)
		use NUMA_data_types
#	elif defined(_GENERIC)
		use generic_data_types
#	endif

	use Config
	use Tools_log
	use Tools_statistics
    use Tools_mpi
    use Tools_openmp
    use Tools_parallel_operators

#	if defined(_ASAGI)
		use asagi
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
        enumerator ::   INNER_OLD       = 8 * 0 + 4 * 0 + 2 * 0 + OLD
        enumerator ::   INNER_NEW       = 8 * 0 + 4 * 0 + 2 * 0 + NEW
        enumerator ::   INNER_OLD_BND   = 8 * 0 + 4 * 0 + 2 * 1 + OLD
        enumerator ::   INNER_NEW_BND   = 8 * 0 + 4 * 0 + 2 * 1 + NEW
        enumerator ::   FIRST_OLD       = 8 * 1 + 4 * 0 + 2 * 0 + OLD   !this case is invalid and possible only during adaptivity, when not all boundary edges have been found yet
        enumerator ::   FIRST_NEW       = 8 * 1 + 4 * 0 + 2 * 0 + NEW
        enumerator ::   FIRST_OLD_BND   = 8 * 1 + 4 * 0 + 2 * 1 + OLD
        enumerator ::   FIRST_NEW_BND   = 8 * 1 + 4 * 0 + 2 * 1 + NEW
        enumerator ::   LAST_OLD        = 8 * 0 + 4 * 1 + 2 * 0 + OLD
        enumerator ::   LAST_OLD_BND    = 8 * 0 + 4 * 1 + 2 * 1 + OLD
        enumerator ::   LAST_NEW_BND    = 8 * 0 + 4 * 1 + 2 * 1 + NEW
        enumerator ::   SINGLE_OLD_BND  = 8 * 1 + 4 * 1 + 2 * 1 + OLD
        enumerator ::   SINGLE_NEW_BND  = 8 * 1 + 4 * 1 + 2 * 1 + NEW
    end enum

    enum, bind(c)
        enumerator ::   RED = -1, GREEN = 0
    end enum


	character (len = 5), dimension(RED : GREEN), parameter 				:: color_to_char = [ '  RED', 'GREEN']
	character (len = 1), dimension(K : H), parameter 					:: turtle_type_to_char = [ 'K', 'V', 'H' ]
	character (len = 7), dimension(OLD : NEW_BND), parameter		    :: edge_type_to_char = [ '    OLD', '    NEW', 'OLD_BND', 'NEW_BND']

    integer, parameter                                 		            :: MAX_DEPTH = bit_size(1_GRID_DI) - 4
	real (kind = GRID_SR), parameter									:: PI = 4.0_GRID_SR * atan(1.0_GRID_SR) !< PI

	!********************************
	!Generic scenario data structures
	!********************************

	type fine_triangle
		integer (kind = BYTE)									:: i_entity_types		! encodes the types of edges and nodes (old/new, inner/boundary)							max. 8 bit
		integer (kind = BYTE)									:: i_depth				! grid depth (0 to MAX_DEPTH)																max. 6? bit
		integer (kind = BYTE)									:: refinement			! refinement info (-1: coarsen, 0: keep, 1-4: refine once or multiple times)				max. 3 bit
		integer (kind = BYTE)									:: i_plotter_type		! plotter grammar type for cell orientation (-8 to 8)										max. 3 bit
		integer (kind = BYTE)									:: i_turtle_type		! turtle grammar type for edge/node indexing) (K = 1, V = 2, H = 3)							max. 2 bit

        integer (kind = BYTE)							        :: i_color_edge_color	! color of the color_edge (in -1:0)																	max. 1 bit

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

	type, extends(num_global_data) :: t_global_data
        integer (kind = GRID_DI), dimension(RED : GREEN)	:: start_distance, min_distance, end_distance

        integer (kind = GRID_SI), dimension(RED : GREEN)	:: start_dest_stack, end_dest_stack, min_dest_stack, max_dest_stack
        integer (kind = GRID_DI)	                        :: dest_cells, last_dest_cell
        integer (kind = GRID_DI)                            :: load, partial_load
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
		logical 							:: refine				!< if true, the edge will be split for refinement (initially false, if set to true the edge will be refined)
		logical 							:: coarsen				!< if true, the edge will either be merged or removed for coarsening (initially true, if set to false the edge will not be coarsened)
		logical 							:: remove				!< if true, the edge will be removed for coarsening (initially true, if remove is set to false the edge will not be removed)

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

        integer (kind = GRID_DI)                            :: min_distance             !< edge minimum distance (only defined for boundary edges)
        logical                                             :: owned_locally = .true.   !< if true, the current section owns the edge, defined only for boundary edges!
        logical                                             :: owned_globally = .true.  !< if true, the current rank owns the edge, defined only for boundary edges! owned_globally alwas implies owned_locally.
        integer (kind = BYTE)                               :: depth                    !< edge depth

		type(t_edge_transform_data), pointer				:: transform_data		    !< local edge transform data
	end type

	! Node storage

	!> Node geometry data structure
	type t_node_geometry
#		if defined(_STORE_NODE_COORDS)
			real (kind = GRID_SR)				            :: position(2)
#		endif
	end type

	!> Node stream data
	type, extends(t_node_geometry)							:: t_node_stream_data
		type(num_node_data_pers)							:: data_pers
	end type

	!> Node stack data
	type, extends(t_node_stream_data)						:: t_node_data
		type(num_node_data_temp)							:: data_temp

#		if !defined(_STORE_NODE_COORDS)
			real (kind = GRID_SR)			                :: position(2)
#		endif

        integer (kind = GRID_DI)                            :: distance             !< node distance, defined only for boundary nodes!
        logical                                             :: owned_locally = .true.        !< if true, the current section owns the node, defined only for boundary nodes!
        logical                                             :: owned_globally = .true.       !< if true, the current rank owns the node, defined only for boundary nodes! owned_globally alwas implies owned_locally.
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

		integer (kind = BYTE)									:: index					!< local edge index
		integer (kind = BYTE)									:: orientation				!< local edge orientation: -1: backward 1: forward
	end type

	!> Reference cell transformation data for all
	type t_cell_transform_data
		!integer (kind = BYTE)									:: plotter_type				!< Sierpinski plotter grammar triangle type (appears not to be used anywhere?)
		integer (kind = BYTE)									:: forward					!< true if forward traversal
		integer (kind = BYTE)									:: orientation				!< local orientation: -1: backward 1: forward
		real (kind = GRID_SR)				                    :: jacobian(2, 2)		    !< Jacobian of the reference element transformation from barycentric to cartesian coordinates
		real (kind = GRID_SR)			                        :: jacobian_inv(2, 2)	    !< Inverse of the Jacobian
		real (kind = GRID_SR)								    :: det_jacobian				!< Determinant of the Jacobian

		type(t_edge_transform_data)                             :: edges(3)				!< Reference edge data
	end type

	!> Element-specific data for the generic triangle <-> Reference triangle transformation
	type t_custom_transform_data
		real (kind = GRID_SR), pointer		                :: offset(:)				!< Element offset
		real (kind = GRID_SR)								:: scaling					!< Element scaling
	end type

	!> All data for the generic triangle <-> Reference triangle transformation
	type t_transform_data
		type(t_cell_transform_data), pointer				:: plotter_data				!< Plotter grammar data
		type(t_custom_transform_data)		 				:: custom_data				!< Element-specific custom data

		contains

        procedure, pass :: get_jacobian_2DH => t_transform_data_get_jacobian_2DH
        procedure, pass :: get_jacobian_inv_2DH => t_transform_data_get_jacobian_inv_2DH
	end type

	type(t_cell_transform_data), target	                    :: ref_plotter_data(-8 : 8)			!< Reference plotter grammar data for the 16 possible triangle orientations

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

	function t_transform_data_get_jacobian_2DH(td) result(jacobian_2DH)
        class(t_transform_data), intent(in) :: td

        real (kind = GRID_SR) :: jacobian_2DH(3, 3)

        jacobian_2DH(1:2, 1:2) = td%custom_data%scaling * td%plotter_data%jacobian
        jacobian_2DH(3, 1:2) = td%custom_data%offset
        jacobian_2DH(1:3, 3) = [0.0_GRID_SR, 0.0_GRID_SR, 1.0_GRID_SR]
    end function

	function t_transform_data_get_jacobian_inv_2DH(td)  result(jacobian_inv_2DH)
        class(t_transform_data), intent(in) :: td

        real (kind = GRID_SR) :: jacobian_inv_2DH(3, 3)

        jacobian_inv_2DH(1:2, 1:2) = td%plotter_data%jacobian_inv / td%custom_data%scaling
        jacobian_inv_2DH(3, 1:2) = matmul(td%plotter_data%jacobian_inv, td%custom_data%offset) / td%custom_data%scaling
        jacobian_inv_2DH(1:3, 3) = [0.0_GRID_SR, 0.0_GRID_SR, 1.0_GRID_SR]
    end function

	function t_edge_data_get_c_pointer(array) result(ptr)
        type(t_edge_data), pointer, intent(inout)	:: array(:)
        type(t_edge_data), pointer					:: ptr
        type(t_edge_data), target	:: dummy

        if (.not. associated(array) .or. size(array) .eq. 0) then
            ptr => dummy
        else if (loc(array(1)) < loc(array(size(array)))) then
            ptr => array(1)
        else
            ptr => array(size(array))
        end if
	end function

	function t_node_data_get_c_pointer(array) result(ptr)
        type(t_node_data), pointer, intent(inout)	:: array(:)
        type(t_node_data), pointer					:: ptr
        type(t_node_data), target	:: dummy

        if (.not. associated(array) .or. size(array) .eq. 0) then
            ptr => dummy
        else if (loc(array(1)) < loc(array(size(array)))) then
            ptr => array(1)
        else
            ptr => array(size(array))
        end if
	end function

    function t_global_data_to_string(gd) result(str)
        class(t_global_data), intent(in)		:: gd
		character (len = 256)					:: str

		write(str, '(A, F0.4, X, F0.4, X, F0.4)') " distance RED (start, min, end): ", decode_distance(gd%start_distance(RED)), decode_distance(gd%min_distance(RED)), decode_distance(gd%end_distance(RED))
		write(str, '(A, F0.4, X, F0.4, X, F0.4)') " distance GREEN (start, min, end): ", decode_distance(gd%start_distance(GREEN)), decode_distance(gd%min_distance(GREEN)), decode_distance(gd%end_distance(GREEN))
    end function

	!edge geometry type-bound procedures

	function edge_to_string(edge) result(str)
		class(t_edge_geometry), intent(in)		:: edge
		character (len = 64)					:: str

		write(str, '(3(A, L))') "refine:", edge%refine, " coarsen:", edge%coarsen,  " remove:", edge%remove
	end function

	!cell geometry type-bound procedures

	elemental subroutine cell_get_edge_indices(cell, i_previous_edge, i_color_edge, i_next_edge)
		class(fine_triangle), intent(in)		:: cell
		integer (kind = BYTE), intent(out)			:: i_previous_edge
		integer (kind = BYTE), intent(out)			:: i_color_edge
		integer (kind = BYTE), intent(out)			:: i_next_edge

		integer (kind = BYTE), dimension(3, 3), parameter :: indices = reshape( \
			int([LEFT_EDGE, RIGHT_EDGE, HYPOTENUSE, \
			LEFT_EDGE, HYPOTENUSE, RIGHT_EDGE, \
			HYPOTENUSE, LEFT_EDGE, RIGHT_EDGE], BYTE), [3, 3])

		i_previous_edge = indices(1, cell%i_turtle_type)
		i_color_edge = indices(2, cell%i_turtle_type)
		i_next_edge = indices(3, cell%i_turtle_type)
	end subroutine

	elemental subroutine cell_get_node_indices(cell, i_color_node_out, i_transfer_node, i_color_node_in)
		class(fine_triangle), intent(in)		:: cell
		integer (kind = BYTE), intent(out)			:: i_color_node_out
		integer (kind = BYTE), intent(out)			:: i_transfer_node
		integer (kind = BYTE), intent(out)			:: i_color_node_in

		integer (kind = BYTE), dimension(3, 3), parameter :: indices = reshape( \
			int([TOP_NODE, LEFT_NODE, RIGHT_NODE, \
			LEFT_NODE, TOP_NODE, RIGHT_NODE, \
			LEFT_NODE, RIGHT_NODE, TOP_NODE], BYTE), [3, 3])

		i_color_node_out = indices(1, cell%i_turtle_type)
		i_transfer_node = indices(2, cell%i_turtle_type)
		i_color_node_in = indices(3, cell%i_turtle_type)
	end subroutine

	elemental subroutine cell_get_edge_types(cell, i_previous_edge_type, i_color_edge_type, i_next_edge_type)
		class(fine_triangle), intent(in)		:: cell
		integer (kind = BYTE), intent(out)			:: i_previous_edge_type
		integer (kind = BYTE), intent(out)			:: i_color_edge_type
		integer (kind = BYTE), intent(out)			:: i_next_edge_type

		i_previous_edge_type = OLD
		i_color_edge_type = 0
		i_next_edge_type = NEW
		call mvbits(cell%i_entity_types, 3, 1, i_previous_edge_type, 1)
		call mvbits(cell%i_entity_types, 2, 1, i_next_edge_type, 1)
		call mvbits(cell%i_entity_types, 0, 2, i_color_edge_type, 0)
	end subroutine

	elemental function cell_get_previous_edge_type(cell) result(i_previous_edge_type)
		class(fine_triangle), intent(in)		:: cell
		integer (kind = BYTE)						:: i_previous_edge_type

		i_previous_edge_type = OLD
		call mvbits(cell%i_entity_types, 3, 1, i_previous_edge_type, 1)
	end function

	elemental function cell_get_next_edge_type(cell) result(i_next_edge_type)
		class(fine_triangle), intent(in)		:: cell
		integer (kind = BYTE)						:: i_next_edge_type

		i_next_edge_type = NEW
		call mvbits(cell%i_entity_types, 2, 1, i_next_edge_type, 1)
	end function

	elemental function cell_get_color_edge_type(cell) result(i_color_edge_type)
		class(fine_triangle), intent(in)		:: cell
		integer (kind = BYTE)						:: i_color_edge_type

		i_color_edge_type = 0
		call mvbits(cell%i_entity_types, 0, 2, i_color_edge_type, 0)
	end function

	elemental subroutine cell_set_edge_types(cell, i_previous_edge_type, i_color_edge_type, i_next_edge_type)
		class(fine_triangle), intent(inout)		:: cell
		integer (kind = BYTE), intent(in)			:: i_previous_edge_type
		integer (kind = BYTE), intent(in)			:: i_color_edge_type
		integer (kind = BYTE), intent(in)			:: i_next_edge_type

		cell%i_entity_types = 0
		call mvbits(i_previous_edge_type, 1, 1, cell%i_entity_types, 3)
		call mvbits(i_next_edge_type, 1, 1, cell%i_entity_types, 2)
		call mvbits(i_color_edge_type, 0, 2, cell%i_entity_types, 0)
	end subroutine

	elemental subroutine cell_set_previous_edge_type(cell, i_previous_edge_type)
		class(fine_triangle), intent(inout)		:: cell
		integer (kind = BYTE), intent(in)			:: i_previous_edge_type

		call mvbits(i_previous_edge_type, 1, 1, cell%i_entity_types, 3)
	end subroutine

	elemental subroutine cell_set_next_edge_type(cell, i_next_edge_type)
		class(fine_triangle), intent(inout)		:: cell
		integer (kind = BYTE), intent(in)	    :: i_next_edge_type

		call mvbits(i_next_edge_type, 1, 1, cell%i_entity_types, 2)
	end subroutine

	elemental subroutine cell_set_color_edge_type(cell, i_color_edge_type)
		class(fine_triangle), intent(inout)		:: cell
		integer (kind = BYTE), intent(in)	    :: i_color_edge_type

		call mvbits(i_color_edge_type, 0, 2, cell%i_entity_types, 0)
	end subroutine

	elemental subroutine cell_reverse(cell)
		class(fine_triangle), intent(inout)		:: cell

        !invert edge types: swap previous and next edge and invert old/new bits
        cell%i_entity_types = ishftc(cell%i_entity_types, 1, 4)
        cell%i_entity_types = ishftc(cell%i_entity_types, -1, 3)
        cell%i_entity_types = ieor(cell%i_entity_types, 1_BYTE)

        !invert turtle type: (K, V, H) -> (H, V, K)
        cell%i_turtle_type = 4_BYTE - cell%i_turtle_type

        !invert plotter type
        cell%i_plotter_type = -cell%i_plotter_type
	end subroutine

	elemental subroutine cell_reverse_inner(cell)
		class(fine_triangle), intent(inout)		:: cell

        !invert edge types: invert old/new color edge bit
        cell%i_entity_types = ieor(cell%i_entity_types, 1_BYTE)

        !invert turtle type: (K, V, H) -> (H, V, K)
        cell%i_turtle_type = 4_BYTE - cell%i_turtle_type

        !invert plotter type
        cell%i_plotter_type = -cell%i_plotter_type
	end subroutine

	elemental subroutine cell_reverse_refinement(cell)
		class(fine_triangle), intent(inout)		:: cell
		integer (kind = BYTE), parameter        :: flip_2_3(-1:4) = [-1_BYTE, 0_BYTE, 1_BYTE, 3_BYTE, 2_BYTE, 4_BYTE]

        !invert refinement flag: flip 2 and 3
        cell%refinement = flip_2_3(cell%refinement)
	end subroutine

	!> returns the volume of an element
	elemental function cell_get_scaling(cell) result(scaling)
		class(fine_triangle), intent(in)		:: cell
        integer (kind = BYTE)                      :: i
		real (kind = GRID_SR)					:: scaling
        real (kind = GRID_SR), parameter, dimension(-1 : MAX_DEPTH)		:: r_scalings = [ (0.5_GRID_SR ** i, 0.5_GRID_SR ** i, i = 0, MAX_DEPTH/2) ]

		scaling = r_scalings(cell%i_depth)
	end function

	!> returns the volume of an element
	elemental function cell_get_volume(cell) result(volume)
		class(fine_triangle), intent(in)		:: cell
        integer (kind = BYTE)                      :: i
		real (kind = GRID_SR)					:: volume
        real (kind = GRID_SR), parameter, dimension(0 : MAX_DEPTH)		:: r_volumes = [ (0.5_GRID_SR ** (i + 1), i = 0, MAX_DEPTH) ]

		volume = r_volumes(cell%i_depth)
	end function

	!> returns the leg size of an element
	elemental function cell_get_leg_size(cell) result(edge_size)
		class(fine_triangle), intent(in)		:: cell
		real (kind = GRID_SR)					:: edge_size

        integer (kind = BYTE)                      :: i
        real (kind = GRID_SR), parameter, dimension(0 : MAX_DEPTH)	:: r_leg_sizes = [ (sqrt(0.5_GRID_SR) ** i, i = 0, MAX_DEPTH) ]

		edge_size = r_leg_sizes(cell%i_depth)
	end function

	!> returns the leg size of an element
	elemental function cell_get_hypo_size(cell) result(edge_size)
		class(fine_triangle), intent(in)		:: cell
		real (kind = GRID_SR)					:: edge_size

        integer (kind = BYTE)                      :: i
        real (kind = GRID_SR), parameter, dimension(0 : MAX_DEPTH)	:: r_leg_sizes = [ (sqrt(0.5_GRID_SR) ** i, i = 0, MAX_DEPTH) ]

		edge_size = r_leg_sizes(cell%i_depth - 1)
	end function

	!> returns the edge size of an element
	pure function cell_get_edge_sizes(cell) result(edge_sizes)
		class(fine_triangle), intent(in)		:: cell
		real (kind = GRID_SR)           		:: edge_sizes(3)

        integer (kind = BYTE)                      :: i
        real (kind = GRID_SR), parameter, dimension(-1 : MAX_DEPTH)	:: r_leg_sizes = [ (sqrt(0.5_GRID_SR) ** i, i = -1, MAX_DEPTH) ]

		edge_sizes = [r_leg_sizes(cell%i_depth), r_leg_sizes(cell%i_depth - 1), r_leg_sizes(cell%i_depth)]
	end function

	elemental function get_cell_volume(depth) result(volume)
        integer (kind = BYTE), intent(in)          :: depth
        integer (kind = BYTE)                      :: i
		real (kind = GRID_SR)					:: volume

        real (kind = GRID_SR), parameter, dimension(0 : MAX_DEPTH)		:: r_volumes = [ (0.5_GRID_SR ** (i + 1), i = 0, MAX_DEPTH) ]

		volume = r_volumes(depth)
	end function

	elemental function cell_to_string(cell) result(str)
		class(fine_triangle), intent(in)		:: cell
		character (len = 64)					:: str

		write(str, '(A, 1X, A, 1X, A, 1X, A, 1X, A, 1X, I0)') turtle_type_to_char(cell%i_turtle_type), color_to_char(cell%i_color_edge_color), edge_type_to_char(cell%get_previous_edge_type()), edge_type_to_char(cell%get_color_edge_type()), edge_type_to_char(cell%get_next_edge_type()), cell%i_depth
	end function

    elemental function get_edge_size(depth) result(edge_size)
        integer (kind = BYTE), intent(in)          :: depth
		real (kind = GRID_SR)					:: edge_size

        integer (kind = BYTE)                      :: i
        real (kind = GRID_SR), parameter, dimension(0 : MAX_DEPTH)  :: r_leg_sizes = [ (sqrt(0.5_GRID_SR) ** i, i = 0, MAX_DEPTH) ]

		edge_size = r_leg_sizes(depth)
    end function

    pure function get_edge_sizes(depth) result(edge_sizes)
        integer (kind = BYTE), intent(in)          :: depth
		real (kind = GRID_SR)					:: edge_sizes(3)

		edge_sizes = get_edge_size([depth, depth - 1_1, depth])
    end function

    elemental function encode_edge_size(depth) result(code)
        integer (kind = BYTE), intent(in)          :: depth
		integer (kind = GRID_DI)				:: code

        integer (kind = BYTE)                      :: i
        integer (kind = GRID_DI), parameter, dimension(0 : MAX_DEPTH + 1)  :: codes = [(ishft(1_GRID_DI, MAX_DEPTH / 2 - i), ishft(1_GRID_DI, MAX_DEPTH / 2 - i), i = 0, MAX_DEPTH / 2)]

        code = codes(depth)
    end function

    elemental function decode_distance(code) result(distance)
        integer (kind = GRID_DI), intent(in)    :: code
        real (kind = GRID_SR)                   :: distance

        real (kind = GRID_SR), parameter		:: scaling = 1.0_GRID_SR / real(ishft(1_GRID_DI, MAX_DEPTH / 2), GRID_SR) !< square root of 2

        distance = real(code, GRID_SR) * scaling
    end function

	!**********************************
	!Transformation data initialization
	!**********************************

	subroutine init_transform_data()
		type(t_edge_transform_data), pointer				:: p_edge_data
		integer (kind = BYTE)							    :: i_plotter_type, i, j
		real (kind = GRID_SR)								:: r_angle
		real (kind = GRID_SR)				                :: edge_vectors(2, 3), edge_normals(2, 3)

		!set transformation matrices for the 8 different plotter grammar patterns

        do i_plotter_type = -8, 8
            if (i_plotter_type .ne. 0) then
                associate(p_cell_data => ref_plotter_data(i_plotter_type))
                    !p_cell_data%plotter_type = i_plotter_type

                    !set rotation angle
                    r_angle = PI / 4.0_GRID_SR * real(abs(i_plotter_type), GRID_SR)

                    if (i_plotter_type > 0) then
                        p_cell_data%orientation = 1
                        p_cell_data%jacobian = reshape([cos(r_angle), sin(r_angle), -sin(r_angle), cos(r_angle)], [2, 2])
                    else if (i_plotter_type < 0) then
                        p_cell_data%orientation = -1
                        p_cell_data%jacobian = reshape([-sin(r_angle), cos(r_angle), cos(r_angle), sin(r_angle)], [2, 2])
                    end if

                    !round values to nearest whole numbers
                    !(the matrices should contain only whole numbers, these can be represented as reals without numerical error)
                    forall (i = 1:2, j = 1:2)
                        p_cell_data%jacobian(j, i) = aint(1.5_GRID_SR * p_cell_data%jacobian(j, i), kind=GRID_SR)
                    end forall

                    p_cell_data%det_jacobian = p_cell_data%jacobian(1, 1) * p_cell_data%jacobian(2, 2) - p_cell_data%jacobian(1, 2) * p_cell_data%jacobian(2, 1)
                    p_cell_data%jacobian_inv = 1.0_GRID_SR / p_cell_data%det_jacobian * reshape([ p_cell_data%jacobian(2, 2), -p_cell_data%jacobian(2, 1), -p_cell_data%jacobian(1, 2), p_cell_data%jacobian(1, 1) ], [ 2, 2 ])

                    _log_write(7, '(X, "jacobian: ")')
                    _log_write(7, '(2X, 2(F0.4, X), /, 2X, 2(F0.4, X))') p_cell_data%jacobian
                    _log_write(7, '(X, "inverse: ")')
                    _log_write(7, '(2X, 2(F0.4, X), /, 2X, 2(F0.4, X))') p_cell_data%jacobian_inv
                    _log_write(7, '(X, "determinant: ")')
                    _log_write(7, '(2X, F0.4)') p_cell_data%det_jacobian
                    _log_write(7, '("")')

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
                            p_edge_data%orientation = int(sign(1.0_GRID_SR, dot_product(edge_vectors(:, i), [1.0_GRID_SR, 0.0_GRID_SR])), BYTE)
                        else if (edge_vectors(1, i) == 0.0_GRID_SR) then
                            p_edge_data%orientation = int(sign(1.0_GRID_SR, dot_product(edge_vectors(:, i), [0.0_GRID_SR, 1.0_GRID_SR])), BYTE)
                        else if (edge_vectors(1, i) - edge_vectors(2, i) == 0.0_GRID_SR) then
                            p_edge_data%orientation = int(sign(1.0_GRID_SR, dot_product(edge_vectors(:, i), [1.0_GRID_SR, 1.0_GRID_SR])), BYTE)
                        else if (edge_vectors(1, i) + edge_vectors(2, i) == 0.0_GRID_SR) then
                            p_edge_data%orientation = int(sign(1.0_GRID_SR, dot_product(edge_vectors(:, i), [1.0_GRID_SR, -1.0_GRID_SR])), BYTE)
                        else
                            !something's wrong
                            assert(.false.)
                        end if

                        edge_vectors(:, i) = real(p_edge_data%orientation, GRID_SR) * edge_vectors(:, i)
                        edge_normals(:, i) = 1.0_GRID_SR / sqrt(dot_product(edge_normals(:, i), edge_normals(:, i))) * edge_normals(:, i)

#	    			    if defined(_USE_SKELETON_OP)
                            p_edge_data%normal = edge_normals(:, i)
#	    			    endif

                        _log_write(7, '(X, A, I0)') "edge: ", p_edge_data%index
                        _log_write(7, '(2X, A, I0)') "orientation: ", p_edge_data%orientation
                        _log_write(7, '(2X, A, 2(F0.4, 1X))') "vector: ", edge_vectors(:, i)
                        _log_write(7, '(2X, A, 2(F0.4, 1X))') "normal: ", edge_normals(:, i)
                    end do
                end associate
            end if
        end do
	end subroutine
end MODULE
