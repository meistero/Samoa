! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

!>Implements a stream for cells
module Cell_stream
	use SFC_data_types

#	define _CNT_DATA_TYPE    		type(t_cell_stream_data)
#	define _CNT_TYPE_NAME    		t_cell_stream

#	include "Tools_stream.f90"
end module

!>Implements a stream for crossed edges
module Crossed_edge_stream
	use SFC_data_types

#	define _CNT_DATA_TYPE    		type(t_crossed_edge_stream_data)
#	define _CNT_TYPE_NAME    		t_crossed_edge_stream

#	include "Tools_stream.f90"
end module

!>Implements a stream for color edges
module Color_edge_stream
	use SFC_data_types

#	define _CNT_DATA_TYPE    		type(t_color_edge_stream_data)
#	define _CNT_TYPE_NAME    		t_color_edge_stream

#	include "Tools_stream.f90"
end module

!>Implements a stack for color edges
module Color_edge_stack
	use SFC_data_types

#	define _CNT_DATA_TYPE    		type(t_edge_data)
#	define _CNT_TYPE_NAME    		t_edge_stack

#	include "Tools_stack.f90"
end module

!>Implements a stream for boundary edges
module Boundary_edge_stream
	use SFC_data_types

#	define _CNT_DATA_TYPE    		type(t_edge_data)
#	define _CNT_TYPE_NAME    		t_boundary_edge_stream

#	include "Tools_stream.f90"
end module

!>Implements a stream for nodes
module Node_stream
	use SFC_data_types

#	define _CNT_DATA_TYPE    		type(t_node_stream_data)
#	define _CNT_TYPE_NAME    		t_node_stream

#	include "Tools_stream.f90"
end module


!>Implements a stack for nodes
module Node_stack
	use SFC_data_types

#	define _CNT_DATA_TYPE    		type(t_node_data)
#	define _CNT_TYPE_NAME    		t_node_stack

#	include "Tools_stack.f90"
end module

!>Implements a stack for boundary nodes
module Boundary_node_stream
	use SFC_data_types

#	define _CNT_DATA_TYPE    		type(t_node_data)
#	define _CNT_TYPE_NAME    		t_boundary_node_stream

#	include "Tools_stream.f90"
end module

!>Implements a stack for integers
module Integer_stack
	use SFC_data_types

#	define _CNT_DATA_TYPE    		integer (kind = GRID_SI)
#	define _CNT_TYPE_NAME    		t_integer_stack

#	include "Tools_stack.f90"
end module

!>Implements a list for neighbors
module Integer_list
	use SFC_data_types

	implicit none

#	define _CNT_DATA_TYPE    		integer (kind = GRID_SI)
#	define _CNT_TYPE_NAME    		t_integer_list

#	include "Tools_list.f90"
end module

!>Implements a list for neighbors
module Neighbor_list
	use SFC_data_types

	implicit none

	type t_comm_interface
		integer                         :: local_rank = -1, neighbor_rank = -1			    !< local and comm process rank
		integer                         :: local_section = -1, neighbor_section = -1	    !< local and comm section index
        integer (kind = GRID_DI)        :: min_distance = 0		                            !< minimum distance in the grid (invariant)
		integer (kind = GRID_SI)        :: i_edges = 0		                                !< number of shared edges
		integer (kind = GRID_SI)        :: i_nodes = 0		                                !< number of shared nodes
		integer                         :: send_requests(2) = MPI_REQUEST_NULL              !< mpi send request
		integer                         :: recv_requests(2) = MPI_REQUEST_NULL              !< mpi receive request
		type(t_edge_data), pointer      :: p_local_edges(:) => null(), p_neighbor_edges(:) => null()
		type(t_node_data), pointer      :: p_local_nodes(:) => null(), p_neighbor_nodes(:) => null()

		contains

		procedure, pass :: create_buffer => comm_interface_create_buffer
		procedure, pass :: destroy_buffer => comm_interface_destroy_buffer
		procedure, pass :: to_string => comm_interface_to_string
 		procedure, pass :: eq => comm_interface_eq

        generic :: operator(.eq.) => eq
	end type

	PUBLIC t_comm_interface

#	define _CNT_DATA_TYPE    		type(t_comm_interface)
#	define _CNT_TYPE_NAME    		t_neighbor_list
#	include "Tools_stream.f90"

	subroutine comm_interface_create_buffer(comm)
		class(t_comm_interface), intent(inout)   :: comm
		integer                                         :: i_error

        if (comm%neighbor_rank .ge. 0 .and. comm%neighbor_rank .ne. rank_MPI) then
            assert(.not. associated(comm%p_neighbor_edges))
            assert(.not. associated(comm%p_neighbor_nodes))

            allocate(comm%p_neighbor_edges(comm%i_edges), stat = i_error); assert_eq(i_error, 0)
            allocate(comm%p_neighbor_nodes(comm%i_nodes), stat = i_error); assert_eq(i_error, 0)
        end if
	end subroutine

	subroutine comm_interface_destroy_buffer(comm)
		class(t_comm_interface), intent(inout)   :: comm
		integer                                         :: i_error

        if (comm%neighbor_rank .ge. 0 .and. comm%neighbor_rank .ne. rank_MPI) then
            if (associated(comm%p_neighbor_edges)) then
                deallocate(comm%p_neighbor_edges, stat = i_error); assert_eq(i_error, 0)
            end if

            if (associated(comm%p_neighbor_nodes)) then
                deallocate(comm%p_neighbor_nodes, stat = i_error); assert_eq(i_error, 0)
            end if
        else
            nullify(comm%p_neighbor_edges)
            nullify(comm%p_neighbor_nodes)
        end if
	end subroutine

	elemental function comm_interface_to_string(comm) result(str)
		class(t_comm_interface), intent(in)		:: comm
		character (len = 128)							:: str

		write(str, '(A, I0, 2X, A, I0, 2X, A, F0.4, 2X, A, I0, 2X, A, I0)') "rank: ", comm%neighbor_rank,  "section: ", comm%neighbor_section, "min. distance: ", decode_distance(comm%min_distance), "edges: ", comm%i_edges, "nodes: ", comm%i_nodes
	end function

	elemental function comm_interface_eq(comm1, comm2)
		class(t_comm_interface), intent(in)		:: comm1, comm2
		logical									:: comm_interface_eq

		comm_interface_eq = (comm1%neighbor_rank .eq. comm2%neighbor_rank) .and. (comm1%neighbor_section .eq. comm2%neighbor_section)
	end function
end module

module Section_info_list
    use SFC_data_types

	type t_grid_info
        integer (kind = GRID_DI)    :: i_cells = 0				                !< cells
		integer (kind = GRID_DI)    :: i_crossed_edges = 0		                !< crossed edges
		integer (kind = GRID_DI)    :: i_color_edges = 0		                !< color edges
		integer (kind = GRID_DI)    :: i_nodes = 0				                !< nodes

		integer (kind = GRID_SI)    :: i_boundary_edges(RED : GREEN) = 0	    !< red and green boundary edges
		integer (kind = GRID_SI)    :: i_boundary_nodes(RED : GREEN) = 0	    !< red and green boundary nodes

		integer (kind = GRID_SI)    :: i_stack_nodes(RED : GREEN) = 0		    !< red and green stack edges
		integer (kind = GRID_SI)    :: i_stack_edges(RED : GREEN) = 0		    !< red and green stack nodes

		integer (kind = GRID_SI)    :: i_comms(RED : GREEN) = 0		            !< red and green comms
        integer (kind = GRID_SI)    :: i_comms_type(OLD : NEW, RED : GREEN) = 0 !< old and new, red and green comms

		contains

		procedure, private, pass :: add => grid_info_add
		procedure, pass :: reduce => grid_info_reduce
		procedure, pass :: estimate_bounds => grid_info_estimate_bounds
		procedure, pass :: print => grid_info_print

        generic :: operator(+) => add
	end type

 	type, extends(t_grid_info) :: t_section_info
		integer (kind = GRID_SI)                                                    :: index	    !< index of the section

		contains

		!procedure, private, pass :: eq => section_info_eq

        !generic :: operator(.eq.) => eq
	end type

#	define _CNT_DATA_TYPE    		type(t_section_info)
#	define _CNT_TYPE_NAME    		t_section_info_list

	PUBLIC t_section_info

#	include "Tools_list.f90"

	elemental function grid_info_add(gi1, gi2) result(gi)
		class(t_grid_info), intent(in)	:: gi1, gi2
		type(t_grid_info)           	:: gi

        gi%i_cells = gi1%i_cells + gi2%i_cells
		gi%i_crossed_edges = gi1%i_crossed_edges + gi2%i_crossed_edges

		gi%i_color_edges = gi1%i_color_edges + gi2%i_color_edges
		gi%i_boundary_edges = gi1%i_boundary_edges + gi2%i_boundary_edges
		gi%i_stack_edges = max(gi1%i_stack_edges, gi2%i_stack_edges)

		gi%i_nodes = gi1%i_nodes + gi2%i_nodes
		gi%i_boundary_nodes = gi1%i_boundary_nodes + gi2%i_boundary_nodes
		gi%i_stack_nodes = max(gi1%i_stack_nodes, gi2%i_stack_nodes)

        !not correct, but an upper bound
        gi%i_comms = gi1%i_comms + gi2%i_comms
        gi%i_comms_type = gi1%i_comms_type + gi2%i_comms_type
    end function

    subroutine grid_info_reduce(s, v, global)
        class(t_grid_info), intent(inout)	:: s
        class(t_grid_info), intent(in)		:: v(:)
        logical, intent(in)                 :: global

        integer                             :: i_color
        integer(kind = GRID_DI)             :: i_cells(size(v))

        !HACK: for some reason, the compiler gets problems if we pass v%i_cells directly to reduce()??
        i_cells = v%i_cells

        call reduce(s%i_cells, i_cells, MPI_SUM, global)
		call reduce(s%i_crossed_edges, v%i_crossed_edges, MPI_SUM, global)
		call reduce(s%i_color_edges, v%i_color_edges, MPI_SUM, global)
		call reduce(s%i_nodes, v%i_nodes, MPI_SUM, global)

		do i_color = RED, GREEN
            call reduce(s%i_boundary_edges(i_color), v%i_boundary_edges(i_color), MPI_SUM, global)
            call reduce(s%i_boundary_nodes(i_color), v%i_boundary_nodes(i_color), MPI_SUM, global)

            call reduce(s%i_stack_edges(i_color), v%i_stack_edges(i_color), MPI_MAX, global)
            call reduce(s%i_stack_nodes(i_color), v%i_stack_nodes(i_color), MPI_MAX, global)

            !not correct, but an upper bound
            call reduce(s%i_comms(i_color), v%i_comms(i_color), MPI_SUM, global)
            call reduce(s%i_comms_type(OLD, i_color), v%i_comms_type(OLD, i_color), MPI_SUM, global)
            call reduce(s%i_comms_type(NEW, i_color), v%i_comms_type(NEW, i_color), MPI_SUM, global)
        end do
    end subroutine

     !> Estimates number of (boundary) edges and nodes from the number of cells and the stack size
    pure subroutine grid_info_estimate_bounds(grid_info)
        class(t_grid_info), intent(inout)   :: grid_info

        !there are only #cells - 1 inner crossed edges, but we need an
        !extra edge for intermediate storage of the last crossed edge,
        !which will be converted to a boundary edge afterwards

        grid_info%i_crossed_edges = grid_info%i_cells

        !estimate the number of boundary edges and nodes by two upper bounds:
        !1) #cells + 2, sharp if each edge is a boundary edge
        !2) 2 * stack size, around sqrt(#cells) in the average case, if(!) the correct stack size is used.
        grid_info%i_boundary_edges = min(2 * grid_info%i_stack_edges, grid_info%i_cells + 2)
        grid_info%i_boundary_nodes = min(2 * grid_info%i_stack_nodes - 1, grid_info%i_cells + 2)

        !#color edges is at most 1/2 #cells, and at least 1/2 #cells - 1/2 #boundary edges + 1
        grid_info%i_color_edges = grid_info%i_cells / 2
        grid_info%i_nodes = grid_info%i_cells / 2
    end subroutine

 	subroutine grid_info_print(grid_info)
		class(t_grid_info), intent(in)	:: grid_info

		_log_write(0, "(A)")			"  Info:"
		_log_write(0, '(A)')			""
		_log_write(0, "(A, I14)")		"  Cells                         :", grid_info%i_cells
		_log_write(0, "(A, I14)")		"  Crossed edges                 :", grid_info%i_crossed_edges

		_log_write(0, "(A, I14)")		"  Inner color edges             :", grid_info%i_color_edges
		_log_write(0, "(A, 2(I14))")	"  Boundary edges (red/green)    :", grid_info%i_boundary_edges
		_log_write(0, "(A, 2(I14))")	"  Stack edges (red/green)       :", grid_info%i_stack_edges
		_log_write(0, "(A, I14)")		"  Inner nodes                   :", grid_info%i_nodes
		_log_write(0, "(A, 2(I14))")	"  Boundary nodes (red/green)    :", grid_info%i_boundary_nodes
		_log_write(0, "(A, 2(I14))")	"  Stack nodes (red/green)       :", grid_info%i_stack_nodes

		_log_write(0, "(A, 2(I14))")	"  Comms (red/green)             :", grid_info%i_comms
		_log_write(0, *) ""
	end subroutine

	elemental function section_info_eq(si1, si2)
		class(t_section_info), intent(in)		    :: si1, si2
		logical										:: section_info_eq

		section_info_eq = (si1%index .eq. si2%index)
	end function
end module

module Grid_thread
	use SFC_data_types

	use Integer_stack
	use Color_edge_stack
	use Node_stack

	implicit none

	type t_grid_thread

		type(t_edge_stack), dimension(RED : GREEN)									:: edges_stack									!< red/green color edge stacks
		type(t_node_stack), dimension(RED : GREEN)									:: nodes_stack									!< red/green node stacks
		type(t_integer_stack), dimension(RED : GREEN)								:: indices_stack								!< red/green stacks for comm cell indices in adaptive traversal

		contains

		procedure, pass :: create => grid_thread_create
		procedure, pass :: destroy => grid_thread_destroy
		procedure, pass :: eq => grid_thread_eq

        generic :: operator(.eq.) => eq
 	end type

	PUBLIC t_grid_thread

#	define _CNT_DATA_TYPE    		type(t_grid_thread)
#	define _CNT_TYPE_NAME    		t_grid_thread_list

#	include "Tools_list.f90"


	elemental function grid_thread_eq(s1, s2)
		class(t_grid_thread), intent(in)		:: s1, s2
		logical									:: grid_thread_eq

		grid_thread_eq = associated(s1%edges_stack(RED)%elements, s2%edges_stack(RED)%elements)
	end function

	!> Adds a thread to a grid
	subroutine grid_thread_create(thread, i_stack_size)
		class(t_grid_thread), intent(inout)		:: thread
		integer (kind = GRID_SI), intent(in)	:: i_stack_size(RED : GREEN)

		integer (kind = GRID_SI)				:: i_color

		do i_color = RED, GREEN
			call thread%edges_stack(i_color)%create(i_stack_size(i_color))
			call thread%nodes_stack(i_color)%create(i_stack_size(i_color))
			call thread%indices_stack(i_color)%create(i_stack_size(i_color))
		end do
	end subroutine

	!> Removes a thread from the grid
	subroutine grid_thread_destroy(thread)
		class(t_grid_thread), intent(inout)			:: thread
		integer (kind = GRID_SI)					:: i_color

		do i_color = RED, GREEN
			call thread%edges_stack(i_color)%destroy()
			call thread%nodes_stack(i_color)%destroy()
			call thread%indices_stack(i_color)%destroy()
		end do
	end subroutine
end module

module Grid_section
	use SFC_data_types

	use Cell_stream
	use Crossed_edge_stream
	use Color_edge_stream
	use Node_stream

	use Boundary_edge_stream
	use Boundary_node_stream

	use Integer_list
	use Neighbor_list
	use Section_info_list

	implicit none

	type, extends(t_global_data) :: t_grid_section
		type(t_statistics)															:: stats

 		integer (kind = GRID_SI)													:: index					                    !< source rank of the section

		type(t_cell_stream)															:: cells										!< cell geometry + pers data + refinement stream
		type(t_crossed_edge_stream)													:: crossed_edges_in, crossed_edges_out			!< crossed edge geometry + pers data stream
		type(t_color_edge_stream)													:: color_edges_in, color_edges_out				!< input and output red/green edge streams
		type(t_node_stream)															:: nodes_in, nodes_out							!< red/green stacks for comm cell indices in adaptive traversal

		type(t_boundary_edge_stream), dimension(RED : GREEN)						:: boundary_edges								!< red/green local boundary edges
		type(t_boundary_node_stream), dimension(RED : GREEN)						:: boundary_nodes								!< red/green local boundary nodes
		type(t_boundary_edge_stream), dimension(OLD : NEW, RED : GREEN)				:: boundary_type_edges							!< red/green and old/new local boundary edges
		type(t_boundary_node_stream), dimension(OLD : NEW, RED : GREEN)				:: boundary_type_nodes							!< red/green and old/new local boundary nodes

		type(t_neighbor_list), dimension(RED : GREEN)						 		:: comms				                        !< red/green communication info (will be filled during adaptive refinement)
		type(t_neighbor_list), dimension(OLD : NEW, RED : GREEN)					:: comms_type				                    !< red/green and old/new communication info (will be filled during adaptive refinement)

		contains

		procedure, pass :: create => grid_section_create
		procedure, pass :: destroy => grid_section_destroy
		procedure, pass :: reset => grid_section_reset
		procedure, pass :: reverse => grid_section_reverse
		procedure, pass :: traverse_empty => grid_section_traverse_empty
		procedure, pass :: estimate_load => grid_section_estimate_load
		procedure, pass :: print => grid_section_print
		procedure, pass :: get_capacity => grid_section_get_capacity
		procedure, pass :: eq => grid_section_eq

        generic :: operator(.eq.) => eq
 	end type

	PUBLIC t_grid_section, t_grid_info, t_section_info, t_section_info_list

#	define _CNT_DATA_TYPE    		type(t_grid_section)
#	define _CNT_TYPE_NAME    		t_grid_section_list

#	include "Tools_list.f90"


	elemental function grid_section_eq(s1, s2)
		class(t_grid_section), intent(in)		:: s1, s2
		logical									:: grid_section_eq

		grid_section_eq = associated(s1%cells%elements, s2%cells%elements)
	end function

	!> Adds a section to a grid
	subroutine grid_section_create(section, section_info)
		class(t_grid_section), intent(inout)		:: section
		type(t_section_info), intent(in)			:: section_info

		integer (kind = GRID_SI)					:: i_color, i_error

        section%index = section_info%index
        section%load = 0.0_GRID_SR

		call section%cells%resize(section_info%i_cells)
		call section%crossed_edges_in%resize(section_info%i_crossed_edges)
		call section%crossed_edges_out%attach(section%crossed_edges_in%elements)

		call section%color_edges_in%resize(section_info%i_color_edges)
		call section%color_edges_out%attach(section%color_edges_in%elements)

		call section%nodes_in%resize(section_info%i_nodes)
		call section%nodes_out%attach(section%nodes_in%elements)

		do i_color = RED, GREEN
			call section%boundary_edges(i_color)%resize(int(section_info%i_boundary_edges(i_color), GRID_DI))
			call section%boundary_nodes(i_color)%resize(int(section_info%i_boundary_nodes(i_color), GRID_DI))
            call section%comms(i_color)%resize(int(section_info%i_comms(i_color), GRID_DI))

			assert_eq(section_info%i_comms(i_color), section_info%i_comms_type(OLD, i_color) + section_info%i_comms_type(NEW, i_color))

			section%comms_type(OLD, i_color)%elements => section%comms(i_color)%elements(1 : section_info%i_comms_type(OLD, i_color))
			section%comms_type(NEW, i_color)%elements => section%comms(i_color)%elements(section_info%i_comms_type(NEW, i_color) + 1 : size(section%comms(i_color)%elements))
		end do
	end subroutine

	!> Removes a section from the grid
	subroutine grid_section_destroy(section)
		class(t_grid_section), intent(inout)			:: section
		integer (kind = GRID_SI)					    :: i, j, i_error

		call section%cells%clear()
		call section%crossed_edges_out%unattach()
		call section%crossed_edges_in%clear()

		call section%color_edges_out%unattach()
		call section%color_edges_in%clear()

		call section%nodes_out%unattach()
		call section%nodes_in%clear()

		do i = RED, GREEN
			call section%boundary_edges(i)%clear()
			call section%boundary_nodes(i)%clear()

			do j = 1, size(section%comms(i)%elements)
                call section%comms(i)%elements(j)%destroy_buffer()
            end do

			call section%comms(i)%clear()
		end do
	end subroutine

    elemental subroutine grid_section_reset(section)
		class(t_grid_section), intent(inout)		        :: section
		integer (kind = 1)                                  :: i_color

		!reverse cell stream
		call section%cells%reset()

		!reverse crossed edge stream
		call section%crossed_edges_in%reset()
		call section%crossed_edges_out%reset()

		!reverse color edge streams
		call section%color_edges_in%reset()
		call section%color_edges_out%reset()
		call section%color_edges_out%reset()
		call section%boundary_edges%reset()
		call section%boundary_type_edges%reset()

		!reverse node streams
		call section%nodes_in%reset()
		call section%nodes_out%reset()
		call section%boundary_nodes%reset()
		call section%boundary_type_nodes%reset()
	end subroutine

	subroutine grid_section_reverse(section)
		class(t_grid_section), intent(inout)	:: section
		type(t_boundary_edge_stream)            :: temp_edge_stream(RED : GREEN)
		type(t_boundary_node_stream)            :: temp_node_stream(RED : GREEN)
		type(t_neighbor_list)                   :: temp_comms(RED : GREEN)
		integer (kind = GRID_DI)                :: temp_distance(RED : GREEN)
		integer (kind = 1)                      :: i_color

        temp_distance = section%start_distance
		section%start_distance = section%end_distance
		section%end_distance = temp_distance

		temp_edge_stream = section%boundary_type_edges(OLD,:)
		section%boundary_type_edges(OLD,:) = section%boundary_type_edges(NEW,:)
		section%boundary_type_edges(NEW,:) = temp_edge_stream

		temp_node_stream = section%boundary_type_nodes(OLD,:)
		section%boundary_type_nodes(OLD,:) = section%boundary_type_nodes(NEW,:)
		section%boundary_type_nodes(NEW,:) = temp_node_stream

		temp_comms = section%comms_type(OLD,:)
		section%comms_type(OLD,:) = section%comms_type(NEW,:)
		section%comms_type(NEW,:) = temp_comms

		!reverse cell stream
		call section%cells%reverse()

		!reverse crossed edge stream
		call section%crossed_edges_in%reverse()
		call section%crossed_edges_out%reverse()

		!reverse color edge streams
		call section%color_edges_in%reverse()
		call section%color_edges_out%reverse()
		call section%boundary_edges%reverse()
		call section%boundary_type_edges%reverse()

		!reverse node streams
		call section%nodes_in%reverse()
		call section%nodes_out%reverse()
		call section%boundary_nodes%reverse()
		call section%boundary_type_nodes%reverse()

		call section%comms%reverse()
		call section%comms_type%reverse()
	end subroutine

	elemental subroutine grid_section_estimate_load(section)
		class(t_grid_section), intent(inout)		        :: section

		!section%load = section%stats%r_computation_time
		section%load = section%dest_cells
    end subroutine

	subroutine grid_section_traverse_empty(section)
		class(t_grid_section), intent(inout)		        :: section

		!reverse cell data
		call section%cells%elements%reverse()

		!reverse section
        call section%reverse()
	end subroutine

	elemental function grid_section_get_capacity(section) result(info)
		class(t_grid_section), intent(in)	    :: section
		type(t_section_info)           	        :: info

        info%index = section%index
        info%i_cells = size(section%cells%elements, kind=GRID_DI)
        info%i_crossed_edges = size(section%crossed_edges_in%elements)

		info%i_color_edges = size(section%color_edges_in%elements)
		info%i_nodes = size(section%nodes_in%elements)
		info%i_boundary_edges = [size(section%boundary_edges(RED)%elements), size(section%boundary_edges(GREEN)%elements)]
		info%i_boundary_nodes = [size(section%boundary_nodes(RED)%elements), size(section%boundary_nodes(GREEN)%elements)]

		info%i_stack_edges = (info%i_boundary_nodes - 1) / 2
		info%i_stack_nodes = (info%i_boundary_nodes - 1) / 2

		info%i_comms = [size(section%comms(RED)%elements), size(section%comms(GREEN)%elements)]
		info%i_comms_type(OLD, :) = [size(section%comms_type(OLD, RED)%elements), size(section%comms_type(OLD, GREEN)%elements)]
		info%i_comms_type(NEW, :) = [size(section%comms_type(NEW, RED)%elements), size(section%comms_type(NEW, GREEN)%elements)]
	end function

	subroutine grid_section_print(section)
		class(t_grid_section), intent(in)	:: section

        integer (KIND = GRID_SI)						:: i, i_color, i_comm

		_log_write(0, "(3X, A)")		"Global data:"
		_log_write(0, '(4X, A)')	    trim(section%to_string())

		_log_write(0, "(3X, A)")		        "Stream/Stack state:"
		_log_write(0, '(4X, A, A)')		        "Cells            ", section%cells%to_string()
		_log_write(0, '(4X, A, A, A)')		    "Crossed Edges    ", section%crossed_edges_in%to_string(), section%crossed_edges_out%to_string()
		_log_write(0, '(4X, A, A, A)')	        "Color Edges      ", section%color_edges_in%to_string(), section%color_edges_out%to_string()
		_log_write(0, '(4X, A, A, A)')	        "Nodes            ", section%nodes_in%to_string(), section%nodes_out%to_string()
		_log_write(0, '(4X, A, 2(A))')	        "Boundary Edges   ", section%boundary_edges%to_string()
		_log_write(0, '(4X, A, 2(A))')	        "Boundary Nodes   ", section%boundary_nodes%to_string()

        _log_write(0, '(3X, A)')	"Neighbors:"
		do i_color = RED, GREEN
            _log_write(0, '(4X, A, A)')	trim(color_to_char(i_color)), ":"

            do i_comm = 1, size(section%comms(i_color)%elements)
                _log_write(0, '(5X, A)') trim(section%comms(i_color)%elements(i_comm)%to_string())
            end do
        end do
	end subroutine
end module

module Grid
	use SFC_data_types
	use Grid_section
	use Grid_thread

    implicit none

    public

	type, extends(t_global_data) :: t_grid
		type(t_grid_section_list)									:: sections
		type(t_grid_thread_list)									:: threads
		type(t_statistics)											:: stats

		contains

		procedure, pass :: create => grid_create
		procedure, pass :: destroy => grid_destroy
		procedure, pass :: print => grid_print

		procedure, pass :: get_capacity => grid_get_capacity
		procedure, pass :: get_local_sections => grid_get_local_sections
		procedure, pass :: get_local_sections_in_traversal_order => grid_get_local_sections_in_traversal_order
		procedure, pass :: reset => grid_reset
		procedure, pass :: reverse => grid_reverse
		procedure, pass :: move => grid_move
	end type

	contains

	!> Allocates all streams and stacks
	subroutine grid_create(grid, dest_grid_sections, i_stack_size)
		class(t_grid), target, intent(inout)		    :: grid
		type(t_section_info_list), intent(inout)	    :: dest_grid_sections
		integer (kind = GRID_SI), intent(in)	        :: i_stack_size(RED : GREEN)

		integer (kind = GRID_SI)                        :: i_section, i_first_local_section, i_last_local_section

        !init distances with 0
        !$omp single
        grid%sections = t_grid_section_list()
        call grid%sections%resize(size(dest_grid_sections%elements))
        call grid%threads%resize(omp_get_max_threads())
        !$omp end single

        call grid%get_local_sections(i_first_local_section, i_last_local_section)

        do i_section = i_first_local_section, i_last_local_section
            call grid%sections%elements_alloc(i_section)%create(dest_grid_sections%elements(i_section))
        end do

        call grid%threads%elements(1 + omp_get_thread_num())%create(i_stack_size)

        !$omp barrier

        !$omp single
        !delete section descriptors
        call dest_grid_sections%clear()
        !$omp end single
	end subroutine

	!> Deallocates all streams and stacks
	subroutine grid_destroy(grid)
		class(t_grid), intent(inout)				:: grid
		integer (KIND = GRID_SI)	                :: i_section, i_first_local_section, i_last_local_section
		integer                                     :: i_error

        if (associated(grid%sections%elements_alloc)) then
            call grid%get_local_sections(i_first_local_section, i_last_local_section)

            do i_section = i_first_local_section, i_last_local_section
                call grid%sections%elements_alloc(i_section)%destroy()
            end do
        end if

        call grid%threads%elements(1 + omp_get_thread_num())%destroy()

        !$omp barrier

        !$omp single
		call grid%sections%clear()
		call grid%threads%clear()
		!$omp end single
	end subroutine

    !< Moves a grid to the destination grid, clearing the source grid in the process
	subroutine grid_move(src_grid, dest_grid)
		class(t_grid), intent(inout)    :: src_grid
		type(t_grid), intent(inout)     :: dest_grid

        dest_grid = src_grid

        src_grid%sections%elements => null()
        src_grid%threads%elements => null()
    end subroutine

	subroutine grid_print(grid)
		class(t_grid), intent(in)	:: grid

		_log_write(0, "(3X, A)")		                "Global data:"
		_log_write(0, '(4X, A)')		                trim(grid%t_global_data%to_string())
	end subroutine

	function grid_get_capacity(grid, global) result(grid_info)
        class(t_grid), intent(in)						    :: grid
        logical, intent(in)                                 :: global
		type(t_section_info)           	                    :: grid_info

        !$omp single
        call grid_info%reduce(grid%sections%elements_alloc%get_capacity(), global)
        !$omp end single copyprivate(grid_info)
	end function

	!t_grid_section

	!> Resets a grid
	elemental subroutine grid_reset(grid)
		class(t_grid), intent(inout)				    :: grid

		call grid%sections%elements%reset()
	end subroutine

	!> Distributes the grid sections among threads **in allocation order** and returns a subset for the local thread
	subroutine grid_get_local_sections(grid, i_first_local_section, i_last_local_section)
		class(t_grid), intent(in)				        :: grid
		integer (kind = GRID_SI), intent(out)			:: i_first_local_section, i_last_local_section
		integer (kind = GRID_SI)			            :: i_thread, i_threads

        !TODO: a more sophisticated scheduler could use the (prefix sum over the)
        !number of cells per section to decide which threads get which sections
        !this is required only if the sections are not of uniform size.

        i_thread = omp_get_thread_num()
        i_threads = omp_get_num_threads()

        i_first_local_section = 1 + (i_thread * size(grid%sections%elements_alloc)) / i_threads
        i_last_local_section = ((1 + i_thread) * size(grid%sections%elements_alloc)) / i_threads

        assert_ge(i_first_local_section, 1)
        assert_le(i_last_local_section, size(grid%sections%elements_alloc))
	end subroutine

	!> Distributes the grid sections among threads **in traversal order** and returns a subset for the local thread
	subroutine grid_get_local_sections_in_traversal_order(grid, i_first_local_section, i_last_local_section)
		class(t_grid), intent(in)				        :: grid
		integer (kind = GRID_SI), intent(out)			:: i_first_local_section, i_last_local_section

		integer (kind = GRID_SI)			            :: i_first_local_section_alloc, i_last_local_section_alloc

        call grid%get_local_sections(i_first_local_section_alloc, i_last_local_section_alloc)

        if (grid%sections%forward) then
            i_first_local_section = i_first_local_section_alloc
            i_last_local_section = i_last_local_section_alloc
        else
            i_first_local_section = size(grid%sections%elements_alloc) + 1 - i_last_local_section_alloc
            i_last_local_section = size(grid%sections%elements_alloc) + 1 - i_first_local_section_alloc
        end if
	end subroutine

	!> Attaches a section to a grid
	subroutine grid_reverse(grid)
		class(t_grid), intent(inout)				        :: grid
		integer (kind = GRID_DI), dimension(RED : GREEN)    :: temp_distance
		integer (kind = GRID_SI)                            :: i_section, i_first_local_section, i_last_local_section

        call grid%get_local_sections(i_first_local_section, i_last_local_section)

		do i_section = i_first_local_section, i_last_local_section
            assert_eq(i_section, grid%sections%elements_alloc(i_section)%index)
            call grid%sections%elements_alloc(i_section)%reverse()
		end do

        !$omp single
        call grid%sections%reverse()

        temp_distance = grid%start_distance
		grid%start_distance = grid%end_distance
		grid%end_distance = temp_distance
		!$omp end single nowait
	end subroutine
end module

