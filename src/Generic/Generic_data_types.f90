! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_GENERIC)
	module generic_data_types
        use, intrinsic :: iso_c_binding

		public

		!data precision

		!data precision
#       if defined(_SINGLE_PRECISION)
            integer, PARAMETER :: GRID_SR = kind(1.0e0)
#       elif defined(_DOUBLE_PRECISION)
            integer, PARAMETER :: GRID_SR = kind(1.0d0)
#       elif defined(_QUAD_PRECISION)
            integer, PARAMETER :: GRID_SR = kind(1.0q0)
#       else
#           error "No floating point precision is chosen!"
#       endif

		integer, PARAMETER :: BYTE = selected_int_kind(1)
		integer, PARAMETER :: SHORT = selected_int_kind(4)
		integer, PARAMETER :: GRID_SI = selected_int_kind(8)
		integer, PARAMETER :: GRID_DI = selected_int_kind(16)

        integer, PARAMETER :: SR = GRID_SR
        integer, PARAMETER :: SI = GRID_SI
        integer, PARAMETER :: DI = GRID_DI

		!*********************************************
		!Persistent Entity data (geometric association)
		!**********************************************

		!> persistent, scenario specific data on a node
		type num_node_data_pers
            integer :: index
		end type num_node_data_pers

		!> persistent, scenario specific data on an edge
		type num_edge_data_pers
            integer :: index
		end type num_edge_data_pers

		!> persistent, scenario specific data on a cell
		type num_cell_data_pers
            integer :: index
		end type num_cell_data_pers

		!*********************************************
		!Temporary Entity data (geometric association)
		!*********************************************

		!> temporary, scenario specific data on a node (deleted after each traversal)
		type num_node_data_temp

		end type num_node_data_temp

		!> temporary, scenario specific data on an edge (deleted after each traversal)
		type num_edge_data_temp

		end type num_edge_data_temp

		!> temporary, scenario specific data on a cell (deleted after each traversal)
		type num_cell_data_temp

		end type num_cell_data_temp

		!*************************
		!Global data
		!*************************

		!> Base data type for the scenario configuration
		type num_global_data
            integer (kind = c_long_long)                :: i_cells, i_edges, i_nodes
            integer (kind = c_long_long)                :: i_bcells, i_bedges, i_bnodes
 			integer (kind = c_long_long), allocatable   :: cells_to_edges_map(:, :), cells_to_nodes_map(:, :), edges_to_nodes_map(:, :)
 			integer (kind = c_long_long), allocatable   :: bcells_to_cells_map(:), bedges_to_edges_map(:), bnodes_to_nodes_map(:)
 			real (kind = c_double), allocatable         :: coords(:, :)
		end type
	end module Generic_data_types
#endif
