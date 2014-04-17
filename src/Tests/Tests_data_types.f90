! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_TESTS)
	MODULE Tests_data_types
		PUBLIC

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

		!*********************************************
		!Persistent Entity data (geometric association)
		!**********************************************

		!> persistent, scenario specific data on a node
		type num_node_data_pers
			real (kind = GRID_SR), DIMENSION(2)							:: r_dummy
			integer (kind = GRID_SI)									:: i_ref_cnt
		END type num_node_data_pers

		!> persistent, scenario specific data on an edge
		type num_edge_data_pers
			real (kind = GRID_SR), DIMENSION(2)							:: r_dummy
			integer (kind = GRID_SI)									:: i_ref_cnt
		END type num_edge_data_pers

		!> persistent, scenario specific data on a cell
		type num_cell_data_pers
			real (kind = GRID_SR), DIMENSION(2)							:: r_dummy
		END type num_cell_data_pers

		!*********************************************
		!Temporary Entity data (geometric association)
		!*********************************************

		!> temporary, scenario specific data on a node (deleted after each traversal)
		type num_node_data_temp
			real (kind = GRID_SR), DIMENSION(2)							:: r_dummy
		END type num_node_data_temp

		!> temporary, scenario specific data on an edge (deleted after each traversal)
		type num_edge_data_temp
			real (kind = GRID_SR), DIMENSION(2)							:: r_dummy
		END type num_edge_data_temp

		!> temporary, scenario specific data on a cell (deleted after each traversal)
		type num_cell_data_temp
			real (kind = GRID_SR), DIMENSION(2)							:: r_dummy
		END type num_cell_data_temp

		!*************************
		!Global data
		!*************************

		!> Base data type for the scenario configuration
		type num_global_data
            contains

            procedure, pass :: init => num_global_data_init
            procedure, pass :: reduce_num_global_data => num_global_data_reduce

            generic :: reduce => reduce_num_global_data
		end type


		elemental subroutine num_global_data_init(gd)
            class(num_global_data), intent(inout)		:: gd

            !nothing happens here
        end subroutine

		elemental subroutine num_global_data_reduce(gd1, gd2)
            class(num_global_data), intent(inout)	:: gd1
            type(num_global_data), intent(in)	    :: gd2

            !nothing happens here
        end subroutine
	END MODULE Tests_data_types
#endif

