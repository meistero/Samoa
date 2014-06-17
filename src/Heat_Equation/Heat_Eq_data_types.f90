! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_HEAT_EQ)
	MODULE Heat_Eq_data_types
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

        integer, PARAMETER :: SR = GRID_SR
        integer, PARAMETER :: SI = GRID_SI
        integer, PARAMETER :: DI = GRID_DI

		!***********************
		!Entity data
		!***********************

		!> persistent, Heat_Eq specific data on a node
		type num_node_data_pers
			real (kind = GRID_SR), DIMENSION(_HEAT_EQ_NODE_SIZE)						:: T						!< temperature
			real (kind = GRID_SR), DIMENSION(_HEAT_EQ_NODE_SIZE)						:: T_temp					!< temporary storage for 2-step solvers
		END type num_node_data_pers

		!> persistent, Heat_Eq specific data on an edge
		type num_edge_data_pers
			real (kind = GRID_SR), DIMENSION(_HEAT_EQ_EDGE_SIZE)						:: T						!< temperature
			real (kind = GRID_SR), DIMENSION(_HEAT_EQ_EDGE_SIZE)						:: T_temp					!< temporary storage for 2-step solvers
		END type num_edge_data_pers

		!> persistent, Heat_Eq specific data on a cell
		type num_cell_data_pers
			real (kind = GRID_SR), DIMENSION(_HEAT_EQ_CELL_SIZE)						:: T						!< temperature
			real (kind = GRID_SR), DIMENSION(_HEAT_EQ_CELL_SIZE)						:: T_temp					!< temporary storage for 2-step solvers
			real (kind = GRID_SR)														:: heat_conductivity		!< heat conductivity
		END type num_cell_data_pers

		!*************************
		!Temporary per-Entity data
		!*************************

		!> temporary, Heat_Eq specific data on a node (deleted after each traversal)
		type num_node_data_temp
			real (kind = GRID_SR), DIMENSION(_HEAT_EQ_NODE_SIZE)						:: r						!< residual
			real (kind = GRID_SR), DIMENSION(_HEAT_EQ_NODE_SIZE)						:: mat_mass_diagonal		!< diagonal entry of the global mass matrix for the node
		END type num_node_data_temp

		!> temporary, Heat_Eq specific data on an edge (deleted after each traversal)
		type num_edge_data_temp
			real (kind = GRID_SR), DIMENSION(_HEAT_EQ_EDGE_SIZE)						:: r						!< residual
			real (kind = GRID_SR), DIMENSION(_HEAT_EQ_EDGE_SIZE)						:: mat_mass_diagonal		!< diagonal entry of the global mass matrix for the node
		END type num_edge_data_temp

		!> temporary, Heat_Eq specific data on a cell (deleted after each traversal)
		type num_cell_data_temp
			real (kind = GRID_SR), DIMENSION(_HEAT_EQ_CELL_SIZE)						:: r						!< residual
			real (kind = GRID_SR), DIMENSION(_HEAT_EQ_CELL_SIZE)						:: mat_mass_diagonal		!< diagonal entry of the global mass matrix for the node
		END type num_cell_data_temp

		!***********************
		!Global data
		!***********************

		!> Base data type for the scenario configuration
		type num_global_data
            integer (kind = BYTE)					            :: i_min_depth, i_max_depth !< minimum and maximum grid depth
			integer (kind = GRID_SI)						:: i_iterations				!< number of iterations
			real (kind = GRID_SR)							:: r_dt						!< time step

			real (kind = GRID_SR)							:: r_time					!< simulation time
			real (kind = GRID_SR)							:: r_laser_rps				!< rotations per second of the laser heat source

            contains

            procedure, pass :: init => num_global_data_init
            procedure, pass :: reduce_num_global_data => num_global_data_reduce

            generic :: reduce => reduce_num_global_data
		end type

        contains

		elemental subroutine num_global_data_init(gd)
            class(num_global_data), intent(inout)		:: gd

            !nothing happens here
        end subroutine

		elemental subroutine num_global_data_reduce(gd1, gd2)
            class(num_global_data), intent(inout)	:: gd1
            type(num_global_data), intent(in)	    :: gd2

            !nothing happens here
        end subroutine
	END MODULE Heat_Eq_data_types
#endif
