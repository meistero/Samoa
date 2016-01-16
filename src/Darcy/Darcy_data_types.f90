! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_DARCY)
	MODULE Darcy_data_types
		implicit none

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
        integer, parameter :: GRID_SL = BYTE

        integer, PARAMETER :: SR = GRID_SR
        integer, PARAMETER :: SI = GRID_SI
        integer, PARAMETER :: DI = GRID_DI
        integer, parameter :: SL = GRID_SL

		!*********************************************
		!Persistent Entity data (geometric association)
		!**********************************************

		!> persistent, scenario specific data on a node
		type num_node_data_pers
            real (kind = GRID_SR)   :: p(_DARCY_LAYERS + 1), rhs(_DARCY_LAYERS + 1)
            real (kind = GRID_SR)   :: A_d(_DARCY_LAYERS + 1), d(_DARCY_LAYERS + 1), r(_DARCY_LAYERS + 1)

            real (kind = GRID_SR)   :: saturation(_DARCY_LAYERS + 1)    !< wetting phase saturation
			integer (kind = SI)     :: boundary_condition(1)
		END type

		!> persistent, scenario specific data on an edge
		type num_edge_data_pers
		END type

		!> persistent, scenario specific data on a cell
		type num_cell_data_pers
#           if (_DARCY_LAYERS > 0)
                real (kind = GRID_SR)   :: base_permeability(_DARCY_LAYERS, 2)      !< horizontal and vertical permeability (diagonal tensor K = [[k_h, 0, 0], [0, k_h, 0], [0, 0, k_v]])
                real (kind = GRID_SR)   :: lambda_t(_DARCY_LAYERS, 7)               !< total mobility in local coordinates
                real (kind = GRID_SR)   :: porosity(_DARCY_LAYERS)                  !< element porosity
#           else
                real (kind = GRID_SR)   :: base_permeability                        !< permeability (uniform tensor K = k * I)
                real (kind = GRID_SR)   :: lambda_t(2)                              !< total mobility in local coordinates
                real (kind = GRID_SR)   :: porosity                                 !< element porosity
#           endif
		END type

		!*********************************************
		!Temporary Entity data (geometric association)
		!*********************************************

		!> temporary, scenario specific data on a node (deleted after each traversal)
		type num_node_data_temp
			real (kind = GRID_SR)       :: mat_diagonal(_DARCY_LAYERS + 1)

			real (kind = GRID_SR)       :: flux(_DARCY_LAYERS + 1)
			real (kind = GRID_SR)		:: volume(_DARCY_LAYERS + 1)
		END type num_node_data_temp

		!> temporary, scenario specific data on an edge (deleted after each traversal)
		type num_edge_data_temp
		END type num_edge_data_temp

		!> temporary, scenario specific data on a cell (deleted after each traversal)
		type num_cell_data_temp
		END type num_cell_data_temp

		!*************************
		!Global data
		!*************************

		!> Base data type for the scenario configuration
		type num_global_data
			real (kind = GRID_SR)       :: r_time					!< simulation time
			real (kind = GRID_SR)       :: r_dt						!< time step
			real (kind = GRID_SR)       :: u_max			        !< maximum velocity

            real (kind = GRID_SR)       :: prod_w(-_DARCY_PRODUCER_WELLS : _DARCY_INJECTOR_WELLS) = 0.0_SR
            real (kind = GRID_SR)       :: prod_n(-_DARCY_PRODUCER_WELLS : _DARCY_INJECTOR_WELLS) = 0.0_SR
            real (kind = GRID_SR)       :: prod_w_acc(-_DARCY_PRODUCER_WELLS : _DARCY_INJECTOR_WELLS) = 0.0_SR
            real (kind = GRID_SR)       :: prod_n_acc(-_DARCY_PRODUCER_WELLS : _DARCY_INJECTOR_WELLS) = 0.0_SR
            real (Kind = GRID_SR)       :: p_bh(_DARCY_INJECTOR_WELLS) = 0.0_SR
		END type
	END MODULE Darcy_data_types
#endif
