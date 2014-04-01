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

		integer, PARAMETER :: GRID_SI = selected_int_kind(8)
		integer, PARAMETER :: GRID_DI = selected_int_kind(16)

		integer, PARAMETER :: GRID_SL = 1

		!*********************************************
		!Persistent Entity data (geometric association)
		!**********************************************

		!> persistent, scenario specific data on a node
		type num_node_data_pers
			real (kind = GRID_SR)   :: p(_DARCY_P_NODE_SIZE)
			real (kind = GRID_SR)   :: A_d(_DARCY_P_NODE_SIZE), d(_DARCY_P_NODE_SIZE), r(_DARCY_P_NODE_SIZE)

			real (kind = GRID_SR)   :: saturation(_DARCY_FLOW_NODE_SIZE)    !< water saturation

#           if (_DARCY_U_NODE_SIZE > 0)
			real (kind = GRID_SR)   :: u(2, _DARCY_U_NODE_SIZE)
#           endif
		END type

		!> persistent, scenario specific data on an edge
		type num_edge_data_pers
#           if (_DARCY_P_EDGE_SIZE > 0)
			real (kind = GRID_SR)   :: p(_DARCY_P_EDGE_SIZE)
			real (kind = GRID_SR)   :: A_d(_DARCY_P_EDGE_SIZE), d(_DARCY_P_EDGE_SIZE), r(_DARCY_P_EDGE_SIZE)
#           endif

#           if (_DARCY_FLOW_EDGE_SIZE > 0)
			real (kind = GRID_SR)   :: saturation(_DARCY_FLOW_EDGE_SIZE)    !< water saturation
#           endif

#           if (_DARCY_U_EDGE_SIZE > 0)
			real (kind = GRID_SR)   :: u(2, _DARCY_U_EDGE_SIZE)		        !< velocity
#           endif
		END type

		!> persistent, scenario specific data on a cell
		type num_cell_data_pers
#           if (_DARCY_P_CELL_SIZE > 0)
			real (kind = GRID_SR)   :: p(_DARCY_P_CELL_SIZE)
			real (kind = GRID_SR)   :: A_d(_DARCY_P_CELL_SIZE), d(_DARCY_P_CELL_SIZE), r(_DARCY_P_CELL_SIZE)
#           endif

#           if (_DARCY_FLOW_CELL_SIZE > 0)
			real (kind = GRID_SR)   :: saturation(_DARCY_FLOW_CELL_SIZE)    !< water saturation
#           endif

#           if (_DARCY_U_CELL_SIZE > 0)
			real (kind = GRID_SR)   :: u(2, _DARCY_U_CELL_SIZE)			    !< velocity
#           endif

			real (kind = GRID_SR)   :: base_permeability
			real (kind = GRID_SR)   :: permeability
		END type

		!*********************************************
		!Temporary Entity data (geometric association)
		!*********************************************

		!> temporary, scenario specific data on a node (deleted after each traversal)
		type num_node_data_temp
			real (kind = GRID_SR), DIMENSION(_DARCY_P_NODE_SIZE)		:: r
			real (kind = GRID_SR), DIMENSION(_DARCY_P_NODE_SIZE)		:: mat_diagonal

			real (kind = GRID_SR), DIMENSION(_DARCY_FLOW_NODE_SIZE)		:: flux
			real (kind = GRID_SR), DIMENSION(_DARCY_FLOW_NODE_SIZE)		:: volume
			logical 									                :: is_dirichlet_boundary(1)
		END type num_node_data_temp

		!> temporary, scenario specific data on an edge (deleted after each traversal)
		type num_edge_data_temp
#           if (_DARCY_P_EDGE_SIZE > 0)
			real (kind = GRID_SR), DIMENSION(_DARCY_P_EDGE_SIZE)		:: r
			real (kind = GRID_SR), DIMENSION(_DARCY_P_EDGE_SIZE)		:: mat_diagonal
#           endif

#           if (_DARCY_FLOW_EDGE_SIZE > 0)
			real (kind = GRID_SR), DIMENSION(_DARCY_FLOW_EDGE_SIZE)		:: flux
			real (kind = GRID_SR), DIMENSION(_DARCY_FLOW_EDGE_SIZE)		:: volume
#           endif
		END type num_edge_data_temp

		!> temporary, scenario specific data on a cell (deleted after each traversal)
		type num_cell_data_temp
#           if (_DARCY_P_CELL_SIZE > 0)
			real (kind = GRID_SR), DIMENSION(_DARCY_P_CELL_SIZE)		:: r
			real (kind = GRID_SR), DIMENSION(_DARCY_P_CELL_SIZE)		:: mat_diagonal
#           endif

#           if (_DARCY_FLOW_CELL_SIZE > 0)
			real (kind = GRID_SR), DIMENSION(_DARCY_FLOW_CELL_SIZE)		:: flux
			real (kind = GRID_SR), DIMENSION(_DARCY_FLOW_CELL_SIZE)		:: volume
#           endif
		END type num_cell_data_temp

		!*************************
		!Global data
		!*************************

		!> Base data type for the scenario configuration
		type num_global_data
			real (kind = GRID_SR)						:: r_time					!< simulation time
			real (kind = GRID_SR)						:: r_dt						!< time step
			real (kind = GRID_SR)						:: u_max			        !< maximum velocity
		END type
	END MODULE Darcy_data_types
#endif
