! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_NUMA)
	MODULE NUMA_data_types
		implicit none

		PUBLIC

		!data precision

		integer, PARAMETER :: GRID_SR = selected_real_kind(14,40)
		integer, PARAMETER :: GRID_DR = selected_real_kind(28,80)

		integer, PARAMETER :: GRID_SI = selected_int_kind(8)
		integer, PARAMETER :: GRID_DI = selected_int_kind(16)

		integer, PARAMETER :: GRID_SL = 1
 
		real (kind = GRID_SR), parameter					:: g = 9.80665_GRID_SR		!< gravitational constant


		!***********************
		!Global data
		!***********************

		!> Data type for the scenario configuration
		type num_global_data
			real (kind = GRID_SR)							:: r_dt	= 1e-3					!< time step
			real (kind = GRID_SR)							:: u_max					!< maximum wave velocity for cfl condition
			integer (kind = 1)								:: d_max					!< current maximum grid depth

			CHARACTER(64)									:: s_file_stamp				!< output file stamp

			real (kind = GRID_SR)							:: r_time					!< simulation time

 			integer					 						:: afh_displacement			!< asagi file handle to displacement data
 			integer					 						:: afh_bathymetry			!< asagi file handle to bathymetry data
		end type

		!***********************
		!Entity data
		!***********************

		!> state vector of DoFs, either as absoulte values or updates
		type t_dof_state
			real (kind = GRID_SR)												:: h						
			real (kind = GRID_SR), dimension(2)										:: p						
			real (kind = GRID_SR)												:: e
			real (kind = GRID_SR)												:: h_ref	!cell constant initialised in the beginning
			real (kind = GRID_SR), dimension(2)										:: p_ref	!cell constant initialised in the beginning				
			real (kind = GRID_SR)												:: e_ref	!cell constant initialised in the beginning
		end type

		!> cell state vector including bathymetry
		type, extends(t_dof_state) :: t_state
		end type

		!> update vector
		type, extends(t_dof_state) :: t_update
			real (kind = GRID_SR)													:: max_wave_speed			!< maximum wave speed required to compute the CFL condition
		end type

		interface operator(+)
			module procedure update_add
			module procedure state_add
			module procedure dof_state_add
		end interface

		interface operator(-)
			module procedure dof_state_inv
		end interface

		interface operator(*)
			module procedure dof_state_scale
		end interface

		!> persistent scenario data on a node
		type num_node_data_pers
			integer (kind = 1)															:: dummy					!< no data	
		END type num_node_data_pers

		!> persistent scenario data on an edge
		type num_edge_data_pers
			integer (kind = 1), dimension(0)											:: dummy					!< no data
		END type num_edge_data_pers

		!> persistent scenario data on a cell
		type num_cell_data_pers
			type(t_state), DIMENSION(_NUMA_CELL_SIZE)									:: Q						!< cell status vector
		END type num_cell_data_pers

		!> Cell representation on an edge, this would typically be everything required from a cell to compute the flux function on an edge
		type num_cell_rep
			type(t_state), DIMENSION(_NUMA_EDGE_SIZE)									:: Q						!< cell representation
		end type

		!> Cell update, this would typically be a flux function
		type num_cell_update
			type(t_update), DIMENSION(_NUMA_EDGE_SIZE)									:: flux						!< cell update
!			type(t_update), DIMENSION(_NUMA_EDGE_SIZE,_NUM_RUNGE_KUTTA_STAGES)									:: flux						!< cell update
		end type

		!*************************
		!Temporary per-Entity data
		!*************************

		!> temporary scenario data on a node (deleted after each traversal)
		type num_node_data_temp
			integer (kind = 1), dimension(0)										:: dummy					!< no data
		END type num_node_data_temp

		!> temporary scenario data on an edge (deleted after each traversal)
		type num_edge_data_temp
			integer (kind = 1), dimension(0)										:: dummy					!< no data
		END type num_edge_data_temp

		!> temporary scenario data on a cell (deleted after each traversal)
		type num_cell_data_temp
			integer (kind = 1), dimension(0)										:: dummy					!< no data
		END type num_cell_data_temp

		contains

		!adds two state vectors
		elemental function state_add(Q1, Q2)	result(Q_out)
			type (t_state), intent(in)		:: Q1, Q2
			type (t_state)					:: Q_out
			
			Q_out = t_state(Q1%h + Q2%h, Q1%p + Q2%p, Q1%e + Q2%e , Q1%h_ref + Q2%h_ref, Q1%p_ref + Q2%p_ref, Q1%e_ref + Q2%e_ref)
		end function

		!adds two update vectors
		elemental function update_add(f1, f2)	result(f_out)
			type (t_update), intent(in)		:: f1, f2
			type (t_update)					:: f_out
			
			f_out = t_update(f1%h + f2%h, f1%p + f2%p, f1%e + f2%e , f1%h_ref + f2%h_ref, f1%p_ref + f2%p_ref, f1%e_ref + f2%e_ref, f1%max_wave_speed + f2%max_wave_speed)
		end function

		!adds two dof state vectors
		elemental function dof_state_add(Q1, Q2)	result(Q_out)
			type (t_dof_state), intent(in)		:: Q1, Q2
			type (t_dof_state)					:: Q_out
			
			Q_out = t_dof_state(Q1%h + Q2%h, Q1%p + Q2%p, Q1%e + Q2%e , Q1%h_ref+ Q2%h_ref, Q1%p_ref + Q2%p_ref, Q1%e_ref + Q2%e_ref)
		end function

		!inverts a dof state vector
		elemental function dof_state_inv(f)	result(f_out)
			type (t_dof_state), intent(in)		:: f
			type (t_dof_state)					:: f_out
			
			f_out = t_dof_state(-f%h, -f%p, -f%e, -f%h_ref, -f%p_ref, -f%e_ref)
		end function


		elemental function t_update_state_inv(f)	result(f_out)
			type (t_update), intent(in)		:: f
			type (t_update)					:: f_out
			
			f_out = t_update(-f%h, -f%p, -f%e, -f%h_ref, -f%p_ref, -f%e_ref,f%max_wave_speed)
		end function


		!multiplies a scalar with a dof state vector
		elemental function dof_state_scale(s, f)	result(f_out)
			real (kind = GRID_SR), intent(in)		:: s
			type (t_dof_state), intent(in)		:: f
			type (t_dof_state)					:: f_out
			
			f_out = t_dof_state(s * f%h, s * f%p, s * f%e , s * f%h_ref, s * f%p_ref, s * f%e_ref)
		end function
	END MODULE NUMA_data_types
#endif
