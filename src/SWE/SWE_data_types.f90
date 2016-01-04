! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_SWE)
	MODULE SWE_data_types
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

        integer, PARAMETER :: SR = GRID_SR
        integer, PARAMETER :: SI = GRID_SI
        integer, PARAMETER :: DI = GRID_DI

		real (kind = GRID_SR), parameter					:: g = 9.80665_GRID_SR		!< gravitational constant


		!***********************
		!Entity data
		!***********************

		!> state vector of DoFs, either as absoulte values or updates
		type t_dof_state
			real (kind = GRID_SR)													:: h						!< water change
			real (kind = GRID_SR), dimension(2)										:: p						!< momentum change

            contains

            procedure, pass :: add => dof_state_add
			procedure, pass :: inv => dof_state_inv
			procedure, pass :: scale => dof_state_scale

            generic :: operator(+) => add
            generic :: operator(-) => inv
            generic :: operator(*) => scale
		end type

		!> cell state vector including bathymetry
		type, extends(t_dof_state) :: t_state
			real (kind = GRID_SR)													:: b						!< bathymetry

            contains

            procedure, pass :: add_state => state_add
            generic :: operator(+) => add_state
		end type

		!> update vector
		type, extends(t_dof_state) :: t_update
			real (kind = GRID_SR)													:: max_wave_speed			!< maximum wave speed required to compute the CFL condition

            contains

            procedure, pass :: add_update => update_add
            generic :: operator(+) => add_update
		end type

		!> persistent scenario data on a node
		type num_node_data_pers
			integer (kind = BYTE)															:: dummy					!< no data
		END type num_node_data_pers

		!> persistent scenario data on an edge
		type num_edge_data_pers
			integer (kind = BYTE), dimension(0)											:: dummy					!< no data
		END type num_edge_data_pers

		!> persistent scenario data on a cell
		type num_cell_data_pers
			type(t_state), DIMENSION(_SWE_CELL_SIZE)									:: Q						!< cell status vector
		END type num_cell_data_pers

		!> Cell representation on an edge, this would typically be everything required from a cell to compute the flux function on an edge
		type num_cell_rep
			type(t_state), DIMENSION(_SWE_EDGE_SIZE)									:: Q						!< cell representation
		end type

		!> Cell update, this would typically be a flux function
		type num_cell_update
			type(t_update), DIMENSION(_SWE_EDGE_SIZE)									:: flux						!< cell update
		end type

		!*************************
		!Temporary per-Entity data
		!*************************

		!> temporary scenario data on a node (deleted after each traversal)
		type num_node_data_temp
			integer (kind = BYTE), dimension(0)										:: dummy					!< no data
		END type num_node_data_temp

		!> temporary scenario data on an edge (deleted after each traversal)
		type num_edge_data_temp
			integer (kind = BYTE), dimension(0)										:: dummy					!< no data
		END type num_edge_data_temp

		!> temporary scenario data on a cell (deleted after each traversal)
		type num_cell_data_temp
			integer (kind = BYTE), dimension(0)										:: dummy					!< no data
		END type num_cell_data_temp

		!***********************
		!Global data
		!***********************

		!> Data type for the scenario configuration
		type num_global_data
			real (kind = GRID_SR)							:: r_time					!< simulation time
			real (kind = GRID_SR)							:: r_dt						!< time step
			real (kind = GRID_SR)							:: r_dt_new					!< new time step for the next iteration
		end type

		contains

		!adds two state vectors
		elemental function state_add(Q1, Q2)	result(Q_out)
			class (t_state), intent(in)		:: Q1
			type (t_state), intent(in)		:: Q2
			type (t_state)					:: Q_out

			Q_out = t_state(Q1%h + Q2%h, Q1%p + Q2%p, Q1%b + Q2%b)
		end function

		!adds two update vectors
		elemental function update_add(f1, f2)	result(f_out)
			class (t_update), intent(in)		:: f1
			type (t_update), intent(in)		    :: f2
			type (t_update)					    :: f_out

			f_out = t_update(f1%h + f2%h, f1%p + f2%p, max_wave_speed = max(f1%max_wave_speed, f2%max_wave_speed))
		end function

		!adds two dof state vectors
		elemental function dof_state_add(Q1, Q2)	result(Q_out)
			class (t_dof_state), intent(in)		:: Q1
			type (t_dof_state), intent(in)		:: Q2
			type (t_dof_state)					:: Q_out

			Q_out = t_dof_state(Q1%h + Q2%h, Q1%p + Q2%p)
		end function

		!inverts a dof state vector
		elemental function dof_state_inv(f)	result(f_out)
			class (t_dof_state), intent(in)		:: f
			type (t_dof_state)					:: f_out

			f_out = t_dof_state(-f%h, -f%p)
		end function

		!multiplies a scalar with a dof state vector
		elemental function dof_state_scale(f, s)	result(f_out)
			class (t_dof_state), intent(in)		:: f
			real (kind = GRID_SR), intent(in)		:: s
			type (t_dof_state)					:: f_out

			f_out = t_dof_state(s * f%h, s * f%p)
		end function
	END MODULE SWE_data_types
#endif
