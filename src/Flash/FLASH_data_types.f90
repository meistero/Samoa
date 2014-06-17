#include "Compilation_control.f90"

#if defined(_FLASH)
  MODULE FLASH_data_types
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

    real (kind = GRID_SR), parameter  :: g = 9.80665_GRID_SR        !< gravitational constant


    !***********************
    !Entity data
    !***********************

    !> state vector of DoFs, either as absoulte values or updates
    type t_dof_state
      !real (kind = GRID_SR), dimension(_TIME_STEP_STAGES)      :: h        !< water change
      !real (kind = GRID_SR), dimension(2,_TIME_STEP_STAGES)    :: p        !< momentum change
      real (kind = GRID_SR)                                     :: h        !< water change
      real (kind = GRID_SR), dimension(2)                       :: p        !< momentum change
      real (kind = GRID_SR)                                     :: h_old    !< old water change
      real (kind = GRID_SR), dimension(2)                       :: p_old    !< old momentum change

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
      real (kind = GRID_SR)                                     :: b        !< bathymetry

      contains

        procedure, pass :: add_state => state_add
        generic :: operator(+) => add_state
    end type

    !> update vector
    type, extends(t_dof_state) :: t_update
      real (kind = GRID_SR)                                     :: max_wave_speed   !< maximum wave speed required to compute the CFL condition

      contains

        procedure, pass :: add_update => update_add
        generic :: operator(+) => add_update
    end type

    !> persistent scenario data on a node
    type num_node_data_pers
      integer (kind = 1)                                        :: dummy    !< no data
    END type num_node_data_pers

    !> persistent scenario data on an edge
    type num_edge_data_pers
      integer (kind = 1), dimension(0)                          :: dummy    !< no data
    END type num_edge_data_pers

    !> persistent scenario data on a cell
    type num_cell_data_pers
      type(t_state), DIMENSION(_FLASH_CELL_SIZE)                :: Q        !< cell status vector
    END type num_cell_data_pers

    !> Cell representation on an edge, this would typically be everything required from a cell to compute the flux function on an edge
    type num_cell_rep
      type(t_state), DIMENSION(_FLASH_CELL_SIZE)                :: Q        !< cell representation
    end type

    !> Cell update, this would typically be a flux function
    type num_cell_update
      type(t_update), DIMENSION(_FLASH_EDGE_SIZE)               :: flux     !< cell update
    end type

    !*************************
    !Temporary per-Entity data
    !*************************

    !> temporary scenario data on a node (deleted after each traversal)
    type num_node_data_temp
      integer (kind = 1), dimension(0)                          :: dummy    !< no data
    END type num_node_data_temp

    !> temporary scenario data on an edge (deleted after each traversal)
    type num_edge_data_temp
      integer (kind = 1), dimension(0)                          :: dummy    !< no data
    END type num_edge_data_temp

    !> temporary scenario data on a cell (deleted after each traversal)
    type num_cell_data_temp
      integer (kind = 1), dimension(0)                          :: dummy    !< no data
    END type num_cell_data_temp

    !***********************
    !Global data
    !***********************

    !> Data type for the scenario configuration
    type num_global_data
      real (kind = GRID_SR)                                     :: r_dt             !< time step
      real (kind = GRID_SR)                                     :: u_max            !< maximum wave velocity for cfl condition
      integer (kind = 1)                                        :: d_max            !< current maximum grid depth

      real (kind = GRID_SR)                                     :: r_time           !< simulation time

      contains

        procedure, pass :: init => num_global_data_init
        procedure, pass :: reduce_num_global_data => num_global_data_reduce

        generic :: reduce => reduce_num_global_data
    end type

    contains

    !adds two state vectors
    elemental function state_add(Q1, Q2)    result(Q_out)
      class (t_state), intent(in)                   :: Q1
      type (t_state), intent(in)		            :: Q2
      type (t_state)                                :: Q_out

      Q_out = t_state(Q1%h + Q2%h, Q1%p + Q2%p, Q1%h_old + Q2%h_old, Q1%p_old + Q2%p_old, Q1%b + Q2%b)
    end function

    !adds two update vectors
    elemental function update_add(f1, f2)   result(f_out)
      class (t_update), intent(in)                    :: f1
      type (t_update), intent(in)                     :: f2
      type (t_update)                                 :: f_out

      f_out = t_update(f1%h + f2%h, f1%p + f2%p, f1%h_old + f2%h_old, f1%p_old + f2%p_old, max_wave_speed = max(f1%max_wave_speed, f2%max_wave_speed))
    end function

    !adds two dof state vectors
    elemental function dof_state_add(Q1, Q2)        result(Q_out)
      class (t_dof_state), intent(in)                 :: Q1
      type (t_dof_state), intent(in)                  :: Q2
      type (t_dof_state)                              :: Q_out

      Q_out = t_dof_state(Q1%h + Q2%h, Q1%p + Q2%p, Q1%h_old + Q2%h_old, Q1%p_old + Q2%p_old)
    end function

    !inverts a dof state vector
    elemental function dof_state_inv(f)     result(f_out)
      class (t_dof_state), intent(in)                 :: f
      type (t_dof_state)                              :: f_out

      f_out = t_dof_state(-f%h, -f%p, -f%h_old, -f%p_old)
    end function

    !multiplies a scalar with a dof state vector
    elemental function dof_state_scale(f, s)        result(f_out)
      class (t_dof_state), intent(in)                 :: f
      real (kind = GRID_SR), intent(in)               :: s
      type (t_dof_state)                              :: f_out

      f_out = t_dof_state(s * f%h, s * f%p, s * f%h_old, s * f%p_old)
    end function


    elemental subroutine num_global_data_init(gd)
      class(num_global_data), intent(inout)           :: gd

      gd%u_max = 0
      gd%d_max = 0
    end subroutine

    elemental subroutine num_global_data_reduce(gd1, gd2)
      class(num_global_data), intent(inout)           :: gd1
      type(num_global_data), intent(in)               :: gd2

      gd1%u_max = max(gd1%u_max, gd2%u_max)
      gd1%d_max = max(gd1%d_max, gd2%d_max)
    end subroutine

  END MODULE FLASH_data_types
#endif
