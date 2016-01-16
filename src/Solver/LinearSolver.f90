! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

module linear_solver
    use SFC_edge_traversal

    implicit none

    !> default parameters for linear solvers, more may be defined by the implementation
    enum, bind(c)
        enumerator  :: LS_ABS_ERROR, LS_REL_ERROR, LS_MIN_ITERS, LS_MAX_ITERS
        enumerator  :: LS_CUR_ERROR, LS_CUR_ITERS
    end enum

    type, abstract  :: t_linear_solver
        type(t_statistics)              :: stats

        real (kind = GRID_SR)           :: rel_error
        real (kind = GRID_SR)           :: abs_error
        real (kind = GRID_SR)           :: cur_error
        integer (kind = GRID_SI)        :: min_iterations
        integer (kind = GRID_SI)        :: max_iterations
        integer (kind = GRID_SI)        :: cur_iterations

        contains

        procedure(i_create), deferred, pass         :: create
        procedure(i_destroy), deferred, pass        :: destroy
        procedure(i_get_info), deferred, pass       :: get_info
        procedure(i_set_parameter), deferred, pass  :: set_parameter
        procedure(i_solve), deferred, pass          :: solve
        procedure(i_reduce_stats), deferred, pass   :: reduce_stats
        procedure(i_clear_stats), deferred, pass    :: clear_stats

        procedure, pass                             :: base_create
        procedure, pass                             :: base_destroy
        procedure, pass                             :: base_get_info
        procedure, pass                             :: base_set_parameter
    end type

    interface
        subroutine i_create(solver)
            import
            class(t_linear_solver), intent(inout)   :: solver
        end subroutine

        subroutine i_destroy(solver)
            import
            class(t_linear_solver), intent(inout)   :: solver
        end subroutine

        subroutine i_set_parameter(solver, param_idx, r_value)
            import
            class(t_linear_solver), intent(inout)   :: solver
            integer, intent(in)                     :: param_idx
            real (kind = GRID_SR), intent(in)       :: r_value
        end subroutine

        function i_get_info(solver, param_idx) result(r_value)
            import
            class(t_linear_solver), intent(inout)   :: solver
            integer, intent(in)                     :: param_idx
            real (kind = GRID_SR)                   :: r_value
        end function

        subroutine i_solve(solver, grid)
            import
            class(t_linear_solver), intent(inout)   :: solver
            type(t_grid), intent(inout)			    :: grid
        end subroutine

        subroutine i_reduce_stats(solver, mpi_op, global)
            import
            class(t_linear_solver), intent(inout)   :: solver
            integer, intent(in)                     :: mpi_op
            logical                                 :: global
        end subroutine

        subroutine i_clear_stats(solver)
            import
            class(t_linear_solver), intent(inout)   :: solver
        end subroutine
    end interface

    contains

    subroutine base_create(solver)
        class(t_linear_solver), intent(inout)   :: solver

        solver%rel_error = epsilon(1.0_GRID_SR)
        solver%abs_error = 0.0_GRID_SR
        solver%cur_error = huge(1.0_GRID_SR)
        solver%min_iterations = 0_GRID_SI
        solver%max_iterations = huge(1_GRID_SI)
        solver%cur_iterations = 0
    end subroutine

    subroutine base_destroy(solver)
        class(t_linear_solver), intent(inout)   :: solver

        !do nothing
    end subroutine

    function base_get_info(solver, param_idx) result(r_value)
        class(t_linear_solver), intent(inout)   :: solver
        integer, intent(in)                     :: param_idx
        real (kind = GRID_SR)                   :: r_value

        select case (param_idx)
            case (LS_CUR_ERROR)
                r_value = sqrt(solver%cur_error)
            case (LS_CUR_ITERS)
                r_value = real(solver%cur_iterations, GRID_SR)
            case default
                assert(.false.)
        end select
    end function

    subroutine base_set_parameter(solver, param_idx, r_value)
        class(t_linear_solver), intent(inout)   :: solver
        integer, intent(in)                     :: param_idx
        real (kind = GRID_SR), intent(in)       :: r_value

        select case (param_idx)
            case (LS_ABS_ERROR)
                solver%abs_error = r_value
            case (LS_REL_ERROR)
                solver%rel_error = r_value
            case (LS_MIN_ITERS)
                solver%min_iterations = int(r_value)
            case (LS_MAX_ITERS)
                solver%max_iterations = int(r_value)
            case default
                assert(.false.)
        end select
    end subroutine
end module
