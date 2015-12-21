module linear_solver
    use SFC_edge_traversal

    implicit none

    !> default parameters for linear solvers, more may be defined by the implementation
    enum, bind(c)
        enumerator  :: LS_MAX_ERROR, LS_MIN_ITERS, LS_RESTART_ITERS, LS_MAX_ITERS
    end enum

    type, abstract  :: t_linear_solver
        type(t_statistics)                          :: stats

        contains

        procedure(i_create), deferred, pass         :: create
        procedure(i_destroy), deferred, pass        :: destroy
        procedure(i_set_parameter), deferred, pass  :: set_parameter
        procedure(i_solve), deferred, pass          :: solve
        procedure(i_reduce_stats), deferred, pass   :: reduce_stats
        procedure(i_clear_stats), deferred, pass    :: clear_stats
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

        function i_solve(solver, grid) result(i_iterations)
            import
            class(t_linear_solver), intent(inout)   :: solver
            type(t_grid), intent(inout)			    :: grid
            integer (kind = GRID_SI)	            :: i_iterations
        end function

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
end module
