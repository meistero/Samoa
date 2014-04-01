module linear_solver
    use SFC_edge_traversal

    implicit none

    type, abstract  :: t_linear_solver
        type(t_statistics)              :: stats

        contains

        procedure(i_solve), deferred, pass :: solve
    end type

    interface
        function i_solve(solver, grid) result(i_iterations)
            import
            class(t_linear_solver), intent(inout)   :: solver
            type(t_grid), intent(inout)			    :: grid
            integer (kind = GRID_SI)	            :: i_iterations
        end function
    end interface
end module
