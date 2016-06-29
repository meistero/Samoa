! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE

#include "Compilation_control.f90"

#define _JACOBI								    _solver
#define _JACOBI_USE								_solver_use

#define _PREFIX3(P, X)							_conc3(P,_,X)
#define _T_JACOBI								_PREFIX3(t,_JACOBI)
#define _JACOBI_(X)								_PREFIX3(_JACOBI,X)
#define _T_JACOBI_(X)							_PREFIX3(t,_JACOBI_(X))

#define _gv_size								(3 * _gv_node_size + 3 * _gv_edge_size + _gv_cell_size)

MODULE _JACOBI_(1)
    use SFC_edge_traversal
    use _JACOBI_USE

    implicit none

    type num_traversal_data
        real (kind = GRID_SR)				:: r_sq			!< r^T * r
        real (kind = GRID_SR)				:: alpha        !< update ratio
    end type

    !LSE variables
    type(_gm_A)				        :: gm_A                !< temporary/persistent matrix
    type(_gv_x)						:: gv_x                !< persistent solution

    !if no rhs is defined, we assume rhs = 0
#   if defined(_gv_rhs)
        type(_gv_rhs)			    :: gv_rhs               !< temporary/persistent right hand side
#   endif

    !solver-specific variables
    type(_gv_r)						:: gv_r                 !< temporary/persistent residual
    type(_gv_trace_A)			    :: gv_trace_A           !< temporary/persistent trace of matrix

    !if no Dirichlet boundaries are defined, we assume Neumann boundaries everywhere
#   if defined(_gv_dirichlet)
        type(_gv_dirichlet)		    :: gv_dirichlet
#   endif

#		define _GT_NAME							_T_JACOBI_(traversal)

#		if (_gv_edge_size > 0)
#			define _GT_EDGES
#		endif

#		if (_gv_node_size > 0)
#		    define _GT_NODES
#		endif

#		define _GT_NODES
#		define _GT_NO_COORDS

#		define _GT_PRE_TRAVERSAL_GRID_OP		pre_traversal_grid_op
#		define _GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op
#		define _GT_PRE_TRAVERSAL_OP				pre_traversal_op

#		define _GT_ELEMENT_OP					element_op

#		define _GT_NODE_FIRST_TOUCH_OP			node_first_touch_op
#		define _GT_NODE_LAST_TOUCH_OP			node_last_touch_op
#		define _GT_NODE_REDUCE_OP			    node_reduce_op

        !turn on optimizations if the dirichlet check is a temporary variable
#       if defined(_gv_dirichlet_is_temporary)
#		    define _GT_INNER_NODE_LAST_TOUCH_OP inner_node_last_touch_op
#		    define _GT_INNER_NODE_REDUCE_OP	    inner_node_reduce_op
#       endif

#		define _GT_NODE_MERGE_OP		        node_merge_op

#		include "SFC_generic_traversal_ringbuffer.f90"

    subroutine pre_traversal_grid_op(traversal, grid)
        type(_T_JACOBI_(traversal)), intent(inout)					:: traversal
        type(t_grid), intent(inout)							        :: grid

        call scatter(traversal%alpha, traversal%sections%alpha)
    end subroutine

    subroutine post_traversal_grid_op(traversal, grid)
        type(_T_JACOBI_(traversal)), intent(inout)					:: traversal
        type(t_grid), intent(inout)							        :: grid

        call reduce(traversal%r_sq, traversal%sections%r_sq, MPI_SUM, .true.)
    end subroutine

    subroutine pre_traversal_op(traversal, section)
        type(_T_JACOBI_(traversal)), intent(inout)				:: traversal
        type(t_grid_section), intent(inout)				:: section

        traversal%r_sq = 0.0_GRID_SR
    end subroutine

    !*******************************
    !Geometry operators
    !*******************************

    elemental subroutine node_first_touch_op(traversal, section, node)
        type(_T_JACOBI_(traversal)), intent(in)	    :: traversal
        type(t_grid_section), intent(in)		    :: section
        type(t_node_data), intent(inout)		    :: node

        real(kind = GRID_SR)                        :: r(_gv_node_size)
        real(kind = GRID_SR)                        :: trace_A(_gv_node_size)

        call pre_dof_op(r, trace_A)

        call gv_r%write(node, r)
        call gv_trace_A%write(node, trace_A)
    end subroutine

    subroutine element_op(traversal, section, element)
        type(_T_JACOBI_(traversal)), intent(inout)		:: traversal
        type(t_grid_section), intent(inout)				:: section
        type(t_element_base), intent(inout)		        :: element

        !local variables
        integer :: i
        real(kind = GRID_SR)	:: x(_gv_size), r(_gv_size), trace_A(_gv_size)

        call gv_x%read(element, x)
        call gm_A%apply(element, x, r)
        call gm_A%get_trace(element, trace_A)

        call gv_r%add(element, r)
        call gv_trace_A%add(element, trace_A)
    end subroutine

    elemental subroutine node_last_touch_op(traversal, section, node)
        type(_T_JACOBI_(traversal)), intent(in)			:: traversal
        type(t_grid_section), intent(in)				:: section
        type(t_node_data), intent(inout)				:: node

        logical(kind = GRID_SL) :: is_dirichlet(_gv_node_size)
        real(kind = GRID_SR)    :: dx(_gv_node_size)
        real(kind = GRID_SR)    :: r(_gv_node_size)
        real(kind = GRID_SR)    :: rhs(_gv_node_size)
        real(kind = GRID_SR)    :: trace_A(_gv_node_size)

        call gv_dirichlet%read(node, is_dirichlet)
        call gv_r%read(node, r)
        call gv_trace_A%read(node, trace_A)

#       if defined(_gv_rhs)
            call gv_rhs%read(node, rhs)
#       else
            rhs = 0.0_GRID_SR
#       endif

        call post_dof_op(traversal%alpha, dx, r, rhs, trace_A)

        where (is_dirichlet)
            dx = 0.0_SR
            r = 0.0_SR
        end where

        call gv_x%add(node, dx)
        call gv_r%write(node, r)
    end subroutine

    elemental subroutine inner_node_last_touch_op(traversal, section, node)
        type(_T_JACOBI_(traversal)), intent(in)			:: traversal
        type(t_grid_section), intent(in)				:: section
        type(t_node_data), intent(inout)				:: node

        real(kind = GRID_SR)    :: dx(_gv_node_size)
        real(kind = GRID_SR)    :: r(_gv_node_size)
        real(kind = GRID_SR)    :: rhs(_gv_node_size)
        real(kind = GRID_SR)    :: trace_A(_gv_node_size)

        call gv_r%read(node, r)
        call gv_trace_A%read(node, trace_A)

#       if defined(_gv_rhs)
            call gv_rhs%read(node, rhs)
#       else
            rhs = 0.0_GRID_SR
#       endif

        call post_dof_op(traversal%alpha, dx, r, rhs, trace_A)

        call gv_x%add(node, dx)
        call gv_r%write(node, r)
    end subroutine

    pure subroutine node_reduce_op(traversal, section, node)
        type(_T_JACOBI_(traversal)), intent(inout)  :: traversal
        type(t_grid_section), intent(in)		    :: section
        type(t_node_data), intent(in)				:: node

        logical (kind = GRID_SL)                 :: is_dirichlet(_gv_node_size)
        real (kind = GRID_SR)   :: r(_gv_node_size)
        integer					:: i

        call gv_dirichlet%read(node, is_dirichlet)
        call gv_r%read(node, r)

        do i = 1, _gv_node_size
            if (.not. is_dirichlet(i)) then
                call reduce_dof_op(traversal%r_sq, r(i))
            end if
        end do
    end subroutine

    pure subroutine inner_node_reduce_op(traversal, section, node)
        type(_T_JACOBI_(traversal)), intent(inout)	    :: traversal
        type(t_grid_section), intent(in)			    :: section
        type(t_node_data), intent(in)				    :: node

        real (kind = GRID_SR)   :: r(_gv_node_size)
        integer				    :: i

        call gv_r%read(node, r)

        do i = 1, _gv_node_size
            call reduce_dof_op(traversal%r_sq, r(i))
        end do
    end subroutine

    elemental subroutine node_merge_op(local_node, neighbor_node)
        type(t_node_data), intent(inout)			    :: local_node
        type(t_node_data), intent(in)				    :: neighbor_node

        real (kind = GRID_SR) :: r(_gv_node_size)
        real (kind = GRID_SR) :: trace_A(_gv_node_size)

        call gv_r%read(neighbor_node, r)
        call gv_trace_A%read(neighbor_node, trace_A)

        call gv_r%add(local_node, r)
        call gv_trace_A%add(local_node, trace_A)
    end subroutine

    !*******************************
    !Volume and DoF operators
    !*******************************

    elemental subroutine pre_dof_op(r, trace_A)
        real (kind = GRID_SR), intent(out)			:: r
        real (kind = GRID_SR), intent(out)			:: trace_A

        r = 0.0_GRID_SR
        trace_A = tiny(1.0_GRID_SR)
    end subroutine

    elemental subroutine post_dof_op(alpha, dx, r, rhs, trace_A)
        real(kind = GRID_SR), intent(in)			:: alpha
        real(kind = GRID_SR), intent(out)			:: dx
        real(kind = GRID_SR), intent(inout)			:: r
        real (kind = GRID_SR), intent(in)			:: rhs
        real(kind = GRID_SR), intent(in)			:: trace_A

        !so far, r assembled Ax, so now we set it to D^(-1)(b - Ax)
        r = (rhs - r) / trace_A
        dx = alpha * r
    end subroutine

    pure subroutine reduce_dof_op(r_sq, r)
        real(kind = GRID_SR), intent(inout)			:: r_sq
        real(kind = GRID_SR), intent(in)			:: r

        r_sq = r_sq + (r * r)
    end subroutine
END MODULE

MODULE _JACOBI
    use SFC_edge_traversal

    use linear_solver
    use _JACOBI_(1)

    implicit none

    type, extends(t_linear_solver)      :: _T_JACOBI
        type(_T_JACOBI_(traversal))     :: jacobi

        contains

        procedure, pass :: create
        procedure, pass :: destroy
        procedure, pass :: get_info
        procedure, pass :: set_parameter
        procedure, pass :: solve
        procedure, pass :: reduce_stats
        procedure, pass :: clear_stats
    end type

    private
    public _T_JACOBI

    contains

    subroutine create(solver)
        class(_T_JACOBI), intent(inout)   :: solver

        call solver%jacobi%create()

        call base_create(solver)
    end subroutine

    subroutine destroy(solver)
        class(_T_JACOBI), intent(inout) :: solver

        call solver%jacobi%destroy()

        call base_destroy(solver)
    end subroutine

    function get_info(solver, param_idx) result(r_value)
        class(_T_JACOBI), intent(inout)         :: solver
        integer, intent(in)                     :: param_idx
        real (kind = GRID_SR)                   :: r_value

        r_value = solver%base_get_info(param_idx)
    end function

    subroutine set_parameter(solver, param_idx, r_value)
        class(_T_JACOBI), intent(inout)          :: solver
        integer, intent(in)                     :: param_idx
        real (kind = GRID_SR), intent(in)       :: r_value

        call solver%base_set_parameter(param_idx, r_value)
    end subroutine

    !> Solves a poisson equation using a Jacobi solver
    !> \returns		number of iterations performed
    subroutine solve(solver, grid)
        class(_T_JACOBI), intent(inout)					:: solver
        type(t_grid), intent(inout)							        :: grid

        integer (kind = GRID_SI)									:: i_iteration
        real (kind = GRID_SR)										:: r_t1, r_t2
        real (kind = GRID_SR)										:: r_sq, r_sq_old, max_error_sq

        !$omp master
        _log_write(2, '(2X, "Jacobi solver: max abs error:", ES14.7, ", max rel error:", ES14.7)') solver%abs_error, solver%rel_error
        !$omp end master

        !set step size to some initial value
        solver%jacobi%alpha = 1.0_GRID_SR
        i_iteration = 0

        !do a jacobi step
        call solver%jacobi%traverse(grid)
        r_sq = solver%jacobi%r_sq
        r_sq_old = r_sq

        max_error_sq = max(solver%rel_error * solver%rel_error * r_sq, solver%abs_error * solver%abs_error)

        do
            if (iand(i_iteration, z'3ff') == z'3ff') then
                !$omp master
                _log_write(1, '(A, I0, A, F0.10, A, ES17.10)') "   i: ", i_iteration, ", alpha: ", solver%jacobi%alpha, ", res: ", sqrt(r_sq)
                !$omp end master
            else
                !$omp master
                _log_write(2, '(A, I0, A, F0.10, A, ES17.10)') "   i: ", i_iteration, ", alpha: ", solver%jacobi%alpha, ", res: ", sqrt(r_sq)
                !$omp end master
            end if

            if ((solver%max_iterations .ge. 0 .and. i_iteration .ge. solver%max_iterations) .or. (i_iteration .ge. solver%min_iterations .and. r_sq .le. max_error_sq)) then
                exit
            end if

            !adjust step size
            !$omp single
            if (r_sq > r_sq_old) then
                solver%jacobi%alpha = 0.95_GRID_SR * solver%jacobi%alpha
            else
                solver%jacobi%alpha = solver%jacobi%alpha + 0.0001_GRID_SR
            end if
            !$omp end single


            !do a jacobi step
            r_sq_old = r_sq
            call solver%jacobi%traverse(grid)
            r_sq = solver%jacobi%r_sq

            i_iteration = i_iteration + 1
        end do

        !$omp single
        _log_write(2, '(2X, A, T24, I0)') "Jacobi iterations:", i_iteration
        solver%cur_iterations = i_iteration
        solver%cur_error = r_sq
        !$omp end single
    end subroutine

    subroutine reduce_stats(solver, mpi_op, global)
        class(_T_JACOBI), intent(inout)   :: solver
        integer, intent(in)             :: mpi_op
        logical                         :: global

        call solver%jacobi%reduce_stats(mpi_op, global)
        solver%stats = solver%jacobi%stats
    end subroutine

    subroutine clear_stats(solver)
        class(_T_JACOBI), intent(inout)   :: solver

        call solver%jacobi%clear_stats()
        call solver%stats%clear()
    end subroutine
END MODULE

#undef _solver
#undef _solver_use
#undef _PREFIX
