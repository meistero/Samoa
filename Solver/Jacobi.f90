! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE

#define _JACOBI								    _solver
#define _JACOBI_USE								_solver_use

#define _CONC2(X, Y)							X ## _ ## Y
#define _PREFIX(P, X)							_CONC2(P, X)
#define _T_JACOBI								_PREFIX(t, _JACOBI)
#define _JACOBI_(X)								_PREFIX(_JACOBI, X)
#define _T_JACOBI_(X)							_PREFIX(t, _JACOBI_(X))

#define _gv_size								(3 * _gv_node_size + 3 * _gv_edge_size + _gv_cell_size)

#include "Compilation_control.f90"

MODULE _JACOBI_(1)
    use SFC_edge_traversal
    use _JACOBI_USE

    implicit none

    type num_traversal_data
        real (kind = GRID_SR)				:: r_sq			!< r^T * r
        real (kind = GRID_SR)				:: alpha        !< update ratio
    end type

    !LSE variables
    type(_gm_A)				        :: gm_A
    type(_gv_x)						:: gv_x

    !solver-specific persistent variables
    type(_gv_r)						:: gv_r
    type(_gv_trace_A)			    :: gv_trace_A

    !if no rhs is defined, we assume rhs = 0
#   if defined(_gv_rhs)
        type(_gv_rhs)			    :: gv_rhs
#   endif

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
#		define _GT_INNER_NODE_LAST_TOUCH_OP		inner_node_last_touch_op
#		define _GT_INNER_NODE_REDUCE_OP		    inner_node_reduce_op

#		define _GT_NODE_MERGE_OP		        node_merge_op

#		include "SFC_generic_traversal_ringbuffer.f90"

    subroutine pre_traversal_grid_op(traversal, grid)
        type(_T_JACOBI_(traversal)), intent(inout)					:: traversal
        type(t_grid), intent(inout)							        :: grid

        call scatter(traversal%alpha, traversal%children%alpha)
    end subroutine

    subroutine post_traversal_grid_op(traversal, grid)
        type(_T_JACOBI_(traversal)), intent(inout)					:: traversal
        type(t_grid), intent(inout)							        :: grid

        call reduce(traversal%r_sq, traversal%children%r_sq, MPI_SUM, .true.)
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
        real(kind = GRID_SR)                        :: rhs(_gv_node_size)
        real(kind = GRID_SR)                        :: trace_A(_gv_node_size)

#       if defined(_gv_rhs)
            call gv_rhs%read(node, rhs)
#       else
            rhs = 0.0_GRID_SR
#       endif

        call pre_dof_op(r, rhs, trace_A)

        call gv_r%write(node, r)
        call gv_trace_A%write(node, trace_A)
    end subroutine

    subroutine element_op(traversal, section, element)
        type(_T_JACOBI_(traversal)), intent(inout)		:: traversal
        type(t_grid_section), intent(inout)				:: section
        type(t_element_base), intent(inout), target		:: element

        !local variables
        integer :: i
        real(kind = GRID_SR)	:: x(_gv_size), r(_gv_size), trace_A(_gv_size)
        real(kind = GRID_SR)	:: A(_gv_size, _gv_size)

        call gv_x%read(element, x)
        call gm_A%read(element, A)

        !add up matrix diagonal
        forall (i = 1 : _gv_size)
            trace_A(i) = A(i, i)
        end forall

        r = -matmul(A, x)

        call gv_r%add(element, r)
        call gv_trace_A%add(element, trace_A)
    end subroutine

    elemental subroutine node_last_touch_op(traversal, section, node)
        type(_T_JACOBI_(traversal)), intent(in)			:: traversal
        type(t_grid_section), intent(in)				:: section
        type(t_node_data), intent(inout)				:: node

        logical :: is_dirichlet(1)

        call gv_dirichlet%read(node, is_dirichlet)

        if (.not. any(is_dirichlet)) then
            call inner_node_last_touch_op(traversal, section, node)
        else
            call gv_r%write(node, spread(0.0_GRID_SR, 1, _gv_node_size))
        end if
    end subroutine

    elemental subroutine inner_node_last_touch_op(traversal, section, node)
        type(_T_JACOBI_(traversal)), intent(in)			:: traversal
        type(t_grid_section), intent(in)				:: section
        type(t_node_data), intent(inout)				:: node

        real(kind = GRID_SR)                        :: x(_gv_node_size)
        real(kind = GRID_SR)                        :: r(_gv_node_size)
        real(kind = GRID_SR)                        :: trace_A(_gv_node_size)

        call gv_r%read(node, r)
        call gv_trace_A%read(node, trace_A)

        call post_dof_op(traversal%alpha, x, r, trace_A)

        call gv_x%add(node, x)
        call gv_r%write(node, r)
    end subroutine

    pure subroutine node_reduce_op(traversal, section, node)
        type(_T_JACOBI_(traversal)), intent(inout)  :: traversal
        type(t_grid_section), intent(in)		    :: section
        type(t_node_data), intent(in)				:: node

        logical :: is_dirichlet(1)

        call gv_dirichlet%read(node, is_dirichlet)

        if (.not. any(is_dirichlet)) then
            call inner_node_reduce_op(traversal, section, node)
        end if
    end subroutine

    pure subroutine inner_node_reduce_op(traversal, section, node)
        type(_T_JACOBI_(traversal)), intent(inout)	    :: traversal
        type(t_grid_section), intent(in)			    :: section
        type(t_node_data), intent(in)				    :: node

        integer											:: i

        real (kind = GRID_SR) :: r(_gv_node_size)

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

    elemental subroutine pre_dof_op(r, rhs, trace_A)
        real (kind = GRID_SR), intent(out)			:: r
        real (kind = GRID_SR), intent(in)			:: rhs
        real (kind = GRID_SR), intent(out)			:: trace_A

        r = rhs
        trace_A = tiny(1.0_GRID_SR)
    end subroutine

    elemental subroutine post_dof_op(alpha, x, r, trace_A)
        real(kind = GRID_SR), intent(in)			:: alpha
        real(kind = GRID_SR), intent(out)			:: x
        real(kind = GRID_SR), intent(inout)			:: r
        real(kind = GRID_SR), intent(in)			:: trace_A

        !Jacobi-step

        !add r / trace(A) to the unknown x
        r = r / trace_A
        x = alpha * r
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
        real (kind = GRID_SR)           :: max_error
        type(_T_JACOBI_(traversal))     :: jacobi

        contains

        procedure, pass :: solve
    end type

    interface _T_JACOBI
        module procedure init_solver
    end interface

    private
    public _T_JACOBI

    contains

    function init_solver(max_error) result(solver)
        real (kind = GRID_SR)   :: max_error
        type(_T_JACOBI) :: solver

        solver%max_error = max_error
    end function

    !> Solves a poisson equation using a Jacobi solver
    !> \returns		number of iterations performed
    function solve(solver, grid) result(i_iteration)
        class(_T_JACOBI), intent(inout)					:: solver
        type(t_grid), intent(inout)							        :: grid

        integer (kind = GRID_SI)									:: i_iteration
        real (kind = GRID_SR)										:: r_t1, r_t2
        real (kind = GRID_SR)										:: r_sq, r_sq_old

        !$omp master
        _log_write(3, '(A, ES14.7)') "  Jacobi solver, max residual error:", solver%max_error
        !$omp end master

        !set step size to some initial value
        solver%jacobi%alpha = 1.0_GRID_SR
        r_sq_old = 1.0_GRID_SR

        do i_iteration = 1, huge(1_GRID_SI)
            !$omp master
            _log_write(2, '(A, I0, A, F7.4, A, ES17.10)') "   i: ", i_iteration, ", alpha: ", solver%jacobi%alpha, ", res: ", sqrt(r_sq)
            !$omp end master

            !do a jacobi step
            call solver%jacobi%traverse(grid)
            r_sq = solver%jacobi%r_sq

            !adjust step size
            if (r_sq > r_sq_old) then
                solver%jacobi%alpha = 0.95_GRID_SR * solver%jacobi%alpha
            else
                solver%jacobi%alpha = solver%jacobi%alpha + 0.0001_GRID_SR
            end if

            if (sqrt(r_sq) < solver%max_error) then
                exit
            end if

            r_sq_old = r_sq
        end do

        !$omp master
        _log_write(3, '(A, T30, I0)') "  Jacobi iterations:", i_iteration
        solver%stats = solver%jacobi%stats
        !$omp end master
    end function
END MODULE

#undef _solver
#undef _solver_use
