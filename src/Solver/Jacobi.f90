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
#       define _GT_EDGES
#       define _GT_EDGES_TEMP
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

        !write(*,*) 'qp of element in jacobi:', x
        assert_eq(x(1),x(1))
        assert_eq(x(2),x(2))
        assert_eq(x(3),x(3))

        !write (*,*) 'node 1: ' ,element%nodes(1)%ptr%position(1) ,','  ,element%nodes(1)%ptr%position(2)

         !write (*,*) 'node 2: ' ,element%nodes(2)%ptr%position(1) ,','  ,element%nodes(2)%ptr%position(2)

         !write (*,*) 'node 3: ' ,element%nodes(3)%ptr%position(1) ,','  ,element%nodes(3)%ptr%position(2)
         !if (element%transform_data%plotter_data%orientation /= element%cell%data_pers%original_lse_orientation(1)) then
         !write (*,*) 'orientation changed '
         !endif

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

        procedure, pass :: create
        procedure, pass :: destroy
        procedure, pass :: solve
        procedure, pass :: reduce_stats
        procedure, pass :: clear_stats
    end type

    private
    public _T_JACOBI

    contains

    subroutine create(solver, max_error)
        class(_T_JACOBI), intent(inout) :: solver
        real (kind = GRID_SR)           :: max_error

        solver%max_error = max_error

        call solver%jacobi%create()
    end subroutine

    subroutine destroy(solver)
        class(_T_JACOBI), intent(inout) :: solver

        call solver%jacobi%destroy()
    end subroutine

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
            !do a jacobi step
            call solver%jacobi%traverse(grid)
            r_sq = solver%jacobi%r_sq

            if (iand(i_iteration, z'3ff') == z'3ff') then
                !$omp master
                _log_write(1, '(A, I0, A, F0.10, A, ES17.10)') "   i: ", i_iteration, ", alpha: ", solver%jacobi%alpha, ", res: ", sqrt(r_sq)
                !$omp end master
            else
                !$omp master
                _log_write(2, '(A, I0, A, F0.10, A, ES17.10)') "   i: ", i_iteration, ", alpha: ", solver%jacobi%alpha, ", res: ", sqrt(r_sq)
                !$omp end master
            end if
            write(*,*) 'i: ',i_iteration, ' res:', sqrt(r_sq)

            if (sqrt(r_sq) < solver%max_error) then
                exit
            end if

            !adjust step size
            if (r_sq > r_sq_old) then
                solver%jacobi%alpha = 0.95_GRID_SR * solver%jacobi%alpha
            else
                solver%jacobi%alpha = solver%jacobi%alpha + 0.0001_GRID_SR
            end if

            r_sq_old = r_sq
        end do

        !$omp master
        _log_write(2, '(A, T30, I0)') "  Jacobi iterations:", i_iteration
        !$omp end master
    end function

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
