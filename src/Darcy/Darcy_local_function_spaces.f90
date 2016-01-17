! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_DARCY)
	!**************
	!pressure space
	!**************

#	define _GV_CELL_SIZE		0
#	define _GV_EDGE_SIZE		0
#	define _GV_NODE_SIZE		(_DARCY_LAYERS + 1)

	MODULE Darcy_gv_p_mod
		use SFC_data_types

#		define _GV_TYPE_NAME		darcy_gv_p
#		define _GV_TYPE				real (kind = GRID_SR)
#		define _GV_NAME				p
#		define _GV_PERSISTENT		1

#		include "Tools_grid_variable.f90"
	END MODULE

	MODULE Darcy_gv_r_mod
		use SFC_data_types

#		define _GV_TYPE_NAME		darcy_gv_r
#		define _GV_TYPE				real (kind = GRID_SR)
#		define _GV_NAME				r
#		define _GV_PERSISTENT		1

#		include "Tools_grid_variable.f90"
	END MODULE

	MODULE Darcy_gv_rhs_mod
		use SFC_data_types

#		define _GV_TYPE_NAME		darcy_gv_rhs
#		define _GV_TYPE				real (kind = GRID_SR)
#		define _GV_NAME				rhs
#		define _GV_PERSISTENT		1

#		include "Tools_grid_variable.f90"
	END MODULE

	MODULE Darcy_gv_d_mod
		use SFC_data_types

#		define _GV_TYPE_NAME		darcy_gv_d
#		define _GV_TYPE				real (kind = GRID_SR)
#		define _GV_NAME				d
#		define _GV_PERSISTENT		1

#		include "Tools_grid_variable.f90"
	END MODULE

	MODULE Darcy_gv_A_d_mod
		use SFC_data_types

#		define _GV_TYPE_NAME		darcy_gv_A_d
#		define _GV_TYPE				real (kind = GRID_SR)
#		define _GV_NAME				A_d
#		define _GV_PERSISTENT		1

#		include "Tools_grid_variable.f90"
	END MODULE

	MODULE Darcy_gv_mat_diagonal_mod
		use SFC_data_types

#		define _GV_TYPE_NAME		darcy_gv_mat_diagonal
#		define _GV_TYPE				real (kind = GRID_SR)
#		define _GV_NAME				mat_diagonal
#		define _GV_PERSISTENT		0

#		include "Tools_grid_variable.f90"
	END MODULE

    !> this is a wrapper module for the is_dirichlet variable required for the linear solvers
    !> only reads are allowed
	MODULE darcy_gv_is_dirichlet_mod
		use SFC_data_types

        private
        public :: darcy_gv_is_dirichlet

        type darcy_gv_is_dirichlet
            contains

            procedure, pass :: read_from_node

            generic :: read => read_from_node
        end type

        contains

        !> Fills the array is_dirichlet with true or false, depending on the boundary condition
        pure subroutine read_from_node(gv, node, is_dirichlet)
            class(darcy_gv_is_dirichlet), intent(in)            :: gv
            type(t_node_data), intent(in)                       :: node
            logical (kind = SL), intent(inout)	                :: is_dirichlet(:)

            is_dirichlet = logical(node%data_pers%boundary_condition(1) < 0, SL)
        end subroutine
	END MODULE

	!**********
	!flow space
	!**********

	MODULE Darcy_gv_saturation_mod
		use SFC_data_types

#		define _GV_TYPE_NAME		darcy_gv_saturation
#		define _GV_TYPE				real (kind = GRID_SR)
#		define _GV_NAME				saturation
#		define _GV_PERSISTENT		1

#		include "Tools_grid_variable.f90"
	END MODULE

	MODULE Darcy_gv_flux_mod
		use SFC_data_types

#		define _GV_TYPE_NAME		darcy_gv_flux
#		define _GV_TYPE				real (kind = GRID_SR)
#		define _GV_NAME				flux
#		define _GV_PERSISTENT		0

#		include "Tools_grid_variable.f90"
	END MODULE

	MODULE Darcy_gv_volume_mod
		use SFC_data_types

#		define _GV_TYPE_NAME		darcy_gv_volume
#		define _GV_TYPE				real (kind = GRID_SR)
#		define _GV_NAME				volume
#		define _GV_PERSISTENT		0

#		include "Tools_grid_variable.f90"
	END MODULE

#	undef _GV_NODE_SIZE
#	define _GV_NODE_SIZE		    1

	MODULE Darcy_gv_boundary_condition_mod
		use SFC_data_types

#		define _GV_TYPE_NAME		darcy_gv_boundary_condition
#		define _GV_TYPE				integer (kind = SI)
#		define _GV_NAME				boundary_condition
#		define _GV_PERSISTENT		1
#       define _GV_ADD_OP(x, y)     add_boundary_conditions(x, y)

#		include "Tools_grid_variable.f90"

        !> merge two boundary condition evaluations
        elemental function add_boundary_conditions(bc1, bc2) result(bc)
            integer, intent(in) :: bc1, bc2
            integer             :: bc

            if (bc1 .ne. 0 .and. bc2 .ne. 0) then
                !prioritize production over injection wells to avoid having an inflow without an outflow
                bc = min(bc1, bc2)
            else if(bc1 .ne. 0) then
                bc = bc1
            else
                bc = bc2
            end if
        end function
	END MODULE

#	undef _GV_CELL_SIZE
#	undef _GV_EDGE_SIZE
#	undef _GV_NODE_SIZE


	MODULE Darcy_gm_A_mod
		use SFC_data_types
		use darcy_gv_saturation_mod
		use darcy_gv_p_mod
		use Darcy_gv_boundary_condition_mod
		use Samoa

		implicit none

		type(darcy_gv_saturation)           :: gv_s
		type(darcy_gv_p)                    :: gv_p
		type(Darcy_gv_boundary_condition)   :: gv_boundary_condition

		private
		public :: darcy_gm_A

        type darcy_gm_A
            contains

            procedure, pass :: apply
            procedure, pass :: get_trace
        end type

        contains

        subroutine apply(gm_A, element, x, r)
            class(darcy_gm_A), intent(in)	        :: gm_A
            type(t_element_base), intent(in)        :: element

            integer :: i, boundary_condition(3)

#           if (_DARCY_LAYERS > 0)
                real (kind = GRID_SR), contiguous, intent(in)       :: x(:)
                real (kind = GRID_SR), contiguous, intent(inout)    :: r(:)

                if (element%transform_data%plotter_data%orientation > 0) then
                    call apply3D( &
                        x(1 : _DARCY_LAYERS + 1), x(_DARCY_LAYERS + 1 + 1 : 2 * (_DARCY_LAYERS + 1)), x(2 * (_DARCY_LAYERS + 1) + 1 : 3 * (_DARCY_LAYERS + 1)), &
                        r(1 : _DARCY_LAYERS + 1), r(_DARCY_LAYERS + 1 + 1 : 2 * (_DARCY_LAYERS + 1)), r(2 * (_DARCY_LAYERS + 1) + 1 : 3 * (_DARCY_LAYERS + 1)), &
                        element%cell%data_pers%lambda_t)
                else
                    call apply3D( &
                        x(2 * (_DARCY_LAYERS + 1) + 1 : 3 * (_DARCY_LAYERS + 1)), x(_DARCY_LAYERS + 1 + 1 : 2 * (_DARCY_LAYERS + 1)), x(1 : _DARCY_LAYERS + 1), &
                        r(2 * (_DARCY_LAYERS + 1) + 1 : 3 * (_DARCY_LAYERS + 1)), r(_DARCY_LAYERS + 1 + 1 : 2 * (_DARCY_LAYERS + 1)), r(1 : _DARCY_LAYERS + 1), &
                        element%cell%data_pers%lambda_t)
                end if
#           else
                real (kind = GRID_SR), contiguous, intent(in)       :: x(:)
                real (kind = GRID_SR), contiguous, intent(inout)    :: r(:)

                if (element%transform_data%plotter_data%orientation > 0) then
                    call apply2D(x(1), x(2), x(3), r(1), r(2), r(3), element%cell%data_pers%lambda_t)
                else
                    call apply2D(x(3), x(2), x(1), r(3), r(2), r(1), element%cell%data_pers%lambda_t)
                end if
#           endif

            !> if a pressure condition has been defined for injection wells,
            !> we allow only scalar changes to the DoF column
            !> hence, add up all residuals and pretend they are a single value.

#           if defined(_ASAGI)
#               if defined(_DARCY_INJ_PRESSURE)
                    call gv_boundary_condition%read_from_element(element, boundary_condition)

                    if (any(boundary_condition > 0)) then
                        do i = 1, 3
                            if (boundary_condition(i) > 0) then
                                r((i - 1) * (_DARCY_LAYERS + 1) + 1 : i * (_DARCY_LAYERS + 1)) = &
                                    sum(r((i - 1) * (_DARCY_LAYERS + 1) + 1 : i * (_DARCY_LAYERS + 1))) / (_DARCY_LAYERS + 1)
                            end if
                        end do
                    end if
#               endif
#           endif
        end subroutine

        !> 14 * 3 * #layers DOPS
        !> 7 * #layers DRWS
        subroutine apply3D(x1, x2, x3, r1, r2, r3, lambda_t)
            real (kind = GRID_SR), contiguous, intent(in)       :: x1(:), x2(:), x3(:)
            real (kind = GRID_SR), contiguous, intent(inout)    :: r1(:), r2(:), r3(:)
            real (kind = GRID_SR), contiguous, intent(in)       :: lambda_t(:, :)

            r1 = 0.0_SR
            r2 = 0.0_SR
            r3 = 0.0_SR

            !bottom horizontal contributions
            r1(1:_DARCY_LAYERS) = r1(1:_DARCY_LAYERS) + lambda_t(:, 1) * (x1(1:_DARCY_LAYERS) -  x2(1: _DARCY_LAYERS))
            r2(1:_DARCY_LAYERS) = r2(1:_DARCY_LAYERS) + lambda_t(:, 1) * (x2(1:_DARCY_LAYERS) -  x1(1: _DARCY_LAYERS))
            r3(1:_DARCY_LAYERS) = r3(1:_DARCY_LAYERS) + lambda_t(:, 2) * (x3(1:_DARCY_LAYERS) -  x2(1: _DARCY_LAYERS))
            r2(1:_DARCY_LAYERS) = r2(1:_DARCY_LAYERS) + lambda_t(:, 2) * (x2(1:_DARCY_LAYERS) -  x3(1: _DARCY_LAYERS))

            !vertical contributions
            r1(1 : _DARCY_LAYERS) = r1(1 : _DARCY_LAYERS) + lambda_t(:, 3) * (x1(1 : _DARCY_LAYERS) -  x1(2 : _DARCY_LAYERS + 1))
            r1(2 : _DARCY_LAYERS + 1) = r1(2 : _DARCY_LAYERS + 1) + lambda_t(:, 3) * (x1(2 : _DARCY_LAYERS + 1) -  x1(1 : _DARCY_LAYERS))

            r2(1 : _DARCY_LAYERS) = r2(1 : _DARCY_LAYERS) + lambda_t(:, 4) * (x2(1 : _DARCY_LAYERS) -  x2(2 : _DARCY_LAYERS + 1))
            r2(2 : _DARCY_LAYERS + 1) = r2(2 : _DARCY_LAYERS + 1) + lambda_t(:, 4) * (x2(2 : _DARCY_LAYERS + 1) -  x2(1 : _DARCY_LAYERS))

            r3(1 : _DARCY_LAYERS) = r3(1 : _DARCY_LAYERS) + lambda_t(:, 5) * (x3(1 : _DARCY_LAYERS) -  x3(2 : _DARCY_LAYERS + 1))
            r3(2 : _DARCY_LAYERS + 1) = r3(2 : _DARCY_LAYERS + 1) + lambda_t(:, 5) * (x3(2 : _DARCY_LAYERS + 1) -  x3(1 : _DARCY_LAYERS))

            !top horizontal contributions
            r1(2 : _DARCY_LAYERS + 1) = r1(2 : _DARCY_LAYERS + 1) + lambda_t(:, 6) * (x1(2 : _DARCY_LAYERS + 1) -  x2(2 : _DARCY_LAYERS + 1))
            r2(2 : _DARCY_LAYERS + 1) = r2(2 : _DARCY_LAYERS + 1) + lambda_t(:, 6) * (x2(2 : _DARCY_LAYERS + 1) -  x1(2 : _DARCY_LAYERS + 1))
            r3(2 : _DARCY_LAYERS + 1) = r3(2 : _DARCY_LAYERS + 1) + lambda_t(:, 7) * (x3(2 : _DARCY_LAYERS + 1) -  x2(2 : _DARCY_LAYERS + 1))
            r2(2 : _DARCY_LAYERS + 1) = r2(2 : _DARCY_LAYERS + 1) + lambda_t(:, 7) * (x2(2 : _DARCY_LAYERS + 1) -  x3(2 : _DARCY_LAYERS + 1))
        end subroutine

        !> 4 * 3 = 12 DOPS
        subroutine apply2D(x1, x2, x3, r1, r2, r3, lambda_t)
            real (kind = GRID_SR), intent(in)                   :: x1, x2, x3
            real (kind = GRID_SR), intent(inout)                :: r1, r2, r3
            real (kind = GRID_SR), contiguous, intent(in)       :: lambda_t(:)

            r1 = 0.0_SR
            r2 = 0.0_SR
            r3 = 0.0_SR

            r1 = r1 + lambda_t(1) * (x1 - x2)
            r2 = r2 + lambda_t(1) * (x2 - x1)
            r3 = r3 + lambda_t(2) * (x3 - x2)
            r2 = r2 + lambda_t(2) * (x2 - x3)
        end subroutine

        subroutine get_trace(gm_A, element, d)
            class(darcy_gm_A), intent(in)	    :: gm_A
            type(t_element_base), intent(in)    :: element

            integer :: i, boundary_condition(3)

#           if (_DARCY_LAYERS > 0)
                real (kind = GRID_SR), contiguous, intent(inout)       :: d(:)

                if (element%transform_data%plotter_data%orientation > 0) then
                    call get_trace3D( &
                        d(1 : _DARCY_LAYERS + 1), d(_DARCY_LAYERS + 1 + 1 : 2 * (_DARCY_LAYERS + 1)), d(2 * (_DARCY_LAYERS + 1) + 1 : 3 * (_DARCY_LAYERS + 1)), &
                        element%cell%data_pers%lambda_t)
                else
                    call get_trace3D( &
                        d(2 * (_DARCY_LAYERS + 1) + 1 : 3 * (_DARCY_LAYERS + 1)), d(_DARCY_LAYERS + 1 + 1 : 2 * (_DARCY_LAYERS + 1)), d(1 : _DARCY_LAYERS + 1), &
                        element%cell%data_pers%lambda_t)
                end if
#           else
                real (kind = GRID_SR), contiguous, intent(inout)       :: d(:)

                if (element%transform_data%plotter_data%orientation > 0) then
                    call get_trace2D(d(1), d(2), d(3), element%cell%data_pers%lambda_t)
                else
                    call get_trace2D(d(3), d(2), d(1), element%cell%data_pers%lambda_t)
                end if
#           endif

            !> if a pressure condition has been defined for injection wells,
            !> we allow only scalar changes to the DoF column
            !> hence, add up all traces and pretend they are a single value.

#           if defined(_ASAGI)
#               if defined(_DARCY_INJ_PRESSURE)
                    call gv_boundary_condition%read_from_element(element, boundary_condition)

                    if (any(boundary_condition > 0)) then
                        do i = 1, 3
                            if (boundary_condition(i) > 0) then
                                d((i - 1) * (_DARCY_LAYERS + 1) + 1 : i * (_DARCY_LAYERS + 1)) = &
                                    sum(d((i - 1) * (_DARCY_LAYERS + 1) + 1 : i * (_DARCY_LAYERS + 1))) / (_DARCY_LAYERS + 1)
                            end if
                        end do
                    end if
#               endif
#           endif
        end subroutine

        !> 14 * #layers FLOPS
        subroutine get_trace3D(d1, d2, d3, lambda_t)
            real (kind = GRID_SR), contiguous, intent(inout)    :: d1(:), d2(:), d3(:)
            real (kind = GRID_SR), contiguous, intent(in)       :: lambda_t(:, :)

            d1 = 0.0_SR
            d2 = 0.0_SR
            d3 = 0.0_SR

            !bottom horizontal contributions
            d1(1:_DARCY_LAYERS) = d1(1:_DARCY_LAYERS) + lambda_t(:, 1)
            d2(1:_DARCY_LAYERS) = d2(1:_DARCY_LAYERS) + lambda_t(:, 1)
            d3(1:_DARCY_LAYERS) = d3(1:_DARCY_LAYERS) + lambda_t(:, 2)
            d2(1:_DARCY_LAYERS) = d2(1:_DARCY_LAYERS) + lambda_t(:, 2)

            !vertical contributions
            d1(1 : _DARCY_LAYERS) = d1(1 : _DARCY_LAYERS) + lambda_t(:, 3)
            d1(2 : _DARCY_LAYERS + 1) = d1(2 : _DARCY_LAYERS + 1) + lambda_t(:, 3)

            d2(1 : _DARCY_LAYERS) = d2(1 : _DARCY_LAYERS) + lambda_t(:, 4)
            d2(2 : _DARCY_LAYERS + 1) = d2(2 : _DARCY_LAYERS + 1) + lambda_t(:, 4)

            d3(1 : _DARCY_LAYERS) = d3(1 : _DARCY_LAYERS) + lambda_t(:, 5)
            d3(2 : _DARCY_LAYERS + 1) = d3(2 : _DARCY_LAYERS + 1) + lambda_t(:, 5)

            !top horizontal contributions
            d1(2 : _DARCY_LAYERS + 1) = d1(2 : _DARCY_LAYERS + 1) + lambda_t(:, 6)
            d2(2 : _DARCY_LAYERS + 1) = d2(2 : _DARCY_LAYERS + 1) + lambda_t(:, 6)
            d3(2 : _DARCY_LAYERS + 1) = d3(2 : _DARCY_LAYERS + 1) + lambda_t(:, 7)
            d2(2 : _DARCY_LAYERS + 1) = d2(2 : _DARCY_LAYERS + 1) + lambda_t(:, 7)
        end subroutine

        !> 4 FLOPS
        subroutine get_trace2D(d1, d2, d3, lambda_t)
            real (kind = GRID_SR), intent(inout)    :: d1, d2, d3
            real (kind = GRID_SR), contiguous, intent(in)       :: lambda_t(:)

            d1 = 0.0_SR
            d2 = 0.0_SR
            d3 = 0.0_SR

            d1 = d1 + lambda_t(1)
            d2 = d2 + lambda_t(1)
            d3 = d3 + lambda_t(2)
            d2 = d2 + lambda_t(2)
        end subroutine
	END MODULE
#endif
