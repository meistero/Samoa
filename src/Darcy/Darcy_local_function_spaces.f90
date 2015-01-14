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

	MODULE Darcy_gv_is_dirichlet_boundary_mod
		use SFC_data_types

#		define _GV_TYPE_NAME		darcy_gv_is_dirichlet_boundary
#		define _GV_TYPE				logical
#		define _GV_NAME				is_dirichlet_boundary
#		define _GV_PERSISTENT		0
#		define _GV_ADD_OP			.or.

#		include "Tools_grid_variable.f90"
	END MODULE
	!undefine macros to avoid compiler warnings
#	undef _GV_CELL_SIZE
#	undef _GV_EDGE_SIZE
#	undef _GV_NODE_SIZE

	!undefine macros to avoid compiler warnings
#	undef _GV_CELL_SIZE
#	undef _GV_EDGE_SIZE
#	undef _GV_NODE_SIZE

	!**********
	!flow space
	!**********

#	define _GV_CELL_SIZE		    0
#	define _GV_EDGE_SIZE		    0
#	define _GV_NODE_SIZE		    (_DARCY_LAYERS + 1)

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

#	undef _GV_CELL_SIZE
#	undef _GV_EDGE_SIZE
#	undef _GV_NODE_SIZE


	MODULE Darcy_gm_A_mod
		use SFC_data_types
		use darcy_gv_saturation_mod
		use darcy_gv_p_mod
		use Samoa

		implicit none

		type(darcy_gv_saturation)   :: gv_s
		type(darcy_gv_p)            :: gv_p

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

#           if (_DARCY_LAYERS > 0)
                real (kind = GRID_SR), intent(in)       :: x(:)
                real (kind = GRID_SR), intent(inout)    :: r(:)

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
                real (kind = GRID_SR), intent(in)       :: x(:)
                real (kind = GRID_SR), intent(inout)    :: r(:)

                if (element%transform_data%plotter_data%orientation > 0) then
                    call apply2D(x(1), x(2), x(3), r(1), r(2), r(3), element%cell%data_pers%lambda_t)
                else
                    call apply2D(x(3), x(2), x(1), r(3), r(2), r(1), element%cell%data_pers%lambda_t)
                end if
#           endif
        end subroutine

        subroutine apply3D(x1, x2, x3, r1, r2, r3, lambda_t)
            real (kind = GRID_SR), intent(in)       :: x1(:), x2(:), x3(:)
            real (kind = GRID_SR), intent(inout)    :: r1(:), r2(:), r3(:)
            real (kind = GRID_SR), intent(in)       :: lambda_t(:, :)

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

        subroutine apply2D(x1, x2, x3, r1, r2, r3, lambda_t)
            real (kind = GRID_SR), intent(in)       :: x1, x2, x3
            real (kind = GRID_SR), intent(inout)    :: r1, r2, r3
            real (kind = GRID_SR), intent(in)       :: lambda_t(:)

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

#           if (_DARCY_LAYERS > 0)
                real (kind = GRID_SR), intent(inout)       :: d(:)
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
                real (kind = GRID_SR), intent(inout)       :: d(:)

                if (element%transform_data%plotter_data%orientation > 0) then
                    call get_trace2D(d(1), d(2), d(3), element%cell%data_pers%lambda_t)
                else
                    call get_trace2D(d(3), d(2), d(1), element%cell%data_pers%lambda_t)
                end if
#           endif
        end subroutine

        subroutine get_trace3D(d1, d2, d3, lambda_t)
            real (kind = GRID_SR), intent(inout)    :: d1(:), d2(:), d3(:)
            real (kind = GRID_SR), intent(in)       :: lambda_t(:, :)

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

        subroutine get_trace2D(d1, d2, d3, lambda_t)
            real (kind = GRID_SR), intent(inout)    :: d1, d2, d3
            real (kind = GRID_SR), intent(in)       :: lambda_t(:)

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
