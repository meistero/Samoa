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
#	define _GV_NODE_SIZE		_DARCY_LAYERS

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

	!**************
	!velocity space
	!**************

#	define _GV_CELL_SIZE		2 * _DARCY_LAYERS
#	define _GV_EDGE_SIZE		0
#	define _GV_NODE_SIZE		0

	MODULE Darcy_gv_u_mod
		use SFC_data_types

#		define _GV_TYPE_NAME		darcy_gv_u
#		define _GV_TYPE				real (kind = GRID_SR)
#		define _GV_NAME				u
#		define _GV_PERSISTENT		1

#		include "Tools_grid_variable.f90"
	END MODULE

	!undefine macros to avoid compiler warnings
#	undef _GV_CELL_SIZE
#	undef _GV_EDGE_SIZE
#	undef _GV_NODE_SIZE

	!**********
	!flow space
	!**********

#	define _GV_CELL_SIZE		    0
#	define _GV_EDGE_SIZE		    0
#	define _GV_NODE_SIZE		    _DARCY_LAYERS

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

		type(darcy_gv_saturation)   :: gv_s
		type(darcy_gv_p)            :: gv_p

		private
		public :: darcy_gm_A

        type darcy_gm_A
            contains

            procedure, pass :: read_from_element

            generic:: read => read_from_element
        end type

        contains

        subroutine read_from_element(gm_A, element, mat)
            class(darcy_gm_A), intent(in)	    :: gm_A
            type(t_element_base), intent(in)    :: element

#           if (_DARCY_LAYERS > 1)
                real (kind = GRID_SR), intent(out)  :: mat(6, 6)
                real (kind = GRID_SR)               :: K_base(2)        !The vector contains horizontal and vertical permeability
                real (kind = GRID_SR)               :: lambda_t(3)

                K_base = element%cell%data_pers%base_permeability

                if (element%transform_data%plotter_data%orientation > 0) then
                    lambda_t = element%cell%data_pers%lambda_t
                else
                    lambda_t(1) = element%cell%data_pers%lambda_t(2)
                    lambda_t(2) = element%cell%data_pers%lambda_t(1)
                    lambda_t(3) = element%cell%data_pers%lambda_t(3)
                end if

                mat = 0.0_SR

                !bottom horizontal contributions
                mat(1:2, 1:2) = mat(1:2, 1:2) + lambda_t(1) * K_base(1) * reshape([0.5_SR, -0.5_SR, -0.5_SR, 0.5_SR], [2, 2])
                mat(2:3, 2:3) = mat(2:3, 2:3) + lambda_t(2) * K_base(1) * reshape([0.5_SR, -0.5_SR, -0.5_SR, 0.5_SR], [2, 2])

                !top horizontal contributions
                mat(4:5, 4:5) = mat(4:5, 4:5) + lambda_t(1) * K_base(1) * reshape([0.5_SR, -0.5_SR, -0.5_SR, 0.5_SR], [2, 2])
                mat(5:6, 5:6) = mat(5:6, 5:6) + lambda_t(2) * K_base(1) * reshape([0.5_SR, -0.5_SR, -0.5_SR, 0.5_SR], [2, 2])

                !vertical contributions
                mat(1:4:3,1:4:3) = mat(1:4:3,1:4:3) + lambda_t(3) * K_base(2) * reshape([0.5_SR, -0.5_SR, -0.5_SR, 0.5_SR], [2, 2])
                mat(2:5:3,2:5:3) = mat(2:5:3,2:5:3) + lambda_t(3) * K_base(2) * reshape([0.5_SR, -0.5_SR, -0.5_SR, 0.5_SR], [2, 2])
                mat(3:6:3,3:6:3) = mat(3:6:3,3:6:3) + lambda_t(3) * K_base(2) * reshape([0.5_SR, -0.5_SR, -0.5_SR, 0.5_SR], [2, 2])

                !Well, we better look at this thing..
                !print '(6(6(ES11.3, X), /))', mat
#           else
                real (kind = GRID_SR), intent(out)  :: mat(3, 3)
                real (kind = GRID_SR)               :: K_base
                real (kind = GRID_SR)               :: lambda_t(2)

                K_base = element%cell%data_pers%base_permeability

                if (element%transform_data%plotter_data%orientation > 0) then
                    lambda_t = element%cell%data_pers%lambda_t
                else
                    lambda_t(1) = element%cell%data_pers%lambda_t(2)
                    lambda_t(2) = element%cell%data_pers%lambda_t(1)
                end if

                mat = 0.0_SR

                mat(1:2, 1:2) = mat(1:2, 1:2) + lambda_t(1) * K_base * reshape([0.5_SR, -0.5_SR, -0.5_SR, 0.5_SR], [2, 2])
                mat(2:3, 2:3) = mat(2:3, 2:3) + lambda_t(2) * K_base * reshape([0.5_SR, -0.5_SR, -0.5_SR, 0.5_SR], [2, 2])

                !Well, we better look at this thing..
                !print '(3(3(ES11.3, X), /))', mat
#           endif
        end subroutine
	END MODULE
#endif
