! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_DARCY)
	!**************
	!pressure space
	!**************

	MODULE Darcy_gm_A_mod
		use SFC_data_types

        type darcy_gm_A
            contains

            procedure, pass :: read_from_element

            generic:: read => read_from_element
        end type

        contains

        pure subroutine read_from_element(gm_A, element, mat)
            class(darcy_gm_A), intent(in)	    :: gm_A
            type(t_element_base), intent(in)    :: element
            real (kind = GRID_SR), intent(out)  :: mat(_DARCY_P_SIZE, _DARCY_P_SIZE)

            real (kind = GRID_SR), parameter    :: mat_const(_DARCY_P_SIZE, _DARCY_P_SIZE) = &
                reshape([ 1.0d0/2.0d0, -1.0d0/2.0d0, 0.0d0, -1.0d0/2.0d0, 1.0d0, -1.0d0/2.0d0, 0.0d0, -1.0d0/2.0d0, 1.0d0/2.0d0], [_DARCY_P_SIZE, _DARCY_P_SIZE])

            mat = element%cell%data_pers%permeability * mat_const
        end subroutine
	END MODULE

#	define _GV_CELL_SIZE		_DARCY_P_CELL_SIZE
#	define _GV_EDGE_SIZE		_DARCY_P_EDGE_SIZE
#	define _GV_NODE_SIZE		_DARCY_P_NODE_SIZE

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

#	define _GV_CELL_SIZE		2 * _DARCY_U_CELL_SIZE
#	define _GV_EDGE_SIZE		2 * _DARCY_U_EDGE_SIZE
#	define _GV_NODE_SIZE		2 * _DARCY_U_NODE_SIZE

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

#	define _GV_CELL_SIZE		    _DARCY_FLOW_CELL_SIZE
#	define _GV_EDGE_SIZE		    _DARCY_FLOW_EDGE_SIZE
#	define _GV_NODE_SIZE		    _DARCY_FLOW_NODE_SIZE

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
#endif
