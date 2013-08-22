! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_DARCY)
	!**************
	!pressure space
	!**************

#	define _LFS_type_NAME		darcy_lfs_pressure
#	define _LFS_CELL_SIZE		_DARCY_P_CELL_SIZE
#	define _LFS_EDGE_SIZE		_DARCY_P_EDGE_SIZE
#	define _LFS_NODE_SIZE		_DARCY_P_NODE_SIZE

	MODULE Darcy_lfs_pressure_mod
		use SFC_data_types

#		define _LFS_type		real (kind = GRID_SR)
#		include "Tools_local_function_space.f90"
	END MODULE

	MODULE Darcy_gv_p_mod
		use SFC_data_types
		use Darcy_lfs_pressure_mod

#		define _GV_type_NAME		darcy_gv_p
#		define _GV_type				real (kind = GRID_SR)
#		define _GV_NAME				p
#		define _GV_PERSISTENT		.true.

#		include "Tools_grid_variable.f90"
	END MODULE

	MODULE Darcy_gv_r_mod
		use SFC_data_types
		use Darcy_lfs_pressure_mod

#		define _GV_type_NAME		darcy_gv_r
#		define _GV_type				real (kind = GRID_SR)
#		define _GV_NAME				r
#		define _GV_PERSISTENT		.true.

#		include "Tools_grid_variable.f90"
	END MODULE

	MODULE Darcy_gv_d_mod
		use SFC_data_types
		use Darcy_lfs_pressure_mod

#		define _GV_type_NAME		darcy_gv_d
#		define _GV_type				real (kind = GRID_SR)
#		define _GV_NAME				d
#		define _GV_PERSISTENT		.true.

#		include "Tools_grid_variable.f90"
	END MODULE

	MODULE Darcy_gv_A_d_mod
		use SFC_data_types
		use Darcy_lfs_pressure_mod

#		define _GV_type_NAME		darcy_gv_A_d
#		define _GV_type				real (kind = GRID_SR)
#		define _GV_NAME				A_d
#		define _GV_PERSISTENT		.true.

#		include "Tools_grid_variable.f90"
	END MODULE

	MODULE Darcy_gv_r_temp_mod
		use SFC_data_types
		use Darcy_lfs_pressure_mod

#		define _GV_type_NAME		darcy_gv_r_temp
#		define _GV_type				real (kind = GRID_SR)
#		define _GV_NAME				r
#		define _GV_PERSISTENT		.false.

#		include "Tools_grid_variable.f90"
	END MODULE

	MODULE Darcy_gv_mat_diagonal_mod
		use SFC_data_types
		use Darcy_lfs_pressure_mod

#		define _GV_type_NAME		darcy_gv_mat_diagonal
#		define _GV_type				real (kind = GRID_SR)
#		define _GV_NAME				mat_diagonal
#		define _GV_PERSISTENT		.false.

#		include "Tools_grid_variable.f90"
	END MODULE

	!undefine macros to avoid compiler warnings
#	undef _LFS_type_NAME
#	undef _LFS_CELL_SIZE
#	undef _LFS_EDGE_SIZE
#	undef _LFS_NODE_SIZE

	!**************
	!velocity space
	!**************

#	define _LFS_type_NAME		darcy_lfs_velocity
#	define _LFS_CELL_SIZE		_DARCY_U_CELL_SIZE
#	define _LFS_EDGE_SIZE		_DARCY_U_EDGE_SIZE
#	define _LFS_NODE_SIZE		_DARCY_U_NODE_SIZE

	MODULE Darcy_lfs_velocity_mod
		use SFC_data_types

#		define _LFS_type		real (kind = GRID_SR)
#		include "Tools_local_function_space.f90"
	END MODULE

	MODULE Darcy_gv_u_mod
		use SFC_data_types
		use Darcy_lfs_velocity_mod

#		define _GV_type_NAME		darcy_gv_u
#		define _GV_type				real (kind = GRID_SR)
#		define _GV_COUNT			2
#		define _GV_NAME				u
#		define _GV_PERSISTENT		.true.

#		include "Tools_grid_variable.f90"
	END MODULE

	!undefine macros to avoid compiler warnings
#	undef _LFS_type_NAME
#	undef _LFS_CELL_SIZE
#	undef _LFS_EDGE_SIZE
#	undef _LFS_NODE_SIZE

	!**********
	!flow space
	!**********

#	define _LFS_type_NAME		darcy_lfs_flow
#	define _LFS_CELL_SIZE		_DARCY_FLOW_CELL_SIZE
#	define _LFS_EDGE_SIZE		_DARCY_FLOW_EDGE_SIZE
#	define _LFS_NODE_SIZE		_DARCY_FLOW_NODE_SIZE

	MODULE Darcy_lfs_flow_mod
		use SFC_data_types

#		define _LFS_type		real (kind = GRID_SR)
#		include "Tools_local_function_space.f90"
	END MODULE

	MODULE Darcy_gv_saturation_mod
		use SFC_data_types
		use Darcy_lfs_flow_mod

#		define _GV_type_NAME		darcy_gv_saturation
#		define _GV_type				real (kind = GRID_SR)
#		define _GV_NAME				saturation
#		define _GV_PERSISTENT		.true.

#		include "Tools_grid_variable.f90"
	END MODULE

	MODULE Darcy_gv_flux_mod
		use SFC_data_types
		use Darcy_lfs_flow_mod

#		define _GV_type_NAME		darcy_gv_flux
#		define _GV_type				real (kind = GRID_SR)
#		define _GV_NAME				flux
#		define _GV_PERSISTENT		.false.

#		include "Tools_grid_variable.f90"
	END MODULE

	MODULE Darcy_gv_volume_mod
		use SFC_data_types
		use Darcy_lfs_flow_mod

#		define _GV_type_NAME		darcy_gv_volume
#		define _GV_type				real (kind = GRID_SR)
#		define _GV_NAME				volume
#		define _GV_PERSISTENT		.false.

#		include "Tools_grid_variable.f90"
	END MODULE

	!undefine macros to avoid compiler warnings
#	undef _LFS_type_NAME
#	undef _LFS_CELL_SIZE
#	undef _LFS_EDGE_SIZE
#	undef _LFS_NODE_SIZE

#endif
