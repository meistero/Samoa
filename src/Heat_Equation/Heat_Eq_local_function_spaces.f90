! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_HEAT_EQ)
	MODULE Heat_Eq_grid_matrix
		use SFC_data_types

		PUBLIC

		real (kind = GRID_SR), DIMENSION(_HEAT_EQ_SIZE, _HEAT_EQ_SIZE)		:: gm_stiffness_matrix, gm_mass_matrix		!< element matrices
	END MODULE Heat_Eq_grid_matrix

	!*****************
	!temperature space
	!*****************

#	define _LFS_type_NAME		heat_eq_lfs_T
#	define _LFS_CELL_SIZE		_HEAT_EQ_CELL_SIZE
#	define _LFS_EDGE_SIZE		_HEAT_EQ_EDGE_SIZE
#	define _LFS_NODE_SIZE		_HEAT_EQ_NODE_SIZE

	MODULE Heat_Eq_lfs_T_mod
		use SFC_data_types

#		define _LFS_type		real (kind = GRID_SR)
#		include "Tools_local_function_space.f90"
	END MODULE

	MODULE Heat_Eq_gv_T_mod
		use SFC_data_types
		use Heat_Eq_lfs_T_mod

#		define _GV_type_NAME		heat_eq_gv_T
#		define _GV_type				real (kind = GRID_SR)
#		define _GV_NAME				T
#		define _GV_PERSISTENT		.true.

#		include "Tools_grid_variable.f90"
	END MODULE

	MODULE Heat_Eq_gv_T_temp_mod
		use SFC_data_types
		use Heat_Eq_lfs_T_mod

#		define _GV_type_NAME		heat_eq_gv_T_temp
#		define _GV_type				real (kind = GRID_SR)
#		define _GV_NAME				T_temp
#		define _GV_PERSISTENT		.true.

#		include "Tools_grid_variable.f90"
	END MODULE

	MODULE Heat_Eq_gv_r_mod
		use SFC_data_types
		use Heat_Eq_lfs_T_mod

#		define _GV_type_NAME		heat_eq_gv_r
#		define _GV_type				real (kind = GRID_SR)
#		define _GV_NAME				r
#		define _GV_PERSISTENT		.false.

#		include "Tools_grid_variable.f90"
	END MODULE

	MODULE Heat_Eq_gv_mat_mass_diagonal_mod
		use SFC_data_types
		use Heat_Eq_lfs_T_mod

#		define _GV_type_NAME		heat_eq_gv_mat_mass_diagonal
#		define _GV_type				real (kind = GRID_SR)
#		define _GV_NAME				mat_mass_diagonal
#		define _GV_PERSISTENT		.false.

#		include "Tools_grid_variable.f90"
	END MODULE

#	undef _LFS_type_NAME
#	undef _LFS_CELL_SIZE
#	undef _LFS_EDGE_SIZE
#	undef _LFS_NODE_SIZE
#endif
