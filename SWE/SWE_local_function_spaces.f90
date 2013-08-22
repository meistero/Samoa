! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_SWE)
	!*****************
	!temperature space
	!*****************

#	define _LFS_type_NAME			t_lfs_Q
#	define _LFS_CELL_SIZE			_SWE_CELL_SIZE
#	define _LFS_EDGE_SIZE			0
#	define _LFS_NODE_SIZE			0

	MODULE SWE_lfs_Q
		use SFC_data_types

#		define _LFS_type			real(kind = GRID_SR)
#		include "Tools_local_function_space.f90"
	END MODULE

	MODULE SWE_gv_Q
		use SFC_data_types
		use SWE_lfs_Q

#		define _GV_type_NAME		t_gv_Q
#		define _GV_type				type(t_state)
#		define _GV_NAME				Q
#		define _GV_PERSISTENT		.true.

#		include "Tools_grid_variable.f90"
	END MODULE

#	undef _LFS_type_NAME
#	undef _LFS_CELL_SIZE
#	undef _LFS_EDGE_SIZE
#	undef _LFS_NODE_SIZE

#	define _LFS_type_NAME			t_lfs_flux
#	define _LFS_CELL_SIZE			0
#	define _LFS_EDGE_SIZE			_SWE_EDGE_SIZE
#	define _LFS_NODE_SIZE			0

	MODULE SWE_lfs_flux
		use SFC_data_types

#		define _LFS_type			real(kind = GRID_SR)
#		include "Tools_local_function_space.f90"
	END MODULE

#	undef _LFS_type_NAME
#	undef _LFS_CELL_SIZE
#	undef _LFS_EDGE_SIZE
#	undef _LFS_NODE_SIZE
#endif
