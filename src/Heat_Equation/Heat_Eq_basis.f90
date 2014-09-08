! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_HEAT_EQ)
	!*****
	!Bases
	!*****

	MODULE Samoa_heat_eq_T_space
		use SFC_data_types

#		define _BF_TYPE_NAME		samoa_basis_T
#		define _BF_ORDER			_HEAT_EQ_ORDER

#		include "Tools_lagrange_basis.f90"

		!undefine macros to avoid compiler warnings
#		undef _BF_TYPE_NAME
#		undef _BF_ORDER

	END MODULE

	!****************
	!Quadrature rules
	!****************

	MODULE Samoa_heat_eq_T_quadrature_rule
		use SFC_data_types
		use Samoa_heat_eq_T_space

#		define _BF_TYPE_NAME		samoa_basis_T
#		define _BF_SIZE				_HEAT_EQ_SIZE
#		define _QR_TYPE_NAME		samoa_qr_T

#		include "Tools_quadrature_rule.f90"

		!undefine macros to avoid compiler warnings
#		undef _BF_TYPE_NAME
#		undef _BF_SIZE
#		undef _QR_TYPE_NAME

	END MODULE

	MODULE Samoa_heat_eq
		use Heat_Eq_grid_matrix

		use Heat_Eq_lfs_T_mod
		use Heat_Eq_gv_T_mod
		use Heat_Eq_gv_T_temp_mod
		use Heat_Eq_gv_r_mod
		use Heat_Eq_gv_mat_mass_diagonal_mod

		use Samoa_heat_eq_T_space
		use Samoa_heat_eq_T_quadrature_rule

		use Samoa

		PUBLIC
	END MODULE
#endif
