! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_DARCY)
#	include "dunavant.f90"

	!*****
	!Bases
	!*****
	MODULE Samoa_darcy_perm_space
		use SFC_data_types

#		define _BF_TYPE_NAME		samoa_basis_perm
#		define _BF_ORDER			0

#		include "Tools_hierarchical_basis.f90"

		!undefine macros to avoid compiler warnings
#		undef _BF_TYPE_NAME
#		undef _BF_ORDER

	END MODULE

	MODULE Samoa_darcy_p_space
		use SFC_data_types

#		define _BF_TYPE_NAME		samoa_basis_p
#		define _BF_ORDER			1

#		include "Tools_hierarchical_basis.f90"

		!undefine macros to avoid compiler warnings
#		undef _BF_TYPE_NAME
#		undef _BF_ORDER

	END MODULE

	MODULE Samoa_darcy_u_space
		use SFC_data_types

#		define _BF_TYPE_NAME		samoa_basis_u
#		define _BF_ORDER			0

#		include "Tools_hierarchical_basis.f90"

		!undefine macros to avoid compiler warnings
#		undef _BF_TYPE_NAME
#		undef _BF_ORDER

	END MODULE


	MODULE Samoa_darcy_flow_space
		use SFC_data_types

#		define _BF_TYPE_NAME		samoa_basis_flow
#		define _BF_ORDER			1

#		include "Tools_lagrange_basis.f90"

		!undefine macros to avoid compiler warnings
#		undef _BF_TYPE_NAME
#		undef _BF_ORDER

	END MODULE

	MODULE Samoa_darcy
		use Tools_log

        use Darcy_gm_A_mod
		use Darcy_gv_p_mod
		use Darcy_gv_rhs_mod
		use Darcy_gv_r_mod
		use Darcy_gv_d_mod
		use Darcy_gv_A_d_mod
		use Darcy_gv_mat_diagonal_mod
		use Darcy_gv_is_pressure_dirichlet_boundary_mod
		use Darcy_gv_is_saturation_dirichlet_boundary_mod
		use Darcy_gv_saturation_mod
		use Darcy_gv_flux_mod
		use Darcy_gv_volume_mod

		use Samoa_darcy_p_space
		use Samoa_darcy_u_space
		use Samoa_darcy_flow_space
		use Samoa_darcy_perm_space

		use Samoa

		public
	END MODULE
#endif
