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

#		define _BF_type_NAME		samoa_basis_perm
#		define _BF_ORDER			0

#		include "Tools_hierarchical_basis.f90"

		!undefine macros to avoid compiler warnings
#		undef _BF_type_NAME
#		undef _BF_ORDER

	END MODULE

	MODULE Samoa_darcy_p_space
		use SFC_data_types

#		define _BF_type_NAME		samoa_basis_p
#		define _BF_ORDER			_DARCY_P_ORDER

#		include "Tools_hierarchical_basis.f90"

		!undefine macros to avoid compiler warnings
#		undef _BF_type_NAME
#		undef _BF_ORDER

	END MODULE

	MODULE Samoa_darcy_u_space
		use SFC_data_types

#		define _BF_type_NAME		samoa_basis_u
#		define _BF_ORDER			_DARCY_U_ORDER

#		include "Tools_hierarchical_basis.f90"

		!undefine macros to avoid compiler warnings
#		undef _BF_type_NAME
#		undef _BF_ORDER

	END MODULE


	MODULE Samoa_darcy_flow_space
		use SFC_data_types

#		define _BF_type_NAME		samoa_basis_flow
#		define _BF_ORDER			_DARCY_FLOW_ORDER

#		include "Tools_lagrange_basis.f90"

		!undefine macros to avoid compiler warnings
#		undef _BF_type_NAME
#		undef _BF_ORDER

	END MODULE

	!****************
	!Quadrature rules
	!****************

	MODULE Samoa_darcy_p_quadrature_rule
		use SFC_data_types
		use Samoa_darcy_p_space

#		define _BF_type_NAME		samoa_basis_p
#		define _BF_SIZE				_DARCY_P_SIZE
#		define _QR_type_NAME		samoa_qr_p

#		include "Tools_quadrature_rule.f90"

		!undefine macros to avoid compiler warnings
#		undef _BF_type_NAME
#		undef _BF_SIZE
#		undef _QR_type_NAME

	END MODULE

	MODULE Samoa_darcy_u_quadrature_rule
		use SFC_data_types
		use Samoa_darcy_u_space


#		define _BF_type_NAME		samoa_basis_u
#		define _BF_SIZE				_DARCY_U_SIZE
#		define _QR_type_NAME		samoa_qr_u

#		include "Tools_quadrature_rule.f90"

		!undefine macros to avoid compiler warnings
#		undef _BF_type_NAME
#		undef _BF_SIZE
#		undef _QR_type_NAME

	END MODULE

	MODULE Samoa_darcy
		use Tools_log

        use Darcy_gm_A_mod
		use Darcy_gv_p_mod
		use Darcy_gv_r_mod
		use Darcy_gv_d_mod
		use Darcy_gv_A_d_mod
		use Darcy_gv_mat_diagonal_mod
		use Darcy_gv_is_dirichlet_boundary_mod
		use Darcy_gv_u_mod
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
