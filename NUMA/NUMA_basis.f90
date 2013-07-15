! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_NUMA)
#	include "dunavant.f90"

	!*****
	!Bases
	!*****

	MODULE NUMA_Q_space
		use SFC_data_types

#		define _BF_type_NAME		t_basis_Q
#		define _BF_ORDER			_NUMA_ORDER

#		include "Tools_lagrange_basis.f90"
	END MODULE

	MODULE NUMA_flux_space
		use SFC_data_types

#		define _BF_type_NAME		t_basis_flux
#		define _BF_ORDER			_NUMA_ORDER

#		include "Tools_boundary_basis.f90"
	END MODULE

	!****************
	!Quadrature rules
	!****************

	MODULE NUMA_Q_quadrature_rule
		use SFC_data_types
		use NUMA_Q_space

#		define _BF_type_NAME		t_basis_Q
#		define _BF_SIZE				_NUMA_CELL_SIZE
#		define _QR_type_NAME		t_qr_Q

#		include "Tools_quadrature_rule.f90"
	END MODULE

	MODULE Samoa_NUMA
		use NUMA_lfs_Q
		use NUMA_lfs_flux
		use NUMA_gv_Q

		use NUMA_Q_space
		use NUMA_Q_quadrature_rule
		use NUMA_flux_space

		use Samoa

		type(t_qr_Q)			:: qr_Q

		PUBLIC
	END MODULE
#endif
