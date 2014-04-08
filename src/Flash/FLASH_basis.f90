#include "Compilation_control.f90"

#if defined(_FLASH)
#	include "dunavant.f90"

	!*****
	!Bases
	!*****

	MODULE FLASH_Q_space
		use SFC_data_types

#		define _BF_TYPE_NAME		t_basis_Q
#		define _BF_ORDER			_FLASH_ORDER

#		include "Tools_lagrange_basis.f90"
	END MODULE

	MODULE FLASH_flux_space
		use SFC_data_types

#		define _BF_TYPE_NAME		t_basis_flux
#		define _BF_ORDER			_FLASH_ORDER

#		include "Tools_boundary_basis.f90"
	END MODULE

	!****************
	!Quadrature rules
	!****************

	MODULE FLASH_Q_quadrature_rule
		use SFC_data_types
		use FLASH_Q_space

#		define _BF_TYPE_NAME		t_basis_Q
#		define _BF_SIZE				_FLASH_CELL_SIZE
#		define _QR_TYPE_NAME		t_qr_Q

#		include "Tools_quadrature_rule.f90"
	END MODULE

	MODULE Samoa_FLASH
		use FLASH_lfs_flux
		use FLASH_gv_Q

		use FLASH_Q_space
		use FLASH_Q_quadrature_rule
		use FLASH_flux_space

		use Samoa

		type(t_qr_Q)			:: qr_Q

		PUBLIC
	END MODULE
#endif
