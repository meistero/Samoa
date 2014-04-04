! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_TESTS)
	MODULE Tests_space_l0
		use SFC_data_types

#		define _BF_TYPE_NAME		Tests_basis_l0
#		define _BF_ORDER			0

#		include "Tools_lagrange_basis.f90"
	END MODULE

	MODULE Tests_space_l1
		use SFC_data_types

#		define _BF_TYPE_NAME		Tests_basis_l1
#		define _BF_ORDER			1

#		include "Tools_lagrange_basis.f90"
	END MODULE

	MODULE Tests_space_l2
		use SFC_data_types

#		define _BF_TYPE_NAME		Tests_basis_l2
#		define _BF_ORDER			2

#		include "Tools_lagrange_basis.f90"
	END MODULE

	MODULE Tests_space_h0
		use SFC_data_types

#		define _BF_TYPE_NAME		Tests_basis_h0
#		define _BF_ORDER			0

#		include "Tools_hierarchical_basis.f90"
	END MODULE

	MODULE Tests_space_h1
		use SFC_data_types

#		define _BF_TYPE_NAME		Tests_basis_h1
#		define _BF_ORDER			1

#		include "Tools_hierarchical_basis.f90"
	END MODULE

	MODULE Tests_space_h2
		use SFC_data_types

#		define _BF_TYPE_NAME		Tests_basis_h2
#		define _BF_ORDER			2

#		include "Tools_hierarchical_basis.f90"
	END MODULE

	MODULE Tests_space_b0
		use SFC_data_types

#		define _BF_TYPE_NAME		Tests_basis_b0
#		define _BF_ORDER			0

#		include "Tools_boundary_basis.f90"
	END MODULE

	MODULE Tests_space_b1
		use SFC_data_types

#		define _BF_TYPE_NAME		Tests_basis_b1
#		define _BF_ORDER			1

#		include "Tools_boundary_basis.f90"
	END MODULE

	MODULE Tests_space_b2
		use SFC_data_types

#		define _BF_TYPE_NAME		Tests_basis_b2
#		define _BF_ORDER			2

#		include "Tools_boundary_basis.f90"
	END MODULE

	MODULE Tests_basis_functions_mod
		use SFC_data_types

		use Tests_space_l0
		use Tests_space_l1
		use Tests_space_l2

		use Tests_space_h0
		use Tests_space_h1
		use Tests_space_h2

		use Tests_space_b0
		use Tests_space_b1
		use Tests_space_b2

		contains

		subroutine tests_basis_functions()
			implicit none

			_log_write(0, *) "Tests: running basis function tests.."

			call Tests_basis_l0_test()
			call Tests_basis_l1_test()
			call Tests_basis_l2_test()

			call Tests_basis_h0_test()
			call Tests_basis_h1_test()
			call Tests_basis_h2_test()

			call Tests_basis_b0_test()
			call Tests_basis_b1_test()
			call Tests_basis_b2_test()
		end subroutine

	END MODULE
#endif
