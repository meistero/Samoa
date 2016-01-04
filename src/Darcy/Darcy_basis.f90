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
		use Darcy_gv_is_dirichlet_mod
		use Darcy_gv_boundary_condition_mod
		use Darcy_gv_saturation_mod
		use Darcy_gv_flux_mod
		use Darcy_gv_volume_mod

		use Samoa_darcy_p_space
		use Samoa_darcy_u_space
		use Samoa_darcy_flow_space
		use Samoa_darcy_perm_space

		use Samoa
		use Config

		public

		contains

#       if (_DARCY_LAYERS > 0)
            pure subroutine compute_base_fluxes_3D(p, base_permeability, dx, dy, dz, Ax, Ay, Az, g_local, u_w, u_n)
                real (kind = GRID_SR), intent(in)       :: p(:, :), base_permeability(:, :), dx, dy, dz, Ax, Ay, Az, g_local(:)
                real (kind = GRID_SR), intent(out)      :: u_w(:, :), u_n(:, :)

                call compute_base_flux_1D(dx, 0.5_SR * Ax, base_permeability(:, 1), p(1 : _DARCY_LAYERS, 2), p(1 : _DARCY_LAYERS, 1), u_w(:, 1), u_n(:, 1), g_local(1))
                call compute_base_flux_1D(dy, 0.5_SR * Ay, base_permeability(:, 1), p(1 : _DARCY_LAYERS, 2), p(1 : _DARCY_LAYERS, 3), u_w(:, 2), u_n(:, 2), g_local(2))

                call compute_base_flux_1D(dz, 0.25_SR * Az, base_permeability(:, 2), p(1 : _DARCY_LAYERS, 1), p(2 : _DARCY_LAYERS + 1, 1), u_w(:, 3), u_n(:, 3), g_local(3))
                call compute_base_flux_1D(dz, 0.50_SR * Az, base_permeability(:, 2), p(1 : _DARCY_LAYERS, 2), p(2 : _DARCY_LAYERS + 1, 2), u_w(:, 4), u_n(:, 4), g_local(3))
                call compute_base_flux_1D(dz, 0.25_SR * Az, base_permeability(:, 2), p(1 : _DARCY_LAYERS, 3), p(2 : _DARCY_LAYERS + 1, 3), u_w(:, 5), u_n(:, 5), g_local(3))

                call compute_base_flux_1D(dx, 0.5_SR * Ax, base_permeability(:, 1), p(2 : _DARCY_LAYERS + 1, 2), p(2 : _DARCY_LAYERS + 1, 1), u_w(:, 6), u_n(:, 6), g_local(1))
                call compute_base_flux_1D(dy, 0.5_SR * Ay, base_permeability(:, 1), p(2 : _DARCY_LAYERS + 1, 2), p(2 : _DARCY_LAYERS + 1, 3), u_w(:, 7), u_n(:, 7), g_local(2))
            end subroutine

            pure subroutine setup_lse_3D(saturation, p, base_permeability, dx, dy, dz, Ax, Ay, Az, g_local, rhs, matrix)
                real (kind = GRID_SR), intent(in)       :: saturation(:, :), p(:, :), base_permeability(:, :), dx, dy, dz, Ax, Ay, Az, g_local(:)
                real (kind = GRID_SR), intent(inout)    :: rhs(:, :), matrix(:, :)

                real (kind = GRID_SR)   :: lambda_w(_DARCY_LAYERS + 1, 3), lambda_n(_DARCY_LAYERS + 1, 3)
                integer :: i

                lambda_w = l_w(saturation)
                lambda_n = l_n(saturation)

                call setup_lse_1D(dx, 0.5_SR * Ax, base_permeability(:, 1), p(1 : _DARCY_LAYERS, 2), p(1 : _DARCY_LAYERS, 1), lambda_w(1 : _DARCY_LAYERS, 2), lambda_w(1 : _DARCY_LAYERS, 1), lambda_n(1 : _DARCY_LAYERS, 2), lambda_n(1 : _DARCY_LAYERS, 1), g_local(1), rhs(1 : _DARCY_LAYERS, 2), rhs(1 : _DARCY_LAYERS, 1), matrix(:, 1))
                call setup_lse_1D(dy, 0.5_SR * Ay, base_permeability(:, 1), p(1 : _DARCY_LAYERS, 2), p(1 : _DARCY_LAYERS, 3), lambda_w(1 : _DARCY_LAYERS, 2), lambda_w(1 : _DARCY_LAYERS, 3), lambda_n(1 : _DARCY_LAYERS, 2), lambda_n(1 : _DARCY_LAYERS, 3), g_local(2), rhs(1 : _DARCY_LAYERS, 2), rhs(1 : _DARCY_LAYERS, 3), matrix(:, 2))

                !this must be done in a loop due to aliasing of rhs:

                do i = 1, _DARCY_LAYERS
                    call setup_lse_1D(dz, 0.25_SR * Az, base_permeability(i, 2), p(i, 1), p(i + 1, 1), lambda_w(i, 1), lambda_w(i + 1, 1), lambda_n(i, 1), lambda_n(i + 1, 1), g_local(3), rhs(i, 1), rhs(i + 1, 1), matrix(i, 3))
                    call setup_lse_1D(dz, 0.50_SR * Az, base_permeability(i, 2), p(i, 2), p(i + 1, 2), lambda_w(i, 2), lambda_w(i + 1, 2), lambda_n(i, 2), lambda_n(i + 1, 2), g_local(3), rhs(i, 2), rhs(i + 1, 2), matrix(i, 4))
                    call setup_lse_1D(dz, 0.25_SR * Az, base_permeability(i, 2), p(i, 3), p(i + 1, 3), lambda_w(i, 3), lambda_w(i + 1, 3), lambda_n(i, 3), lambda_n(i + 1, 3), g_local(3), rhs(i, 3), rhs(i + 1, 3), matrix(i, 5))
                end do

                call setup_lse_1D(dx, 0.5_SR * Ax, base_permeability(:, 1), p(2 : _DARCY_LAYERS + 1, 2), p(2 : _DARCY_LAYERS + 1, 1), lambda_w(2 : _DARCY_LAYERS + 1, 2), lambda_w(2 : _DARCY_LAYERS + 1, 1), lambda_n(2 : _DARCY_LAYERS + 1, 2), lambda_n(2 : _DARCY_LAYERS + 1, 1), g_local(1), rhs(2 : _DARCY_LAYERS + 1, 2), rhs(2 : _DARCY_LAYERS + 1, 1), matrix(:, 6))
                call setup_lse_1D(dy, 0.5_SR * Ay, base_permeability(:, 1), p(2 : _DARCY_LAYERS + 1, 2), p(2 : _DARCY_LAYERS + 1, 3), lambda_w(2 : _DARCY_LAYERS + 1, 2), lambda_w(2 : _DARCY_LAYERS + 1, 3), lambda_n(2 : _DARCY_LAYERS + 1, 2), lambda_n(2 : _DARCY_LAYERS + 1, 3), g_local(2), rhs(2 : _DARCY_LAYERS + 1, 2), rhs(2 : _DARCY_LAYERS + 1, 3), matrix(:, 7))
            end subroutine

            pure subroutine compute_wave_speeds_3D(saturation, u_w, u_n, xi_w, xi_n)
                real (kind = GRID_SR), intent(in)       :: saturation(:, :), u_w(:, :), u_n(:, :)
                real (kind = GRID_SR), intent(inout)    :: xi_w(:, :), xi_n(:, :)

                real (kind = GRID_SR)   :: lambda_w(_DARCY_LAYERS + 1, 3), lambda_n(_DARCY_LAYERS + 1, 3)
                integer :: i

                lambda_w = dl_w_dS(saturation)
                lambda_n = dl_n_dS(saturation)

                call compute_wave_speed_1D(lambda_w(1 : _DARCY_LAYERS, 2), lambda_w(1 : _DARCY_LAYERS, 1), lambda_n(1 : _DARCY_LAYERS, 2), lambda_n(1 : _DARCY_LAYERS, 1), u_w(:, 1), u_n(:, 1), xi_w(1 : _DARCY_LAYERS, 2), xi_w(1 : _DARCY_LAYERS, 1), xi_n(1 : _DARCY_LAYERS, 2), xi_n(1 : _DARCY_LAYERS, 1))
                call compute_wave_speed_1D(lambda_w(1 : _DARCY_LAYERS, 2), lambda_w(1 : _DARCY_LAYERS, 3), lambda_n(1 : _DARCY_LAYERS, 2), lambda_n(1 : _DARCY_LAYERS, 3), u_w(:, 2), u_n(:, 2), xi_w(1 : _DARCY_LAYERS, 2), xi_w(1 : _DARCY_LAYERS, 3), xi_n(1 : _DARCY_LAYERS, 2), xi_n(1 : _DARCY_LAYERS, 3))

                !this must be done in a loop due to aliasing of xi_w and of xi_n:

                do i = 1, _DARCY_LAYERS
                    call compute_wave_speed_1D(lambda_w(i, 1), lambda_w(i + 1, 1), lambda_n(i, 1), lambda_n(i + 1, 1), u_w(i, 3), u_n(i, 3), xi_w(i, 1), xi_w(i + 1, 1), xi_n(i, 1), xi_n(i + 1, 1))
                    call compute_wave_speed_1D(lambda_w(i, 2), lambda_w(i + 1, 2), lambda_n(i, 2), lambda_n(i + 1, 2), u_w(i, 4), u_n(i, 4), xi_w(i, 2), xi_w(i + 1, 2), xi_n(i, 2), xi_n(i + 1, 2))
                    call compute_wave_speed_1D(lambda_w(i, 3), lambda_w(i + 1, 3), lambda_n(i, 3), lambda_n(i + 1, 3), u_w(i, 5), u_n(i, 5), xi_w(i, 3), xi_w(i + 1, 3), xi_n(i, 3), xi_n(i + 1, 3))
                end do

                call compute_wave_speed_1D(lambda_w(2 : _DARCY_LAYERS + 1, 2), lambda_w(2 : _DARCY_LAYERS + 1, 1), lambda_n(2 : _DARCY_LAYERS + 1, 2), lambda_n(2 : _DARCY_LAYERS + 1, 1), u_w(:, 6), u_n(:, 6), xi_w(2 : _DARCY_LAYERS + 1, 2), xi_w(2 : _DARCY_LAYERS + 1, 1), xi_n(2 : _DARCY_LAYERS + 1, 2), xi_n(2 : _DARCY_LAYERS + 1, 1))
                call compute_wave_speed_1D(lambda_w(2 : _DARCY_LAYERS + 1, 2), lambda_w(2 : _DARCY_LAYERS + 1, 3), lambda_n(2 : _DARCY_LAYERS + 1, 2), lambda_n(2 : _DARCY_LAYERS + 1, 3), u_w(:, 7), u_n(:, 7), xi_w(2 : _DARCY_LAYERS + 1, 2), xi_w(2 : _DARCY_LAYERS + 1, 3), xi_n(2 : _DARCY_LAYERS + 1, 2), xi_n(2 : _DARCY_LAYERS + 1, 3))
           end subroutine


            pure subroutine compute_flux_vector_3D(saturation, u_w, u_n, flux_w, flux_n)
                real (kind = GRID_SR), intent(in)          :: saturation(:, :), u_w(:, :), u_n(:, :)
                real (kind = GRID_SR), intent(inout)	   :: flux_w(:, :), flux_n(:, :)

                real (kind = GRID_SR)   :: dummy(1 : _DARCY_LAYERS), lambda_w(_DARCY_LAYERS + 1, 3), lambda_n(_DARCY_LAYERS + 1, 3)

                integer :: i

                lambda_w = l_w(saturation)
                lambda_n = l_n(saturation)

                call compute_upwind_flux_1D(lambda_w(1 : _DARCY_LAYERS, 2), lambda_w(1 : _DARCY_LAYERS, 1), lambda_n(1 : _DARCY_LAYERS, 2), lambda_n(1 : _DARCY_LAYERS, 1), u_w(:, 1), u_n(:, 1), flux_w(1 : _DARCY_LAYERS, 1), dummy, flux_n(1 : _DARCY_LAYERS, 1), dummy)
                call compute_upwind_flux_1D(lambda_w(1 : _DARCY_LAYERS, 2), lambda_w(1 : _DARCY_LAYERS, 3), lambda_n(1 : _DARCY_LAYERS, 2), lambda_n(1 : _DARCY_LAYERS, 3), u_w(:, 2), u_n(:, 2), flux_w(1 : _DARCY_LAYERS, 2), dummy, flux_n(1 : _DARCY_LAYERS, 2), dummy)

                call compute_upwind_flux_1D(lambda_w(1 : _DARCY_LAYERS, 1), lambda_w(2 : _DARCY_LAYERS + 1, 1), lambda_n(1 : _DARCY_LAYERS, 1), lambda_n(2 : _DARCY_LAYERS + 1, 1), u_w(:, 3), u_n(:, 3), flux_w(1 : _DARCY_LAYERS, 3), dummy, flux_n(1 : _DARCY_LAYERS, 3), dummy)
                call compute_upwind_flux_1D(lambda_w(1 : _DARCY_LAYERS, 2), lambda_w(2 : _DARCY_LAYERS + 1, 2), lambda_n(1 : _DARCY_LAYERS, 2), lambda_n(2 : _DARCY_LAYERS + 1, 2), u_w(:, 4), u_n(:, 4), flux_w(1 : _DARCY_LAYERS, 3), dummy, flux_n(1 : _DARCY_LAYERS, 3), dummy)
                call compute_upwind_flux_1D(lambda_w(1 : _DARCY_LAYERS, 3), lambda_w(2 : _DARCY_LAYERS + 1, 3), lambda_n(1 : _DARCY_LAYERS, 3), lambda_n(2 : _DARCY_LAYERS + 1, 3), u_w(:, 5), u_n(:, 5), flux_w(1 : _DARCY_LAYERS, 3), dummy, flux_n(1 : _DARCY_LAYERS, 3), dummy)

                call compute_upwind_flux_1D(lambda_w(2 : _DARCY_LAYERS + 1, 2), lambda_w(2 : _DARCY_LAYERS + 1, 1), lambda_n(2 : _DARCY_LAYERS + 1, 2), lambda_n(2 : _DARCY_LAYERS + 1, 1), u_w(:, 6), u_n(:, 6), flux_w(1 : _DARCY_LAYERS, 1), dummy, flux_n(1 : _DARCY_LAYERS, 1), dummy)
                call compute_upwind_flux_1D(lambda_w(2 : _DARCY_LAYERS + 1, 2), lambda_w(2 : _DARCY_LAYERS + 1, 3), lambda_n(2 : _DARCY_LAYERS + 1, 2), lambda_n(2 : _DARCY_LAYERS + 1, 3), u_w(:, 7), u_n(:, 7), flux_w(1 : _DARCY_LAYERS, 2), dummy, flux_n(1 : _DARCY_LAYERS, 2), dummy)
            end subroutine

            pure subroutine compute_fluxes_3D(saturation, u_w, u_n, flux_w, flux_n)
                real (kind = GRID_SR), intent(in)          :: saturation(:, :), u_w(:, :), u_n(:, :)
                real (kind = GRID_SR), intent(inout)	   :: flux_w(:, :), flux_n(:, :)

                real (kind = GRID_SR)   :: lambda_w(_DARCY_LAYERS + 1, 3), lambda_n(_DARCY_LAYERS + 1, 3)

                integer :: i

                lambda_w = l_w(saturation)
                lambda_n = l_n(saturation)

                call compute_upwind_flux_1D(lambda_w(1 : _DARCY_LAYERS, 2), lambda_w(1 : _DARCY_LAYERS, 1), lambda_n(1 : _DARCY_LAYERS, 2), lambda_n(1 : _DARCY_LAYERS, 1), u_w(:, 1), u_n(:, 1), flux_w(1 : _DARCY_LAYERS, 2), flux_w(1 : _DARCY_LAYERS, 1), flux_n(1 : _DARCY_LAYERS, 2), flux_n(1 : _DARCY_LAYERS, 1))
                call compute_upwind_flux_1D(lambda_w(1 : _DARCY_LAYERS, 2), lambda_w(1 : _DARCY_LAYERS, 3), lambda_n(1 : _DARCY_LAYERS, 2), lambda_n(1 : _DARCY_LAYERS, 3), u_w(:, 2), u_n(:, 2), flux_w(1 : _DARCY_LAYERS, 2), flux_w(1 : _DARCY_LAYERS, 3), flux_n(1 : _DARCY_LAYERS, 2), flux_n(1 : _DARCY_LAYERS, 3))

                !this must be done in a loop due to aliasing of flux_w and of flux_n:

                do i = 1, _DARCY_LAYERS
                    call compute_upwind_flux_1D(lambda_w(i, 1), lambda_w(i + 1, 1), lambda_n(i, 1), lambda_n(i + 1, 1), u_w(i, 3), u_n(i, 3), flux_w(i, 1), flux_w(i + 1, 1), flux_n(i, 1), flux_n(i + 1, 1))
                    call compute_upwind_flux_1D(lambda_w(i, 2), lambda_w(i + 1, 2), lambda_n(i, 2), lambda_n(i + 1, 2), u_w(i, 4), u_n(i, 4), flux_w(i, 2), flux_w(i + 1, 2), flux_n(i, 2), flux_n(i + 1, 2))
                    call compute_upwind_flux_1D(lambda_w(i, 3), lambda_w(i + 1, 3), lambda_n(i, 3), lambda_n(i + 1, 3), u_w(i, 5), u_n(i, 5), flux_w(i, 3), flux_w(i + 1, 3), flux_n(i, 3), flux_n(i + 1, 3))
                end do

                call compute_upwind_flux_1D(lambda_w(2 : _DARCY_LAYERS + 1, 2), lambda_w(2 : _DARCY_LAYERS + 1, 1), lambda_n(2 : _DARCY_LAYERS + 1, 2), lambda_n(2 : _DARCY_LAYERS + 1, 1), u_w(:, 6), u_n(:, 6), flux_w(2 : _DARCY_LAYERS + 1, 2), flux_w(2 : _DARCY_LAYERS + 1, 1), flux_n(2 : _DARCY_LAYERS + 1, 2), flux_n(2 : _DARCY_LAYERS + 1, 1))
                call compute_upwind_flux_1D(lambda_w(2 : _DARCY_LAYERS + 1, 2), lambda_w(2 : _DARCY_LAYERS + 1, 3), lambda_n(2 : _DARCY_LAYERS + 1, 2), lambda_n(2 : _DARCY_LAYERS + 1, 3), u_w(:, 7), u_n(:, 7), flux_w(2 : _DARCY_LAYERS + 1, 2), flux_w(2 : _DARCY_LAYERS + 1, 3), flux_n(2 : _DARCY_LAYERS + 1, 2), flux_n(2 : _DARCY_LAYERS + 1, 3))
            end subroutine
#       else
            pure subroutine compute_base_fluxes_2D(p, base_permeability, dx, dy, Ax, Ay, g_local, u_w, u_n)
                real (kind = GRID_SR), intent(in)       :: p(:), base_permeability, dx, dy, Ax, Ay, g_local(:)
                real (kind = GRID_SR), intent(out)      :: u_w(:), u_n(:)

                call compute_base_flux_1D(dx, Ax, base_permeability, p(2), p(1), u_w(1), u_n(1), g_local(1))
                call compute_base_flux_1D(dy, Ay, base_permeability, p(2), p(3), u_w(2), u_n(2), g_local(2))
            end subroutine

            pure subroutine setup_lse_2D(saturation, p, base_permeability, dx, dy, Ax, Ay, g_local, rhs, matrix)
                real (kind = GRID_SR), intent(in)       :: saturation(:), p(:), base_permeability, dx, dy, Ax, Ay, g_local(:)
                real (kind = GRID_SR), intent(inout)    :: rhs(:), matrix(:)

                real (kind = GRID_SR)   :: lambda_w(3), lambda_n(3)

                lambda_w = l_w(saturation)
                lambda_n = l_n(saturation)

                call setup_lse_1D(dx, Ax, base_permeability, p(2), p(1), lambda_w(2), lambda_w(1), lambda_n(2), lambda_n(1), g_local(1), rhs(2), rhs(1), matrix(1))
                call setup_lse_1D(dy, Ay, base_permeability, p(2), p(3), lambda_w(2), lambda_w(3), lambda_n(2), lambda_n(3), g_local(2), rhs(2), rhs(3), matrix(2))
            end subroutine

            pure subroutine compute_wave_speeds_2D(saturation, u_w, u_n, xi_w, xi_n)
                real (kind = GRID_SR), intent(in)       :: saturation(:), u_w(:), u_n(:)
                real (kind = GRID_SR), intent(inout)    :: xi_w(:), xi_n(:)

                real (kind = GRID_SR)   :: lambda_w(3), lambda_n(3)

                lambda_w = dl_w_dS(saturation)
                lambda_n = dl_n_dS(saturation)

                call compute_wave_speed_1D(lambda_w(2), lambda_w(1), lambda_n(2), lambda_n(1), u_w(1), u_n(1), xi_w(2), xi_w(1), xi_n(2), xi_n(1))
                call compute_wave_speed_1D(lambda_w(2), lambda_w(3), lambda_n(2), lambda_n(3), u_w(2), u_n(2), xi_w(2), xi_w(3), xi_n(2), xi_n(3))
           end subroutine

            pure subroutine compute_flux_vector_2D(saturation, u_w, u_n, flux_w, flux_n)
                real (kind = GRID_SR), intent(in)          :: saturation(:), u_w(:), u_n(:)
                real (kind = GRID_SR), intent(inout)	   :: flux_w(:), flux_n(:)

                real (kind = GRID_SR)   :: dummy, lambda_w(3), lambda_n(3)

                lambda_w = l_w(saturation)
                lambda_n = l_n(saturation)

                call compute_upwind_flux_1D(lambda_w(2), lambda_w(1), lambda_n(2), lambda_n(1), u_w(1), u_n(1), flux_w(1), dummy, flux_n(1), dummy)
                call compute_upwind_flux_1D(lambda_w(2), lambda_w(3), lambda_n(2), lambda_n(3), u_w(2), u_n(2), flux_w(2), dummy, flux_n(2), dummy)
            end subroutine

            pure subroutine compute_fluxes_2D(saturation, u_w, u_n, flux_w, flux_n)
                real (kind = GRID_SR), intent(in)          :: saturation(:), u_w(:), u_n(:)
                real (kind = GRID_SR), intent(inout)	   :: flux_w(:), flux_n(:)

                real (kind = GRID_SR)   :: lambda_w(3), lambda_n(3)

                lambda_w = l_w(saturation)
                lambda_n = l_n(saturation)

                call compute_upwind_flux_1D(lambda_w(2), lambda_w(1), lambda_n(2), lambda_n(1), u_w(1), u_n(1), flux_w(2), flux_w(1), flux_n(2), flux_n(1))
                call compute_upwind_flux_1D(lambda_w(2), lambda_w(3), lambda_n(2), lambda_n(3), u_w(2), u_n(2), flux_w(2), flux_w(3), flux_n(2), flux_n(3))
            end subroutine
#       endif

        elemental subroutine compute_base_flux_1D(dx, area, base_permeability, pL, pR, u_w, u_n, g_local)
            real (kind = GRID_SR), intent(in)       :: dx, area, base_permeability, pL, pR, g_local
            real (kind = GRID_SR), intent(out)      :: u_w, u_n

            u_w = area * base_permeability * (-(pR - pL) / dx + cfg%r_rho_w * g_local)
            u_n = area * base_permeability * (-(pR - pL) / dx + cfg%r_rho_n * g_local)
        end subroutine

        elemental subroutine setup_lse_1D(dx, area, base_permeability, pL, pR, lambda_wL, lambda_wR, lambda_nL, lambda_nR, g_local, rhs_l, rhs_r, mat_l)
            real (kind = GRID_SR), intent(in)       :: area, dx, base_permeability, pL, pR, lambda_wL, lambda_wR, lambda_nL, lambda_nR, g_local
            real (kind = GRID_SR), intent(inout)    :: rhs_l, rhs_r, mat_l

            real (kind = SR)    :: u_w, u_n, lambda_w_local, lambda_n_local, rhs_local

            u_w = area * base_permeability * (-(pR - pL) / dx + cfg%r_rho_w * g_local)
            u_n = area * base_permeability * (-(pR - pL) / dx + cfg%r_rho_n * g_local)

            if (u_w > 0) then
                lambda_w_local = lambda_wL
            else
                lambda_w_local = lambda_wR
            end if

            if (u_n > 0) then
                lambda_n_local = lambda_nL
            else
                lambda_n_local = lambda_nR
            end if

            rhs_local = area * base_permeability * (lambda_w_local * cfg%r_rho_w + lambda_n_local * cfg%r_rho_n) * g_local

            rhs_l = rhs_l - rhs_local
            rhs_r = rhs_r + rhs_local
            mat_l = area / dx * base_permeability * (lambda_w_local + lambda_n_local)
        end subroutine

        elemental subroutine compute_wave_speed_1D(lambda_w_l, lambda_w_r, lambda_n_l, lambda_n_r, u_w, u_n, xi_wl, xi_wr, xi_nl, xi_nr)
            real (kind = GRID_SR), intent(in)       :: lambda_w_l, lambda_w_r, lambda_n_l, lambda_n_r, u_w, u_n
            real (kind = GRID_SR), intent(inout)    :: xi_wl, xi_wr, xi_nl, xi_nr

            real (kind = GRID_SR) :: flux_w, flux_n

            !Find the wavespeeds by precomputing the fluxes

            flux_w = lambda_w_l * max(u_w, 0.0_SR) + lambda_w_r * min(u_w, 0.0_SR)
            flux_n = lambda_n_l * max(u_n, 0.0_SR) + lambda_n_r * min(u_n, 0.0_SR)

            xi_wl = xi_wl + abs(flux_w)
            xi_wr = xi_wr + abs(flux_w)
            xi_nl = xi_nl + abs(flux_n)
            xi_nr = xi_nr + abs(flux_n)
        end subroutine

        elemental subroutine compute_upwind_flux_1D(lambda_w_l, lambda_w_r, lambda_n_l, lambda_n_r, u_w, u_n, flux_w_l, flux_w_r, flux_n_l, flux_n_r)
            real (kind = GRID_SR), intent(in)       :: lambda_w_l, lambda_w_r, lambda_n_l, lambda_n_r, u_w, u_n
            real (kind = GRID_SR), intent(inout)	:: flux_w_l, flux_w_r, flux_n_l, flux_n_r

            real (kind = GRID_SR) :: flux_w, flux_n

            flux_w = lambda_w_l * max(u_w, 0.0_SR) + lambda_w_r * min(u_w, 0.0_SR)
            flux_n = lambda_n_l * max(u_n, 0.0_SR) + lambda_n_r * min(u_n, 0.0_SR)

            flux_w_l = flux_w_l + flux_w
            flux_w_r = flux_w_r - flux_w
            flux_n_l = flux_n_l + flux_n
            flux_n_r = flux_n_r - flux_n
        end subroutine

        !computes the algebraic wetting phase flux from the total flux and the gravity term
        !with the formula u_w = l_w / (l_w + l_n) * (u_T + l_n * (rho_w - rho_n) * k * g)
        elemental function u_w(lambda_w, lambda_n, u_T, rho_w_minus_rho_n_k_g)
            real (kind = GRID_SR), intent(in)   :: lambda_w, lambda_n, u_T, rho_w_minus_rho_n_k_g
            real (kind = GRID_SR)               :: u_w

            u_w = lambda_w / (lambda_w + lambda_n) * (u_T + lambda_n * rho_w_minus_rho_n_k_g)
        end function

        !computes the algebraic non-wetting phase flux from the total flux and the gravity term
        !with the formula u_n = l_n / (l_w + l_n) * (u_T - l_w * (rho_w - rho_n) * k * g)
        elemental function u_n(lambda_w, lambda_n, u_T, rho_w_minus_rho_n_k_g)
            real (kind = GRID_SR), intent(in)   :: lambda_w, lambda_n, u_T, rho_w_minus_rho_n_k_g
            real (kind = GRID_SR)               :: u_n

            u_n = lambda_n / (lambda_w + lambda_n) * (u_T - lambda_w * rho_w_minus_rho_n_k_g)
        end function

        !computes the wetting phase mobility from the saturation using one of three possible models
        elemental function l_w(S)
            real (kind = GRID_SR), intent(in)   :: S
            real (kind = GRID_SR)               :: l_w
            real, parameter                     :: lambda = 2.0_SR

#           if defined (_DARCY_MOB_LINEAR)
                l_w = S / cfg%r_nu_w
#           elif defined (_DARCY_MOB_QUADRATIC)
                l_w = S * S / cfg%r_nu_w
#           elif defined (_DARCY_MOB_BROOKS_COREY)
                l_w = S * S * (S ** (2.0_SR / lambda + 1.0_SR)) / cfg%r_nu_w
#           endif
        end function

        !computes the term l_w(S) / S
        elemental function dl_w_dS(S)
            real (kind = GRID_SR), intent(in)   :: S
            real (kind = GRID_SR)               :: dl_w_dS
            real, parameter                     :: lambda = 2.0_SR

#           if defined (_DARCY_MOB_LINEAR)
                dl_w_dS = 1.0_SR / cfg%r_nu_w
#           elif defined (_DARCY_MOB_QUADRATIC)
                dl_w_dS = S / cfg%r_nu_w
#           elif defined (_DARCY_MOB_BROOKS_COREY)
                !if N := 2 / lambda + 2 is integer, the result will be \sum_{k = 0}^N S_l^k S_r^{N-k}
                !Since we cannot assume that we estimtate
                !if S_l == S_r, the result will be accurate, otherwise this is only an estimate
                dl_w_dS = S * (S ** (2.0_SR / lambda + 1.0_SR)) / cfg%r_nu_w
#           endif
        end function

        !computes the non-wetting phase mobility from the saturation using one of three possible models
        elemental function l_n(S)
            real (kind = GRID_SR), intent(in)   :: S
            real (kind = GRID_SR)               :: l_n
            real, parameter                     :: lambda = 2.0_SR

#           if defined (_DARCY_MOB_LINEAR)
                l_n = (1.0_SR - S) / cfg%r_nu_n
#           elif defined (_DARCY_MOB_QUADRATIC)
                l_n = (1.0_SR - S) * (1.0_SR - S) / cfg%r_nu_n
#           elif defined (_DARCY_MOB_BROOKS_COREY)
                l_n = (1.0_SR - S) * (1.0_SR - S) * (1.0_SR - S ** (2.0_SR / lambda + 1.0_SR)) / cfg%r_nu_n
#           endif
        end function

        !computes the term l_n(S) / (1 - S)
        elemental function dl_n_dS(S)
            real (kind = GRID_SR), intent(in)   :: S
            real (kind = GRID_SR)               :: dl_n_dS
            real, parameter                     :: lambda = 2.0_SR

#           if defined (_DARCY_MOB_LINEAR)
                dl_n_dS = 1.0_SR / cfg%r_nu_n
#           elif defined (_DARCY_MOB_QUADRATIC)
                dl_n_dS = (1.0_SR - S) / cfg%r_nu_n
#           elif defined (_DARCY_MOB_BROOKS_COREY)
                dl_n_dS = (1.0_SR - S) * (1.0_SR - S ** (2.0_SR / lambda + 1.0_SR)) / cfg%r_nu_n
#           endif
        end function
	END MODULE
#endif
