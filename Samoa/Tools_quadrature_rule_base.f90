! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

MODULE Tools_quadrature_rule_base
	use SFC_data_types

	type t_qr_base
		integer (kind = GRID_SI)									:: i_qpts						!< Number of quadrature points
		real (kind = GRID_SR), DIMENSION(:, :), ALLOCATABLE			:: r_qpt_coords					!< Positions of the quadrature points (barycentric)
		real (kind = GRID_SR), DIMENSION(:, :), ALLOCATABLE			:: r_qpt_normals				!< Normal of the quadrature point if on a boundary (barycentric)

		real (kind = GRID_SR), DIMENSION(:), ALLOCATABLE			:: r_qpt_weights				!< Weights of the quadrature points

		real (kind = GRID_SR), DIMENSION(:, :), ALLOCATABLE			:: r_dof_qpt_psi				!< Form function value for a degree of freedom at a quadrature point
		real (kind = GRID_SR), DIMENSION(:, :), ALLOCATABLE			:: r_dof_qpt_d_psi_d_lambda1	!< Partial derivative of the form function by lambda1 (barycentric) for a degree of freedom at a quadrature point
		real (kind = GRID_SR), DIMENSION(:, :), ALLOCATABLE			:: r_dof_qpt_d_psi_d_lambda2	!< Partial derivative of the form function by lambda2 (barycentric) for a degree of freedom at a quadrature point
	end type

	PUBLIC

	CONTAINS

	!**********************
	!Constructor/Destructor
	!**********************

	subroutine t_qr_base_destroy(qr)
		implicit none
		type(t_qr_base), intent(inout)				:: qr			!< quadrature formula

		integer (kind = GRID_SI)					:: i_error

		deallocate(qr%r_qpt_coords, stat = i_error); assert_eq(i_error, 0)
		deallocate(qr%r_qpt_normals, stat = i_error); assert_eq(i_error, 0)
		deallocate(qr%r_qpt_weights, stat = i_error); assert_eq(i_error, 0)
		deallocate(qr%r_dof_qpt_psi, stat = i_error); assert_eq(i_error, 0)
		deallocate(qr%r_dof_qpt_d_psi_d_lambda1, stat = i_error); assert_eq(i_error, 0)
		deallocate(qr%r_dof_qpt_d_psi_d_lambda2, stat = i_error); assert_eq(i_error, 0)
	end subroutine

	!*****************************
	!Per-Qpt transfomrations
	!*****************************

	!> Returns the derivation of the form function for a degree of freedom at a quadrature point (barycentric coordinates)
	pure function t_qr_base_d_psi_d_lambda(qr, i_dof, i_qpt) result(d_psi_d_lambda)
		implicit none
		type(t_qr_base), intent(in)						:: qr
		integer (kind = GRID_SI), intent(in)			:: i_dof, i_qpt
		real (kind = GRID_SR), DIMENSION(2)				:: d_psi_d_lambda

		d_psi_d_lambda(:) =	(/ qr%r_dof_qpt_d_psi_d_lambda1(i_dof, i_qpt), qr%r_dof_qpt_d_psi_d_lambda2(i_dof, i_qpt) /)
	end function
END MODULE
