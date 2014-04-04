! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


!> Generic quadrature rule for a template basis
!> Warning: this a template module body and requires preprocessor commands before inclusion.
!> Usage: Inside your module definition, insert
!>
!> #define _BF_TYPE_NAME	<value>
!> #define _BF_SIZE			<value>
!> #define _QR_TYPE_NAME	<value>
!>
!> #include <this_file>
!>
!> where
!>
!> @item _BF_TYPE_NAME		Type name of the basis
!> @item _BF_SIZE			Number of DoFs of the basis
!> @item _QR_TYPE_NAME		Type name for the quadrature rule
!>
!> The resulting quadrature rule type is defined as _QR_TYPE_NAME
!> @author Oliver Meister

#define _PREFIX(P, X)			_conc3(P,_,X)

#define _BF_(X)					_PREFIX(_BF_TYPE_NAME,X)
#define _BF						_BF_TYPE_NAME

#define _QR_(X)					_PREFIX(_QR_TYPE_NAME,X)
#define _QR						_QR_TYPE_NAME

use SFC_data_types
use Tools_quadrature_rule_base

type, extends(t_qr_base)	:: _QR
end type

PRIVATE
PUBLIC :: _QR, t_qr_create_dunavant_rule, t_qr_create, t_qr_destroy

interface t_qr_create_dunavant_rule
	module procedure _QR_(create_dunavant_rule)
end interface

interface t_qr_create
	module procedure _QR_(create)
end interface

interface t_qr_destroy
	module procedure _QR_(destroy)
end interface

CONTAINS

!**********************
!Constructor/Destructor
!**********************

!> Creates a dunavant quadrature rule for a local basis
subroutine _QR_(create_dunavant_rule)(qr, i_order)
	implicit none
	type(_QR), intent(inout)						:: qr			!< quadrature formula
	integer (kind = GRID_SI), intent(in)			:: i_order		!< order of precision

	integer (kind = GRID_SI)						:: i_error

	!get the number of quadrature points for the given order
	call dunavant_order_num(i_order, qr%i_qpts)

	allocate(qr%r_qpt_coords(2, qr%i_qpts), stat = i_error); assert_eq(i_error, 0)
	allocate(qr%r_qpt_weights(qr%i_qpts), stat = i_error); assert_eq(i_error, 0)

	!get a quadrature rule for the respective order
	call dunavant_rule(i_order, qr%i_qpts, qr%r_qpt_coords, qr%r_qpt_weights)

	!multiply the weights by 0.5 to account for the template triangle area 0.5
	qr%r_qpt_weights(:) = 0.5_GRID_SR * qr%r_qpt_weights(:)

	!set all normals to 0
	allocate(qr%r_qpt_normals(2, qr%i_qpts), stat = i_error); assert_eq(i_error, 0)
	qr%r_qpt_normals(:, :) = 0.0_GRID_SR

	call _QR_(create)(qr)
end subroutine

!> Creates an arbitrary quadrature rule for a local basis. Number, positions and weights of quadrature points must be provided.
subroutine _QR_(create)(qr)
	implicit none
	type(_QR), intent(inout)						:: qr			!< quadrature formula

	real (kind = GRID_SR), dimension(2)				:: x
	integer (kind = GRID_SI)						:: i_qpt, i_error

	!init arrays
	allocate(qr%r_dof_qpt_psi(_BF_SIZE, qr%i_qpts), stat = i_error); assert_eq(i_error, 0)
	allocate(qr%r_dof_qpt_d_psi_d_lambda1(_BF_SIZE, qr%i_qpts), stat = i_error); assert_eq(i_error, 0)
	allocate(qr%r_dof_qpt_d_psi_d_lambda2(_BF_SIZE, qr%i_qpts), stat = i_error); assert_eq(i_error, 0)

	!evaluate basis functions and derivations at the quadrature points
	do i_qpt = 1, qr%i_qpts
		qr%r_dof_qpt_psi(:, i_qpt) = _BF_(at)(qr%r_qpt_coords(:, i_qpt))
		qr%r_dof_qpt_d_psi_d_lambda1(:, i_qpt) = _BF_(d_dx)(qr%r_qpt_coords(:, i_qpt))
		qr%r_dof_qpt_d_psi_d_lambda2(:, i_qpt) = _BF_(d_dy)(qr%r_qpt_coords(:, i_qpt))
	end do
end subroutine

subroutine _QR_(destroy)(qr)
	type(_QR), intent(inout)						:: qr			!< quadrature formula

	call t_qr_base_destroy(qr%t_qr_base)
end subroutine

#undef _BF
#undef _QR
#undef _BF_TYPE_NAME
#undef _BF_SIZE
#undef _QR_TYPE_NAME

