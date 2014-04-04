! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


!> Generic lagrange basis
!> Warning: this a template module body and requires preprocessor commands before inclusion.
!> Usage: Inside your module definition, insert
!>
!> #define _BF_TYPE_NAME	<Type name for the basis>
!> #define _BF_ORDER		<Order of the lagrange basis>
!> #define _BF_type			<DoF and value type, must be real or complex>
!>
!> #include <this_file>
!>
!> The resulting basis type is defined as _BF_TYPE_NAME
!> @author Oliver Meister

#define _PREFIX(P, X)			_conc3(P,_,X)
#define _BF_(X)					_PREFIX(_BF_TYPE_NAME,X)

#define _TO_STRING(Y)			_stringify(Y)
#define _BF_STRING				_TO_STRING(_BF_TYPE_NAME)

#define _BF						_BF_TYPE_NAME
#define _BF_SIZE				(_BF_ORDER + 2) * (_BF_ORDER + 1) / 2

#if !defined(_BF_type)
#	define _BF_type				real(kind = GRID_SR)
#endif

integer (kind = GRID_SI), parameter									:: i_dofs = _BF_SIZE	!< Number of degrees of freedom

!> Positions of the degrees of freedom
#if (_BF_ORDER == 0)
#	define _PSI(x, y) 			[ 1.0_GRID_SR ]
#	define _DPSI_DX(x, y) 		[ 0.0_GRID_SR ]
#	define _DPSI_DY(x, y) 		[ 0.0_GRID_SR ]

	real (kind = GRID_SR), DIMENSION(2, _BF_SIZE), parameter		:: r_dof_coords = reshape( &
		[	1.0_GRID_SR / 3.0_GRID_SR, 1.0_GRID_SR / 3.0_GRID_SR ], &
		[ 2, _BF_SIZE ] )
#elif (_BF_ORDER == 1)
#	define _PSI(x, y) 			[ x, 1.0_GRID_SR - x - y, y ]
#	define _DPSI_DX(x, y) 		[ 1.0_GRID_SR, -1.0_GRID_SR, 0.0_GRID_SR ]
#	define _DPSI_DY(x, y) 		[ 0.0_GRID_SR, -1.0_GRID_SR, 1.0_GRID_SR ]

	real (kind = GRID_SR), DIMENSION(2, _BF_SIZE), parameter		:: r_dof_coords = reshape( &
		[	1.0_GRID_SR, 0.0_GRID_SR, &
			0.0_GRID_SR, 0.0_GRID_SR, &
			0.0_GRID_SR, 1.0_GRID_SR ], &
		[ 2, _BF_SIZE ] )
#elif (_BF_ORDER == 2)
#	define _PSI(x, y) 			[ 2.0_GRID_SR * x * (x - 0.5_GRID_SR), 2.0_GRID_SR * (0.5_GRID_SR - x - y) * (1.0_GRID_SR - x - y), 2.0_GRID_SR * y * (y - 0.5_GRID_SR), 4.0_GRID_SR * y * (1.0_GRID_SR - x - y), 4.0_GRID_SR * x * y, 4.0_GRID_SR * x * (1.0_GRID_SR - x - y) ]
#	define _DPSI_DX(x, y) 		[ 4.0_GRID_SR * x - 1.0_GRID_SR, 4.0_GRID_SR * x + 4.0_GRID_SR * y - 3.0_GRID_SR, 0.0_GRID_SR, -4.0_GRID_SR * y, 4.0_GRID_SR * y, 4.0_GRID_SR - 8.0_GRID_SR * x - 4.0_GRID_SR * y ]
#	define _DPSI_DY(x, y) 		[ 0.0_GRID_SR, 4.0_GRID_SR * x + 4.0_GRID_SR * y - 3.0_GRID_SR, 4.0_GRID_SR * y - 1.0_GRID_SR, 4.0_GRID_SR - 4.0_GRID_SR * x - 8.0_GRID_SR * y, 4.0_GRID_SR * x, -4.0_GRID_SR * x ]

	real (kind = GRID_SR), DIMENSION(2, _BF_SIZE), parameter		:: r_dof_coords = reshape( &
		[	1.0_GRID_SR, 0.0_GRID_SR, &
			0.0_GRID_SR, 0.0_GRID_SR, &
			0.0_GRID_SR, 1.0_GRID_SR, &
			0.0_GRID_SR, 0.5_GRID_SR, &
			0.5_GRID_SR, 0.5_GRID_SR, &
			0.5_GRID_SR, 0.0_GRID_SR ], &
		[ 2, _BF_SIZE ] )
#else
#	error: Lagrange basis must be of order less than 3
#endif

PRIVATE
PUBLIC _BF_(at), _BF_(d_dx), _BF_(d_dy), _BF_(dofs_to_values), _BF_(values_to_dofs), _BF_(eval), _BF_(merge), _BF_(split), _BF_(gradient), _BF_(get_dof_coords), _BF_(test)

interface _BF_(merge)
	module procedure _BF_(merge_vector)
	module procedure _BF_(merge_scalar)
end interface

interface _BF_(split)
	module procedure _BF_(split_interpolate)
	module procedure _BF_(split_interpolate_scalar)
	module procedure _BF_(split_lookup)
	module procedure _BF_(split_lookup_scalar)
end interface

contains

!>Evaluates the lagrange basis functions at a point x (Barycentric coordinates)
pure function _BF_(at)(x) result(r_psi)
	implicit none
	real (kind = GRID_SR), DIMENSION(2), intent(in)		:: x						!< position in barycentric coordinates
	real (kind = GRID_SR), DIMENSION(_BF_SIZE)			:: r_psi					!< (out) basis function values

	r_psi = _PSI(x(1), x(2))
end function

!>Evaluates the lagrange basis x-derivatives at a point x (Barycentric coordinates)
pure function _BF_(d_dx)(x) result(r_d_psi_d_lambda1)
	implicit none
	real (kind = GRID_SR), DIMENSION(2), intent(in)		:: x						!< position in barycentric coordinates
	real (kind = GRID_SR), DIMENSION(_BF_SIZE)			:: r_d_psi_d_lambda1		!< (out) lambda1-derivatives of basis functions

	r_d_psi_d_lambda1 = _DPSI_DX(x(1), x(2))
end function

!>Evaluates the lagrange basis y-derivatives at a point x (Barycentric coordinates)
pure function _BF_(d_dy)(x) result(r_d_psi_d_lambda2)
	implicit none
	real (kind = GRID_SR), DIMENSION(2), intent(in)		:: x						!< position in barycentric coordinates
	real (kind = GRID_SR), DIMENSION(_BF_SIZE)			:: r_d_psi_d_lambda2		!< (out) lambda2-derivatives of basis functions

	r_d_psi_d_lambda2 = _DPSI_DY(x(1), x(2))
end function

!**********
!Evaluation
!**********

!>Evaluates a function at the DoF positions
pure function _BF_(dofs_to_values)(r_dofs) result(r_values)
	implicit none
	_BF_type, DIMENSION(_BF_SIZE), intent(in)			:: r_dofs			!< DoFs describing the function
	_BF_type, DIMENSION(_BF_SIZE)						:: r_values			!< function value

	!do nothing
	r_values = r_dofs
end function

!>Returns the DoF description of a function evaluated at the DoF positions
pure function _BF_(values_to_dofs)(r_values) result(r_dofs)
	implicit none
	_BF_type, DIMENSION(_BF_SIZE), intent(in)			:: r_values 		!< function values, overwritten by degrees of freedom
	_BF_type, DIMENSION(_BF_SIZE)						:: r_dofs			!< DoFs describing the function

	!do nothing
	r_dofs = r_values
end function

!> Merges two functions of adjacent elements for adaptive coarsening
subroutine _BF_(merge_vector)(r_dofs1, r_dofs2, r_dofs)
	_BF_type, DIMENSION(:), intent(in)					:: r_dofs1, r_dofs2		!< DoFs of the two input elements
	_BF_type, DIMENSION(:), intent(out)					:: r_dofs				!< output DoFs

#	if (_BF_ORDER == 0)
		r_dofs = 0.5_GRID_SR * (r_dofs1 + r_dofs2)
#	elif (_BF_ORDER == 1)
		assert_eq(r_dofs1(3), r_dofs2(1))
		assert_eq(r_dofs1(2), r_dofs2(2))

		r_dofs(1) = r_dofs1(1)
		r_dofs(2) = r_dofs1(3)
		r_dofs(3) = r_dofs2(3)
#	elif (_BF_ORDER == 2)
		assert_eq(r_dofs1(3), r_dofs2(1))
		assert_eq(r_dofs1(2), r_dofs2(2))
		assert_eq(r_dofs1(4), r_dofs2(6))

		r_dofs(1) = r_dofs1(1)
		r_dofs(2) = r_dofs1(3)
		r_dofs(3) = r_dofs2(3)

		r_dofs(4) = r_dofs2(5)
		r_dofs(5) = r_dofs1(2)
		r_dofs(6) = r_dofs1(5)
#	else
#		error: Lagrange basis must be of order less than 3
#	endif
end subroutine

!> Merges two functions of adjacent elements for adaptive coarsening (scalar variant for order 0)
subroutine _BF_(merge_scalar)(r_dofs1, r_dofs2, r_dofs)
	_BF_type, intent(in)			:: r_dofs1, r_dofs2		!< DoFs of the two input elements
	_BF_type, intent(out)			:: r_dofs				!< output DoFs

	r_dofs = 0.5_GRID_SR * (r_dofs1 + r_dofs2)
end subroutine

!> Splits a function into two adjacent elements for adaptive refinement by interpolation
pure subroutine _BF_(split_interpolate)(r_dofs, r_dofs1, r_dofs2)
	_BF_type, DIMENSION(:), intent(in)			:: r_dofs				!< DoFs of the two input elements
	_BF_type, DIMENSION(:), intent(out)			:: r_dofs1, r_dofs2		!< output DoFs

#	if (_BF_ORDER == 0)
		r_dofs1 = r_dofs
		r_dofs2 = r_dofs
#	elif (_BF_ORDER == 1)
		r_dofs1(1) = r_dofs(1)
		r_dofs1(2) = DOT_PRODUCT(_PSI(0.5_GRID_SR, 0.5_GRID_SR), r_dofs)
		r_dofs1(3) = r_dofs(2)

		r_dofs2(1) = r_dofs(2)
		r_dofs2(2) = r_dofs1(2)
		r_dofs2(3) = r_dofs(3)
#	elif (_BF_ORDER == 2)
		r_dofs1(1) = r_dofs(1)
		r_dofs1(2) = r_dofs(5)
		r_dofs1(3) = r_dofs(2)
		r_dofs1(4) = DOT_PRODUCT(_PSI(0.25_GRID_SR, 0.25_GRID_SR), r_dofs)
		r_dofs1(5) = r_dofs(6)
		r_dofs1(6) = DOT_PRODUCT(_PSI(0.75_GRID_SR, 0.25_GRID_SR), r_dofs)

		r_dofs2(1) = r_dofs(2)
		r_dofs2(2) = r_dofs1(2)
		r_dofs2(3) = r_dofs(3)
		r_dofs2(4) = DOT_PRODUCT(_PSI(0.25_GRID_SR, 0.75_GRID_SR), r_dofs)
		r_dofs2(5) = r_dofs(4)
		r_dofs2(6) = r_dofs1(4)
#	else
#		error: Lagrange basis must be of order less than 3
#	endif
end subroutine

!> Splits a function into two adjacent elements for adaptive refinement by interpolation (scalar variant for order 0)
pure subroutine _BF_(split_interpolate_scalar)(r_dofs, r_dofs1, r_dofs2)
	_BF_type, intent(in)			:: r_dofs				!< DoFs of the two input elements
	_BF_type, intent(out)			:: r_dofs1, r_dofs2		!< output DoFs

	r_dofs1 = r_dofs
	r_dofs2 = r_dofs
end subroutine

!> Splits a function into two adjacent elements for adaptive refinement by lookup
subroutine _BF_(split_lookup)(r_dofs, r_dofs1, r_dofs2, r_mat_transform, r_offset_transform, lookup)
	_BF_type, DIMENSION(:), intent(in)			:: r_dofs				!< DoFs of the two input elements
	_BF_type, DIMENSION(_BF_SIZE), intent(out)			:: r_dofs1, r_dofs2		!< output DoFs
	real (kind = GRID_SR), DIMENSION(2, 2), intent(in)				:: r_mat_transform		!< local to global coordinate transformation matrix
	real (kind = GRID_SR), DIMENSION(2), intent(in)					:: r_offset_transform	!< local to global coordinate transformation offset

	interface
		!< Returns a function value for the point x in global coordinates
		function lookup(x) result(r_value)
			import
			real (kind = GRID_SR), dimension(2), intent(in)		:: x						!< position in world coordinates
			_BF_type											:: r_value					!< value
		end function
	end interface

#	if (_BF_ORDER == 0)
		r_dofs1 = lookup(MATMUL(r_mat_transform, [ 0.5_GRID_SR, 1.0_GRID_SR / 6.0_GRID_SR ]) + r_offset_transform)
		r_dofs2 = lookup(MATMUL(r_mat_transform, [ 1.0_GRID_SR / 6.0_GRID_SR, 0.5_GRID_SR ]) + r_offset_transform)
#	elif (_BF_ORDER == 1)
		r_dofs1(1) = r_dofs(2)
		r_dofs1(2) = lookup(MATMUL(r_mat_transform, [ 0.5_GRID_SR, 0.5_GRID_SR ]) + r_offset_transform)
		r_dofs1(3) = r_dofs(1)

		r_dofs2(1) = r_dofs(3)
		r_dofs2(2) = r_dofs1(2)
		r_dofs2(3) = r_dofs(2)
#	elif (_BF_ORDER == 2)
		r_dofs1(1) = r_dofs(1)
		r_dofs1(2) = r_dofs(5)
		r_dofs1(3) = r_dofs(2)
		r_dofs1(4) = lookup(MATMUL(r_mat_transform, [ 0.25_GRID_SR, 0.25_GRID_SR ]) + r_offset_transform)
		r_dofs1(5) = r_dofs(6)
		r_dofs1(6) = lookup(MATMUL(r_mat_transform, [ 0.75_GRID_SR, 0.25_GRID_SR ]) + r_offset_transform)

		r_dofs2(1) = r_dofs(2)
		r_dofs2(2) = r_dofs1(2)
		r_dofs2(3) = r_dofs(3)
		r_dofs2(4) = lookup(MATMUL(r_mat_transform, [ 0.25_GRID_SR, 0.75_GRID_SR ]) + r_offset_transform)
		r_dofs2(5) = r_dofs(4)
		r_dofs2(6) = r_dofs1(4)
#	else
#		error: Lagrange basis must be of order less than 3
#	endif
end subroutine

!> Splits a function into two adjacent elements for adaptive refinement by lookup (scalar variant for order 0)
subroutine _BF_(split_lookup_scalar)(r_dofs, r_dofs1, r_dofs2, r_mat_transform, r_offset_transform, lookup)
	_BF_type, intent(in)								:: r_dofs				!< DoFs of the two input elements
	_BF_type, intent(out)								:: r_dofs1, r_dofs2		!< output DoFs
	real (kind = GRID_SR), DIMENSION(2, 2), intent(in)				:: r_mat_transform		!< local to global coordinate transformation matrix
	real (kind = GRID_SR), DIMENSION(2), intent(in)					:: r_offset_transform	!< local to global coordinate transformation offset

	interface
		!< Returns a function value for the point x in global coordinates
		function lookup(x) result(r_value)
			import
			real (kind = GRID_SR), dimension(2), intent(in)		:: x						!< position in world coordinates
			_BF_type								:: r_value					!< value
		end function
	end interface

	r_dofs1 = lookup(MATMUL(r_mat_transform, [ 0.5_GRID_SR, 1.0_GRID_SR / 6.0_GRID_SR ]) + r_offset_transform)
	r_dofs2 = lookup(MATMUL(r_mat_transform, [ 1.0_GRID_SR / 6.0_GRID_SR, 0.5_GRID_SR ]) + r_offset_transform)
end subroutine

!>Evaluates a function at a point x (Barycentric coordinates)
pure function _BF_(eval)(x, r_dofs) result(r_value)
	implicit none
	real (kind = GRID_SR), DIMENSION(:), intent(in)					:: x				!< position in barycentric coordinates
	real (kind = GRID_SR), DIMENSION(:), intent(in)					:: r_dofs			!< DoFs describing the function
	_BF_type														:: r_value			!< function value

	r_value = DOT_PRODUCT(_PSI(x(1), x(2)), r_dofs)
end function

!>Derives a function at a point x (Barycentric coordinates)
pure function _BF_(gradient)(x, r_dofs) result(r_value)
	implicit none
	real (kind = GRID_SR), DIMENSION(:), intent(in)					:: x				!< position in barycentric coordinates
	_BF_type, DIMENSION(:), intent(in)								:: r_dofs			!< DoFs describing the function
	_BF_type, DIMENSION(2)											:: r_value			!< barycentric gradient of the function

	r_value = [ DOT_PRODUCT(_DPSI_DX(x(1), x(2)), r_dofs), DOT_PRODUCT(_DPSI_DY(x(1), x(2)), r_dofs) ]
end function

!>Returns the DoF coordinates
pure function _BF_(get_dof_coords)(i_dof) result(r_coords)
	implicit none
	integer (kind = GRID_SI), intent(in)							:: i_dof		!< DoF index
	real (kind = GRID_SR), DIMENSION(2)								:: r_coords		!< DoF coords

	r_coords(:) = r_dof_coords(:, i_dof)
end function

subroutine _BF_(test)()
	implicit none
	integer (kind = GRID_SI)								:: i, j
	real (kind = GRID_SR), dimension(2), parameter			:: x = [ 0.25_GRID_SR, 0.0_GRID_SR ]
	real (kind = GRID_SR), dimension(_BF_SIZE), parameter	:: psi_x = _PSI(x(1), x(2))
	real (kind = GRID_SR), dimension(_BF_SIZE)				:: r_dofs = [ (dble(i), i = 1, _BF_SIZE) ]
	real (kind = GRID_SR), dimension(_BF_SIZE)				:: r_dofs1, r_dofs2
	real (kind = GRID_SR), dimension(_BF_SIZE)				:: r_values
	real (kind = GRID_SR)									:: r_value
	real (kind = GRID_SR)									:: t1, t2, t3

	_log_write(0, *) "Test: ", _BF_STRING

	!test for coordinate <-> basis function consistency

	do i = 1, _BF_SIZE
		r_values = _BF_(at)(r_dof_coords(:, i))

		do j = 1, _BF_SIZE
			if (i == j) then
				assert_eq(r_values(j), 1.0_GRID_SR)
			else
				assert_eq(r_values(j), 0.0_GRID_SR)
			endif
		end do
	end do

	!test for split correctness

	call _BF_(split)(r_dofs, r_dofs1, r_dofs2)

	do i = 1, _BF_SIZE
		assert_eq(_BF_(eval)(r_dof_coords(:, i), r_dofs1), _BF_(eval)([ 0.5_GRID_SR + 0.5_GRID_SR * r_dof_coords(1, i) - 0.5_GRID_SR * r_dof_coords(2, i), 0.5_GRID_SR - 0.5_GRID_SR * r_dof_coords(1, i) - 0.5_GRID_SR * r_dof_coords(2, i) ], r_dofs))
		assert_eq(_BF_(eval)(r_dof_coords(:, i), r_dofs2), _BF_(eval)([ 0.5_GRID_SR - 0.5_GRID_SR * r_dof_coords(1, i) - 0.5_GRID_SR * r_dof_coords(2, i), 0.5_GRID_SR - 0.5_GRID_SR * r_dof_coords(1, i) + 0.5_GRID_SR * r_dof_coords(2, i) ], r_dofs))
	end do

	!test for merge correctness

	call _BF_(merge)(r_dofs1, r_dofs2, r_values)

	do i = 1, _BF_SIZE
		assert_eq(r_dofs(i), r_values(i))
	end do

	!timing tests
	_log_write(0, *) " dofs_to_values performance:"

	r_values = r_dofs
	call cpu_time(t1)
	do i = 1, 6*(10**7)
		r_values = _BF_(dofs_to_values)(r_values)
	end do
	call cpu_time(t2)
	_log_write(0, "(A, F0.4, G0.4)") "   dtv                : ", (t2 - t1), sum(r_values)

	_log_write(0, *) " ct-constant position evaluation:"

	r_value = 0.0_GRID_SR
	call cpu_time(t1)
	do i = 1, 6*(10**7)
		r_value = r_value + _BF_(eval)(x, r_dofs)
	end do
	call cpu_time(t2)
	_log_write(0, "(A, F0.4, G0.4)") "   eval(x, dofs)      : ", (t2 - t1), r_value

	r_value = 0.0_GRID_SR
	call cpu_time(t1)
	do i = 1, 6*(10**7)
		r_value = r_value + DOT_PRODUCT(_BF_(at)(x), r_dofs)
	end do
	call cpu_time(t2)
	_log_write(0, "(A, F0.4, G0.4)") "   dot(at(x), dofs)   : ", (t2 - t1), r_value

	r_value = 0.0_GRID_SR
	call cpu_time(t1)
	do i = 1, 6*(10**7)
		r_value = r_value + DOT_PRODUCT(_PSI(x(1), x(2)), r_dofs)
	end do
	call cpu_time(t2)
	_log_write(0, "(A, F0.4, G0.4)") "   dot(_PSI(x), dofs) : ", (t2 - t1), r_value

	r_value = 0.0_GRID_SR
	call cpu_time(t1)
	do i = 1, 6*(10**7)
		r_value = r_value + DOT_PRODUCT(psi_x, r_dofs)
	end do
	call cpu_time(t2)
	_log_write(0, "(A, F0.4, G0.4)") "   dot(psi_x, dofs)   : ", (t2 - t1), r_value

	r_value = 0.0_GRID_SR
	call cpu_time(t1)
	do i = 1, 6*(10**7)
#		if (_BF_ORDER == 0)
			r_value = r_value + r_dofs(1)
#		elif (_BF_ORDER == 1)
			r_value = r_value + (0.75_GRID_SR * r_dofs(2) + 0.25_GRID_SR * r_dofs(1))
#		elif (_BF_ORDER == 2)
			r_value = r_value + (0.375_GRID_SR * r_dofs(2) + 0.75_GRID_SR * r_dofs(6) - 0.125_GRID_SR * r_dofs(1))
#		endif
	end do
	call cpu_time(t2)
	_log_write(0, "(A, F0.4, G0.4)") "   hard-coded         : ", (t2 - t1), r_value
end subroutine

#undef _BF
#undef _BF_SIZE
#undef _PSI
#undef _DPSI_DX
#undef _DPSI_DY

#undef _BF_type
#undef _BF_TYPE_NAME
#undef _BF_ORDER
