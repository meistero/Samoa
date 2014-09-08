! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


!> Generic 1D lagrange basis
!> Warning: this a template module body and requires preprocessor commands before inclusion.
!> Usage: Inside your module definition, insert
!>
!> #define _BF_TYPE_NAME	<value>
!> #define _BF_ORDER		<value>
!>
!> #include <this_file>
!>
!> where
!>
!> @item _BF_TYPE_NAME		Type name for the basis
!> @item _BF_ORDER			Order of the lagrange basis
!>
!> The resulting basis type is defined as _BF_TYPE_NAME
!> @author Oliver Meister

#define _PREFIX(P, X)			_conc3(P,_,X)
#define _BF_(X)					_PREFIX(_BF_TYPE_NAME,X)

#define _BF						_BF_TYPE_NAME
#define _BF_SIZE				_BF_ORDER + 1

integer (kind = GRID_SI), parameter									:: i_dofs = _BF_SIZE	!< Number of degrees of freedom

!> Positions of the degrees of freedom
#if (_BF_ORDER == 0)
#	define _PSI(x)												[1.0_GRID_SR]
#	define _DPSI_DX(x)											[0.0_GRID_SR]
	real (kind = GRID_SR), DIMENSION(_BF_SIZE), parameter		:: r_dof_coords = [0.5_GRID_SR]
#elif (_BF_ORDER == 1)
#	define _PSI(x)												[1.0_GRID_SR - x, x]
#	define _DPSI_DX(x)											[-1.0_GRID_SR, 1.0_GRID_SR]
	real (kind = GRID_SR), DIMENSION(_BF_SIZE), parameter		:: r_dof_coords = [ 0.0_GRID_SR, 1.0_GRID_SR ]
#elif (_BF_ORDER == 2)
#	define _PSI(x)												[2.0_GRID_SR * x * x - 3.0_GRID_SR * x + 1.0_GRID_SR, -4.0 * x * x + 4.0 * x, 2.0 * x * x - x]
#	define _DPSI_DX(x)											[4.0_GRID_SR * x - 3.0_GRID_SR, -8.0_GRID_SR * x + 4.0_GRID_SR, 4.0_GRID_SR * x - 1.0_GRID_SR]
	real (kind = GRID_SR), DIMENSION(_BF_SIZE), parameter		:: r_dof_coords =  [ 0.0_GRID_SR, 0.5_GRID_SR, 1.0_GRID_SR ]
#else
#	error Boundary basis must be of order 2 or less
#endif

PRIVATE
PUBLIC _BF_(at), _BF_(d_dx), _BF_(dofs_to_values), _BF_(values_to_dofs), _BF_(eval), _BF_(merge), _BF_(split), _BF_(gradient), _BF_(get_dof_coords), _BF_(test)

contains

!>Evaluates the basis at a point x (Barycentric coordinates)
pure function _BF_(at)(x) result(r_psi)
	implicit none
	real (kind = GRID_SR), intent(in)							:: x						!< position in barycentric coordinates
	real (kind = GRID_SR), DIMENSION(_BF_SIZE)					:: r_psi					!< (out) basis function values

	r_psi = _PSI(x)
end function

pure function _BF_(d_dx)(x) result(r_dpsi_dx)
	implicit none
	real (kind = GRID_SR), intent(in)							:: x						!< position in barycentric coordinates
	real (kind = GRID_SR), DIMENSION(_BF_SIZE)					:: r_dpsi_dx				!< (out) x-derivatives of basis functions

	r_dpsi_dx = _DPSI_DX(x)
end function

!**********
!Evaluation
!**********

!>Evaluates a function at the position of DoF i
pure function _BF_(dofs_to_values)(r_dofs) result(r_values)
	implicit none
	real (kind = GRID_SR), DIMENSION(:), intent(in)			:: r_dofs			!< DoFs describing the function
	real (kind = GRID_SR), DIMENSION(_BF_SIZE)						:: r_values			!< function value

	!do nothing
	r_values = r_dofs
end function

!>Returns the DoF description of a function evaluated at the DoF positions
pure function _BF_(values_to_dofs)(r_values) result(r_dofs)
	implicit none
	real (kind = GRID_SR), DIMENSION(:), intent(in)			:: r_values 	!< function values, overwritten by degrees of freedom
	real (kind = GRID_SR), DIMENSION(_BF_SIZE)						:: r_dofs			!< DoFs describing the function

	!do nothing
	r_dofs = r_values
end function

!> Merges two functions of adjacent elements for adaptive coarsening
pure subroutine _BF_(merge)(r_dofs1, r_dofs2, r_dofs)
	real (kind = GRID_SR), DIMENSION(:), intent(in)			:: r_dofs1, r_dofs2		!< DoFs of the two input elements
	real (kind = GRID_SR), DIMENSION(:), intent(out)			:: r_dofs				!< output DoFs

#	if (_BF_ORDER == 0)
		r_dofs = 0.5_GRID_SR * (r_dofs1 + r_dofs2)
#	elif (_BF_ORDER == 1)
		r_dofs(1) = r_dofs1(1)
		r_dofs(2) = r_dofs2(2)
#	elif (_BF_ORDER == 2)
		r_dofs(1) = r_dofs1(1)
		r_dofs(2) = 0.5_GRID_SR * (r_dofs1(3) + r_dofs2(1))
		r_dofs(3) = r_dofs2(3)
#	else
#		error Boundary basis must be of order 2 or less
#	endif
end subroutine

!> Splits a function into two adjacent elements for adaptive refinement
pure subroutine _BF_(split)(r_dofs, r_dofs1, r_dofs2)
	real (kind = GRID_SR), DIMENSION(:), intent(in)			:: r_dofs				!< DoFs of the two input elements
	real (kind = GRID_SR), DIMENSION(:), intent(out)			:: r_dofs1, r_dofs2		!< output DoFs

#	if (_BF_ORDER == 0)
		r_dofs1 = r_dofs
		r_dofs2 = r_dofs
#	elif (_BF_ORDER == 1)
		r_dofs1(1) = r_dofs(1)
		r_dofs1(2) = 0.5_GRID_SR * (r_dofs(1) + r_dofs(2))
		r_dofs2(1) = 0.5_GRID_SR * (r_dofs(1) + r_dofs(2))
		r_dofs2(2) = r_dofs(2)
#	elif (_BF_ORDER == 2)
		r_dofs1(1) = r_dofs(1)
		r_dofs1(2) = DOT_PRODUCT(r_dofs, _BF_(at)(0.25_GRID_SR))
		r_dofs1(3) = r_dofs(2)
		r_dofs2(1) = r_dofs(2)
		r_dofs2(2) = DOT_PRODUCT(r_dofs, _BF_(at)(0.75_GRID_SR))
		r_dofs2(3) = r_dofs(3)
#	else
#		error Boundary basis must be of order 2 or less
#	endif
end subroutine

!>Evaluates a function at a point x (Barycentric coordinates)
pure function _BF_(eval)(x, r_dofs) result(r_value)
	implicit none
	real (kind = GRID_SR), intent(in)								:: x				!< position in barycentric coordinates
	real (kind = GRID_SR), DIMENSION(_BF_SIZE), intent(in)			:: r_dofs			!< DoFs describing the function
	real (kind = GRID_SR)											:: r_value			!< function value

	r_value = DOT_PRODUCT(_PSI(x), r_dofs)
end function

!>Evaluates a function at a point x (Barycentric coordinates)
pure function _BF_(gradient)(x, r_dofs) result(r_value)
	implicit none
	real (kind = GRID_SR), intent(in)								:: x				!< position in barycentric coordinates
	real (kind = GRID_SR), DIMENSION(_BF_SIZE), intent(in)			:: r_dofs			!< DoFs describing the function
	real (kind = GRID_SR)											:: r_value			!< function value

	r_value = DOT_PRODUCT(_DPSI_DX(x), r_dofs)
end function

!>Returns the DoF coordinates
pure function _BF_(get_dof_coords)(i_dof) result(r_coords)
	implicit none
	integer (kind = GRID_SI), intent(in)							:: i_dof		!< DoF index
	real (kind = GRID_SR)											:: r_coords		!< DoF coords

	r_coords = r_dof_coords(i_dof)
end function

subroutine _BF_(test)()
	implicit none
	integer (kind = GRID_SI)								:: i, j
	real (kind = GRID_SR), parameter						:: x = 0.5_GRID_SR
	real (kind = GRID_SR), dimension(_BF_SIZE), parameter	:: psi_x = _PSI(x)
	real (kind = GRID_SR), dimension(_BF_SIZE)				:: r_dofs = [ (dble(i), i = 1, _BF_SIZE) ]
	real (kind = GRID_SR), dimension(_BF_SIZE)				:: r_dofs1, r_dofs2
	real (kind = GRID_SR), dimension(_BF_SIZE)				:: r_values
	real (kind = GRID_SR)									:: r_value
	real (kind = GRID_SR)									:: t1, t2, t3

	_log_write(0, *) "Test: ", _BF_STRING

	!test for split -> merge correctness

	call _BF_(split)(r_dofs, r_dofs1, r_dofs2)
	call _BF_(merge)(r_dofs1, r_dofs2, r_values)

	do i = 1, _BF_SIZE
		assert_eq(r_dofs(i), r_values(i))
	end do

	!timing tests (not in debug mode)
	_log_write(0, *) " dofs_to_values performance:"

	r_values = r_dofs
	call cpu_time(t1)
	do i = 1, (10**8)
		r_values = _BF_(dofs_to_values)(r_values)
	end do
	call cpu_time(t2)
	_log_write(0, "(A, F7.4, G0.4)") "   dtv                : ", (t2 - t1), sum(r_values)

	_log_write(0, *) " ct-constant position evaluation:"

	r_value = 0.0_GRID_SR
	call cpu_time(t1)
	do i = 1, (10**8)
		r_value = r_value + _BF_(eval)(x, r_dofs)
	end do
	call cpu_time(t2)
	_log_write(0, "(A, F7.4, G0.4)") "   eval(x, dofs)      : ", (t2 - t1), r_value

	r_value = 0.0_GRID_SR
	call cpu_time(t1)
	do i = 1, (10**8)
		r_value = r_value + DOT_PRODUCT(_BF_(at)(x), r_dofs)
	end do
	call cpu_time(t2)
	_log_write(0, "(A, F7.4, G0.4)") "   dot(at(x), dofs)   : ", (t2 - t1), r_value

	r_value = 0.0_GRID_SR
	call cpu_time(t1)
	do i = 1, (10**8)
		r_value = r_value + DOT_PRODUCT(_PSI(x), r_dofs)
	end do
	call cpu_time(t2)
	_log_write(0, "(A, F7.4, G0.4)") "   dot(_PSI(x), dofs) : ", (t2 - t1), r_value

	r_value = 0.0_GRID_SR
	call cpu_time(t1)
	do i = 1, (10**8)
		r_value = r_value + DOT_PRODUCT(psi_x, r_dofs)
	end do
	call cpu_time(t2)
	_log_write(0, "(A, F7.4, G0.4)") "   dot(psi_x, dofs)   : ", (t2 - t1), r_value

	r_value = 0.0_GRID_SR
	call cpu_time(t1)
	do i = 1, (10**8)
#		if (_BF_ORDER == 0)
			r_value = r_value + r_dofs(1)
#		elif (_BF_ORDER == 1)
			r_value = r_value + (0.5_GRID_SR * r_dofs(1) + 0.5_GRID_SR * r_dofs(2))
#		elif (_BF_ORDER == 2)
			r_value = r_value + r_dofs(2)
#		endif
	end do
	call cpu_time(t2)
	_log_write(0, "(A, F7.4, G0.4)") "   hard-coded         : ", (t2 - t1), r_value
end subroutine

#undef _BF
#undef _BF_SIZE
#undef _PSI
#undef _DPSI_DX
#undef _DPSI_DY

#undef _BF_type
#undef _BF_TYPE_NAME
#undef _BF_ORDER
