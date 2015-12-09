! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

MODULE Samoa
	use SFC_data_types
	use Tools_quadrature_rule_base

    implicit none

	public

	contains

	!*******************************************
	!Barycentric <-> World coordinate conversion
	!*******************************************

	!> Transforms a point from barycentric to world coordinates
	pure function samoa_barycentric_to_world_point(td, x) result(r_pos)
		type(t_transform_data), intent(in)						:: td
		real (kind = GRID_SR), DIMENSION(:), intent(in)			:: x
		real (kind = GRID_SR), DIMENSION(2)						:: r_pos

		!> computes s A x + b

		r_pos = td%custom_data%scaling * MATMUL(td%plotter_data%jacobian, x) + td%custom_data%offset
	end function

	!> Transforms a vector from barycentric to world coordinates
	pure function samoa_barycentric_to_world_vector(td, x) result(r_pos)
		type(t_transform_data), intent(in)						:: td
		real (kind = GRID_SR), DIMENSION(:), intent(in)			:: x
		real (kind = GRID_SR), DIMENSION(2)						:: r_pos

		!> computes s A x
		r_pos = td%custom_data%scaling * MATMUL(td%plotter_data%jacobian, x)
	end function

	!> Transforms a matrix from barycentric to world coordinates
	pure function samoa_barycentric_to_world_matrix(td, m_b) result(m_w)
		type(t_transform_data), intent(in)  :: td
		real (kind = GRID_SR), intent(in)   :: m_b(:, :)
		real (kind = GRID_SR)               :: m_w(2, 2)

		!> computes s A x
		m_w = td%custom_data%scaling * MATMUL(td%plotter_data%jacobian, m_b)
	end function

	!> Transforms a normal from barycentric to world coordinates (length-preserving)
	pure function samoa_barycentric_to_world_normal(td, x) result(r_pos)
		type(t_transform_data), intent(in)						:: td
		real (kind = GRID_SR), DIMENSION(:), intent(in)			:: x
		real (kind = GRID_SR), DIMENSION(2)						:: r_pos

		!> computes A^(-T) / sqrt(|A^-T|) x
		r_pos = MATMUL(x, td%plotter_data%jacobian_inv) * sqrt(abs(td%plotter_data%det_jacobian))
	end function

	!> Transforms a point from world to barycentric coordinates
	pure function samoa_world_to_barycentric_point(td, x) result(r_pos)
		type(t_transform_data), intent(in)						:: td
		real (kind = GRID_SR), DIMENSION(:), intent(in)			:: x
		real (kind = GRID_SR), DIMENSION(2)						:: r_pos

		!> computes 1/s A^(-1) (x - b)
		r_pos = MATMUL(td%plotter_data%jacobian_inv, x - td%custom_data%offset) / td%custom_data%scaling
	end function

	!> Transforms a vector from world to barycentric coordinates
	pure function samoa_world_to_barycentric_vector(td, x) result(r_pos)
		type(t_transform_data), intent(in)						:: td
		real (kind = GRID_SR), DIMENSION(:), intent(in)			:: x
		real (kind = GRID_SR), DIMENSION(2)						:: r_pos

		!> computes 1/s A^(-1) x
		r_pos = MATMUL(td%plotter_data%jacobian_inv, x) / td%custom_data%scaling
	end function

	!> Transforms a matrix from world to barycentric coordinates
	pure function samoa_world_to_barycentric_matrix(td, m_w) result(m_b)
		type(t_transform_data), intent(in)						:: td
		real (kind = GRID_SR), intent(in)			            :: m_w(:, :)
		real (kind = GRID_SR)					                :: m_b(2, 2)

		!> computes 1/s A^(-1) x
		m_b = MATMUL(td%plotter_data%jacobian_inv, m_w) / td%custom_data%scaling
	end function

	!> Transforms a normal from world to barycentric coordinates (length-preserving)
	pure function samoa_world_to_barycentric_normal(td, x) result(r_pos)
		type(t_transform_data), intent(in)						:: td
		real (kind = GRID_SR), DIMENSION(:), intent(in)			:: x
		real (kind = GRID_SR), DIMENSION(2)						:: r_pos

		!> computes A^T / sqrt(|A^T|) x
		r_pos = MATMUL(x, td%plotter_data%jacobian) / sqrt(abs(td%plotter_data%det_jacobian))
	end function

	!*********************
	!Integration Operators
	!*********************

	!> Calculates an element matrix over the reference element
	pure subroutine samoa_calc_element_matrix(td, qr, samoa_mat_op, mat)
		type(t_transform_data), intent(in)						:: td
		type(t_qr_base), intent(in)								:: qr
		real (kind = GRID_SR), DIMENSION(:,:), intent(out)		:: mat

		interface
			!> interface for an element matrix operator
			pure function samoa_mat_op(td, qr, i, j, i_qpt) result(r_value)
				import
				type(t_transform_data), intent(in)				:: td
				type(t_qr_base), intent(in)						:: qr
				integer (kind = GRID_SI), intent(in)			:: i, j
				integer (kind = GRID_SI), intent(in)			:: i_qpt
				real (kind = GRID_SR)							:: r_value
			end function
		end interface

		integer (kind = GRID_SI)								:: i, j
		integer (kind = GRID_SI)								:: i_qpt

		mat = 0.0_GRID_SR

		!loop through quadrature points and DoFs
		do i_qpt = 1, qr%i_qpts
			forall (i = lbound(mat, 1) : ubound(mat, 1), j = lbound(mat, 2) : ubound(mat, 2))
				mat(i, j) = mat(i, j) + qr%r_qpt_weights(i_qpt) * samoa_mat_op(td, qr, i, j, i_qpt)
			end forall
		end do
	end subroutine

	!> Integrates over the reference element
	pure function samoa_calc_element_integral(td, qr, samoa_int_op) result(r_value)
		type(t_transform_data), intent(in)						:: td
		type(t_qr_base), intent(in)								:: qr
		real (kind = GRID_SR)									:: r_value

		interface
			!> interface for an element integral operator
			pure function samoa_int_op(td, qr, i_qpt) result(r_value)
				import
				type(t_transform_data), intent(in)				:: td
				type(t_qr_base), intent(in)						:: qr
				integer (kind = GRID_SI), intent(in)			:: i_qpt
				real (kind = GRID_SR)							:: r_value
			end function
		end interface

		integer (kind = GRID_SI)								:: i_qpt

		r_value = 0.0_GRID_SR

		!loop through quadrature points to compute the element integral
		do i_qpt = 1, qr%i_qpts
			r_value = r_value + qr%r_qpt_weights(i_qpt) * samoa_int_op(td, qr, i_qpt)
		end do
	end function
END MODULE
