! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_SWE)
	!*****************
	!pressure space
	!*****************
#define _GV_CELL_SIZE		1
#define _GV_EDGE_SIZE       0
#define _GV_NODE_SIZE		0

    MODULE SWE_gm_A_mod
        use SFC_data_types

type swe_gm_A
            contains
            procedure, pass :: read_from_element
            procedure, pass :: write_to_element
            procedure, pass :: apply
            procedure, pass :: get_trace

            generic:: read => read_from_element
            generic:: write=> write_to_element
end type

contains
         pure subroutine read_from_element(gm_A, element, mat)
            class(swe_gm_A), intent(in)	    :: gm_A
            type(t_element_base), intent(in)    :: element
            real (kind = GRID_SR), intent(out)  :: mat(3, 3)
            real (kind=GRID_SR)       :: temp(3)

            mat=element%cell%data_pers%A

            if( (element%transform_data%plotter_data%orientation) /= (element%cell%data_pers%original_lse_orientation(1))) then
                !permute rows
                temp(:)=mat(1,:)
                mat(1,:)=mat(3,:)
                mat(3,:)= temp(:)

                !permute columns
                temp(:)=mat(:,1)
                mat(:,1)=mat(:,3)
                mat(:,3)=temp(:)


            endif


        end subroutine


        pure subroutine write_to_element(gm_A,element, mat)
          class (swe_gm_A), intent(in) :: gm_A
          type(t_element_base), intent(inout)    :: element
          real (kind = GRID_SR), intent(in)  :: mat(3, 3)
            integer (kind = GRID_SI)			    :: i,j


          do i=1,3
            do j=1,3
                element%cell%data_pers%A(i,j)=mat(i,j)
            end do
          end do

        end subroutine

        pure subroutine apply(gm_A, element, x, r)
            class(swe_gm_A), intent(in)         :: gm_A
            type(t_element_base), intent(in)    :: element
            real(kind = GRID_SR), intent(in)    :: x(3)
            real(kind = GRID_SR), intent(inout) :: r(3)

            real(kind = GRID_SR)                :: mat(3, 3)

            r = matmul(element%cell%data_pers%A, x)
        end subroutine

        pure subroutine get_trace(gm_A, element, d)
            class(swe_gm_A), intent(in)         :: gm_A
            type(t_element_base), intent(in)    :: element
            real(kind = GRID_SR), intent(inout) :: d(3)

            real(kind = GRID_SR)                :: mat(3, 3)

            d(1) = element%cell%data_pers%A(1,1)
            d(2) = element%cell%data_pers%A(2,2)
            d(3) = element%cell%data_pers%A(3,3)
        end subroutine
    END MODULE



MODULE SWE_gv_original_lse_orientation_mod
 use SFC_data_types


#		define _GV_TYPE_NAME		swe_gv_original_lse_orientation
#		define _GV_TYPE				integer
#		define _GV_NAME				original_lse_orientation
#		define _GV_PERSISTENT		1

#		include "Tools_grid_variable.f90"
	END MODULE

#	define _GV_CELL_SIZE		_SWE_P_CELL_SIZE
#	define _GV_EDGE_SIZE		_SWE_P_EDGE_SIZE
#	define _GV_NODE_SIZE		_SWE_P_NODE_SIZE


    MODULE SWE_gv_qp_mod
		use SFC_data_types

#		define _GV_TYPE_NAME		swe_gv_qp
#		define _GV_TYPE				real (kind = GRID_SR)
#		define _GV_NAME				qp
#		define _GV_PERSISTENT		1

#		include "Tools_grid_variable.f90"
	END MODULE


MODULE SWE_gv_w_mod
		use SFC_data_types

#		define _GV_TYPE_NAME		swe_gv_w
#		define _GV_TYPE				real (kind = GRID_SR)
#		define _GV_NAME				w
#		define _GV_PERSISTENT		1

#		include "Tools_grid_variable.f90"
	END MODULE

 MODULE SWE_gv_div_mod
		use SFC_data_types

#		define _GV_TYPE_NAME		swe_gv_div
#		define _GV_TYPE				real (kind = GRID_SR)
#		define _GV_NAME				div
#		define _GV_PERSISTENT		1

#		include "Tools_grid_variable.f90"
	END MODULE

	   MODULE SWE_gv_div_old_mod
		use SFC_data_types

#		define _GV_TYPE_NAME		swe_gv_div_old
#		define _GV_TYPE				real (kind = GRID_SR)
#		define _GV_NAME				div_old
#		define _GV_PERSISTENT		1

#		include "Tools_grid_variable.f90"
	END MODULE

    MODULE SWE_gv_mat_diagonal_mod
		use SFC_data_types

#		define _GV_TYPE_NAME		swe_gv_mat_diagonal
#		define _GV_TYPE				real (kind = GRID_SR)
#		define _GV_NAME				mat_diagonal
#		define _GV_PERSISTENT		0

#		include "Tools_grid_variable.f90"
	END MODULE

    MODULE SWE_gv_r_mod
		use SFC_data_types

#		define _GV_TYPE_NAME		swe_gv_r
#		define _GV_TYPE				real (kind = GRID_SR)
#		define _GV_NAME				r
#		define _GV_PERSISTENT		1

#		include "Tools_grid_variable.f90"
	END MODULE

    MODULE SWE_gv_d_mod
		use SFC_data_types

#		define _GV_TYPE_NAME		swe_gv_d
#		define _GV_TYPE				real (kind = GRID_SR)
#		define _GV_NAME				d
#		define _GV_PERSISTENT		1

#		include "Tools_grid_variable.f90"
	END MODULE

    MODULE SWE_gv_A_d_mod
		use SFC_data_types

#		define _GV_TYPE_NAME		swe_gv_A_d
#		define _GV_TYPE				real (kind = GRID_SR)
#		define _GV_NAME				A_d
#		define _GV_PERSISTENT		1

#		include "Tools_grid_variable.f90"
	END MODULE

    MODULE SWE_gv_is_dirichlet_boundary_mod
		use SFC_data_types

#		define _GV_TYPE_NAME		swe_gv_is_dirichlet_boundary
#		define _GV_TYPE				logical
#		define _GV_NAME				is_dirichlet_boundary
#		define _GV_PERSISTENT		1
#		define _GV_ADD_OP			.or.

#		include "Tools_grid_variable.f90"
	END MODULE

	MODULE SWE_gv_rhs_mod
        use SFC_data_types
#		define _GV_TYPE_NAME		swe_gv_rhs
#		define _GV_TYPE				real (kind = GRID_SR)
#		define _GV_NAME				rhs
#		define _GV_PERSISTENT		1

#		include "Tools_grid_variable.f90"
    END MODULE

!undefine macros to avoid compiler warnings
#	undef _GV_CELL_SIZE
#	undef _GV_EDGE_SIZE
#	undef _GV_NODE_SIZE

#	define _GV_CELL_SIZE		_SWE_CELL_SIZE
#	define _GV_EDGE_SIZE		0
#	define _GV_NODE_SIZE		0

	MODULE SWE_gv_Q
		use SFC_data_types

#		define _GV_TYPE_NAME		t_gv_Q
#		define _GV_TYPE				type(t_state)
#		define _GV_NAME				Q
#		define _GV_PERSISTENT		1

#		include "Tools_grid_variable.f90"
	END MODULE

	MODULE SWE_gv_h_old
		use SFC_data_types

#		define _GV_TYPE_NAME		t_gv_h_old
#		define _GV_TYPE				real (kind=GRID_SR)
#		define _GV_NAME				h_old
#		define _GV_PERSISTENT		1

#		include "Tools_grid_variable.f90"
	END MODULE

#	define _LFS_type_NAME			t_lfs_flux
#	define _LFS_CELL_SIZE			0
#	define _LFS_EDGE_SIZE			_SWE_EDGE_SIZE
#	define _LFS_NODE_SIZE			0

	MODULE SWE_lfs_flux
		use SFC_data_types

#		define _LFS_type			real(kind = GRID_SR)
#		include "Tools_local_function_space.f90"
	END MODULE



#endif
