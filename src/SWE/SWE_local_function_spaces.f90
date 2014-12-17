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

            generic:: read => read_from_element
            generic:: write=> write_to_element
end type

contains
         pure subroutine read_from_element(gm_A, element, mat)
            class(swe_gm_A), intent(in)	    :: gm_A
            type(t_element_base), intent(in)    :: element
            real (kind = GRID_SR), intent(out)  :: mat(3, 3)


        mat=element%cell%data_pers%A
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
!#define _GV_TYPE_NAME		swe_gm_A
!#define _GV_TYPE			real (kind = GRID_SR), dimension(3,3)
!#define _GV_NAME			A
!#define _GV_PERSISTENT		1

!#include "Tools_grid_variable.f90"
    END MODULE
        !type swe_gm_A
        !contains
         !   procedure, pass :: read_from_element

          !  generic:: read => read_from_element
        !end type

        !contains
         ! pure subroutine read_from_element(gm_A, element, mat)
          !  class(swe_gm_A), intent(in)	    :: gm_A
           ! type(t_element_base), intent(in)    :: element
            !real (kind = GRID_SR), intent(out)  :: mat(3, 3)

            !real (kind=GRID_SR) :: c,dt
            !real(kind=GRID_SR):: h,hu,hv,w
            !c=1 !TODO
            !dt=1 !TODO

            !h= element%cell%data_pers%Q(1)%h
            !hu= element%cell%data_pers%Q(1)%p(1)
            !hv= element%cell%data_pers%Q(1)%p(2)
            !w= element%cell%data_pers%Q(1)%w


            !mat= reshape([c*0.5_GRID_SR*dt*h*h+2*dt*c*c*(1._GRID_SR/3._GRID_SR), -c*dt*0.5_GRID_SR*h*h*dt*0.5_GRID_SR*c*c, (1._GRID_SR/6._GRID_SR)*dt*c*c, &
            !-c*h*h*dt*0.5_GRID_SR+2*dt*c*c*(1._GRID_SR/12._GRID_SR), c*dt*0.5_GRID_SR*h*h+c*dt*0.5_GRID_SR*h*h+dt*c*c, -c*dt*0.5_GRID_SR*h*h+2*dt*c*c*(1._GRID_SR/12._GRID_SR),&
            !(1._GRID_SR/6._GRID_SR)*dt*c*c,-c*dt*0.5_GRID_SR*h*h+2*dt*c*c*0.25_GRID_SR, c*dt*0.5_GRID_SR*h*h+2*dt*c*c*(1._GRID_SR/3._GRID_SR) &
            !],[3,3])
            !real (kind = GRID_SR), parameter    :: mat_const(_DARCY_P_SIZE, _DARCY_P_SIZE) = &
            !    reshape([ 0.5_GRID_SR, -0.5_GRID_SR, 0.0_GRID_SR, -0.5_GRID_SR, 1.0_GRID_SR, -0.5_GRID_SR, 0.0_GRID_SR, -0.5_GRID_SR, 0.5_GRID_SR], [_DARCY_P_SIZE, _DARCY_P_SIZE])

            !mat = element%cell%data_pers%permeability * mat_const

        !end subroutine



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
