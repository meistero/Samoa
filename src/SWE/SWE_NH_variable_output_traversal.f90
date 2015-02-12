

#include "Compilation_control.f90"

#if defined (_SWE)
        MODULE SWE_nh_variable_output_traversal
            use SFC_edge_traversal
            use Samoa_swe

        type num_traversal_data
        end type


#define _GT_NAME t_swe_nh_variable_output_traversal

#define _GT_NODES
#define _GT_EDGES
#define _GT_EDGES_TEMP
#define _GT_ELEMENT_OP element_op

#include "SFC_generic_traversal_ringbuffer.f90"


!*******************************
!Geometry operators
!*******************************

subroutine element_op(traversal, section, element)
    type(t_swe_nh_variable_output_traversal), intent(inout) :: traversal
    type(t_grid_section), intent(inout) :: section
    type(t_element_base), intent(inout), target	:: element

    real(kind=GRID_SR):: h,hu,hv,w,q1,q2,q3,rhs1,rhs2,rhs3

    q1= element%nodes(1)%ptr%data_pers%qp(1)
    q2= element%nodes(2)%ptr%data_pers%qp(1)
    q3= element%nodes(3)%ptr%data_pers%qp(1)

    rhs1= element%nodes(1)%ptr%data_pers%rhs(1)
    rhs2= element%nodes(2)%ptr%data_pers%rhs(1)
    rhs3= element%nodes(3)%ptr%data_pers%rhs(1)


    h= element%cell%data_pers%Q(1)%h
    hu= element%cell%data_pers%Q(1)%p(1)
    hv= element%cell%data_pers%Q(1)%p(2)
    w= element%cell%data_pers%Q(1)%w


    write (*,*) 'node 1: ' ,element%nodes(1)%ptr%position(1) ,','  ,element%nodes(1)%ptr%position(2)
    write (*,*) 'q:'    , q1
    write (*,*) 'rhs1:'    , rhs1

    write (*,*) 'node 2: ' ,element%nodes(2)%ptr%position(1) ,','  ,element%nodes(2)%ptr%position(2)
    write (*,*) 'q:'    , q2
    write (*,*) 'rhs2:'    , rhs2


    write (*,*) 'node 3: ' ,element%nodes(3)%ptr%position(1) ,','  ,element%nodes(3)%ptr%position(2)
    write (*,*) 'q:'    , q3
    write (*,*) 'rhs3:'    , rhs3


    write(*,*) 'element variables:'
    write(*,*) 'h:', h
    write(*,*) 'hu:', hu
    write(*,*) 'hv:', hv
    write(*,*) 'w:', w


    end subroutine
end MODULE
#endif


