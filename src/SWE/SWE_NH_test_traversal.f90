

#include "Compilation_control.f90"

#if defined (_SWE)
        MODULE SWE_nh_test_traversal
            use SFC_edge_traversal
            use Samoa_swe

        type num_traversal_data
        end type

        type(swe_gv_div) :: gv_div
        type(swe_gv_rhs) :: gv_rhs
        type(swe_gv_r) :: gv_r
        type(swe_gv_qp) ::gv_qp
        type(swe_gm_A):: gm_A

#define _GT_NAME t_swe_nh_test_traversal

#define _GT_NODES
#define _GT_EDGES
#define _GT_EDGES_TEMP
#define _GT_ELEMENT_OP element_op
#define	_GT_INNER_ELEMENT_OP inner_element_op
#define _GT_NODE_FIRST_TOUCH_OP node_first_touch_op


#include "SFC_generic_traversal_ringbuffer.f90"


!*******************************
!Geometry operators
!*******************************

elemental subroutine node_first_touch_op(traversal, section, node)
    type(t_swe_NH_test_traversal), intent(in) :: traversal
    type(t_grid_section), intent(in) :: section
    type(t_node_data), intent(inout) :: node
    node%data_pers%div_old=node%data_pers%div
    node%data_pers%div = 0.0_GRID_SR
    node%data_pers%r=node%data_pers%rhs
end subroutine


subroutine element_op(traversal, section, element)
    type(t_swe_NH_test_traversal), intent(inout) :: traversal
    type(t_grid_section), intent(inout) :: section
    type(t_element_base), intent(inout), target	:: element

     real (kind=GRID_SR),dimension(2)::  midpoint12, midpoint13, midpoint23
    real (kind=GRID_SR)             ::hu,hv,c, half_hypo
    real (kind = GRID_SR) :: div(3)

    call inner_element_op(traversal,section,element)

    c=0.5_GRID_SR * element%cell%geometry%get_leg_size()*cfg%scaling
    half_hypo= 0.5_GRID_SR* sqrt(8._GRID_SR)*c

    hu= element%cell%data_pers%Q(1)%p(1)
    hv= element%cell%data_pers%Q(1)%p(2)

    div=0.0_GRID_SR

    midpoint12= (element%nodes(1)%ptr%position +element%nodes(2)%ptr%position) *0.5_GRID_SR
    midpoint23= (element%nodes(2)%ptr%position +element%nodes(3)%ptr%position) *0.5_GRID_SR
    midpoint13= (element%nodes(1)%ptr%position +element%nodes(3)%ptr%position) *0.5_GRID_SR

    !leg on left boundary
!    if (midpoint12(1) .eq. 0.0_GRID_SR) then
!        div(1)=div(1)+ (-1.0_GRID_SR)*c*hu
!        div(2)=div(2)+ (-1.0_GRID_SR)*c*hu
!    !leg on right boundary
!    elseif (midpoint12(1) .eq. 1.0_GRID_SR) then
!        div(1)=div(1) + (1.0_GRID_SR)*c*hu
!        div(2)=div(2) + (1.0_GRID_SR)*c*hu
!    elseif (midpoint12(2) .eq. 0.0_GRID_SR) then
!        div(1)=div(1) + (-1.0_GRID_SR)*c*hv
!        div(2)=div(2) + (-1.0_GRID_SR)*c*hv
!    elseif (midpoint12(2) .eq. 1.0_GRID_SR) then
!        div(1)=div(1) + (1.0_GRID_SR)*c*hv
!        div(2)=div(2) + (1.0_GRID_SR)*c*hv
!    endif
!
!    !leg on left boundary
!    if (midpoint23(1) .eq. 0.0_GRID_SR) then
!        div(2)=div(2)+ (-1.0_GRID_SR)*c*hu
!        div(3)=div(3)+ (-1.0_GRID_SR)*c*hu
!    !leg on right boundary
!    elseif (midpoint23(1) .eq. 1.0_GRID_SR) then
!        div(2)=div(2) + (1.0_GRID_SR)*c*hu
!        div(3)=div(3) + (1.0_GRID_SR)*c*hu
!    elseif (midpoint23(2) .eq. 0.0_GRID_SR) then
!        div(2)=div(2) + (-1.0_GRID_SR)*c*hv
!        div(3)=div(3) + (-1.0_GRID_SR)*c*hv
!    elseif (midpoint23(2) .eq. 1.0_GRID_SR) then
!        div(2)=div(2) + (1.0_GRID_SR)*c*hv
!        div(3)=div(3) + (1.0_GRID_SR)*c*hv
!    endif
!
!
!
!    !leg on left boundary
!    if (midpoint13(1) .eq. 0.0_GRID_SR) then
!        div(1)=div(1)+ (-1.0_GRID_SR)*half_hypo*hu
!        div(3)=div(3)+ (-1.0_GRID_SR)*half_hypo*hu
!    !leg on right boundary
!    elseif (midpoint13(1) .eq. 1.0_GRID_SR) then
!        div(1)=div(1) + (1.0_GRID_SR)*half_hypo*hu
!        div(3)=div(3) + (1.0_GRID_SR)*half_hypo*hu
!    elseif (midpoint13(2) .eq. 0.0_GRID_SR) then
!        div(1)=div(1) + (-1.0_GRID_SR)*half_hypo*hv
!        div(3)=div(3) + (-1.0_GRID_SR)*half_hypo*hv
!    elseif (midpoint13(2) .eq. 1.0_GRID_SR) then
!        div(1)=div(1) + (1.0_GRID_SR)*half_hypo*hv
!        div(3)=div(3) + (1.0_GRID_SR)*half_hypo*hv
!    endif

    call gv_div%add_to_element(element, div)

end subroutine

subroutine inner_element_op(traversal, section, element)
    type(t_swe_NH_test_traversal), intent(inout) :: traversal
    type(t_grid_section), intent(inout) :: section
    type(t_element_base), intent(inout), target	:: element

    real (kind=GRID_SR) :: c
    real(kind=GRID_SR):: hu,hv,w1,w2,w3
    real (kind=GRID_SR),dimension(2):: normal_x, normal_y
    real (kind=GRID_SR),dimension(3,3)::A
    real (kind=GRID_SR),dimension(3)::qp,r

    real (kind = GRID_SR):: nxvn, nxvp, nxhn, nxhp,nyvp,nyvn, nyhn, nyhp, cont1_p,cont1_w, cont2_p,cont2_w, cont3_p,cont3_w

    normal_x(1)=1.0_GRID_SR
    normal_x(2)=0.0_GRID_SR

    normal_y(1)=0.0_GRID_SR
    normal_y(2)=1.0_GRID_SR

    c=0.5_GRID_SR * element%cell%geometry%get_leg_size()*cfg%scaling

    hu= element%cell%data_pers%Q(1)%p(1)
    hv= element%cell%data_pers%Q(1)%p(2)

    w1= element%nodes(1)%ptr%data_pers%w(1)
    w2= element%nodes(2)%ptr%data_pers%w(1)
    w3= element%nodes(3)%ptr%data_pers%w(1)

    normal_x(1:2) = samoa_barycentric_to_world_normal(element%transform_data, normal_x(1:2))
    normal_x(1:2)= normal_x(1:2)* element%transform_data%custom_data%scaling * sqrt(abs(element%transform_data%plotter_data%det_jacobian))

    normal_y(1:2) = samoa_barycentric_to_world_normal(element%transform_data, normal_y(1:2))
    normal_y(1:2)= normal_y(1:2)* element%transform_data%custom_data%scaling * sqrt(abs(element%transform_data%plotter_data%det_jacobian))


    !write (*,*) 'new element'
    !write (*,*) 'node 1: ' ,element%nodes(1)%ptr%position(1) ,','  ,element%nodes(1)%ptr%position(2)

    !write (*,*) 'node 2: ' ,element%nodes(2)%ptr%position(1) ,','  ,element%nodes(2)%ptr%position(2)

    !write (*,*) 'node 3: ' ,element%nodes(3)%ptr%position(1) ,','  ,element%nodes(3)%ptr%position(2)

    !write (*,*) 'normal_x rotated: ',normal_x
    !write (*,*) 'normal_y rotated: ',normal_y

    nxvp=normal_x(1)
    nyvp=normal_x(2)

    nxhp=normal_y(1)
    nyhp=normal_y(2)

    nxvn=-nxvp
    nyvn=-nyvp
    nxhn=-nxhp
    nyhn=-nyhp

    !contribution node 1 -> triangle -> only one boundary integral, normal
    cont1_p=c*(nxvn*hu+nyvn*hv)
    cont1_w=c*c* ((w1/3.0_GRID_SR) +(w2/12.0_GRID_SR) +(w3/12.0_GRID_SR))

    cont2_p=c*((nxvp+nxhp)*hu+(nyvp+nyhp)*hv)
    cont2_w=c*c*((w1/4.0_GRID_SR) +(w2/2.0_GRID_SR) +(w3/4.0_GRID_SR))

    cont3_p=c*(nxhn*hu+nyhn*hv)
    cont3_w=c*c* ((w1/12.0_GRID_SR) +(w2/12.0_GRID_SR) +(w3/3.0_GRID_SR))

    element%nodes(1)%ptr%data_pers%div(1)=element%nodes(1)%ptr%data_pers%div(1)+ cont1_p  + cont1_w
    element%nodes(2)%ptr%data_pers%div(1)=element%nodes(2)%ptr%data_pers%div(1)+ cont2_p + cont2_w
    element%nodes(3)%ptr%data_pers%div(1)=element%nodes(3)%ptr%data_pers%div(1)+ cont3_p + cont3_w

    !compute residual
    call gv_qp%read(element, qp)
    call gm_A%read(element, A)
    r = -matmul(A, qp)

    call gv_r%add(element, r)

    end subroutine
end MODULE
#endif
