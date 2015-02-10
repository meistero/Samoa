

#include "Compilation_control.f90"

#if defined (_SWE)
        MODULE SWE_nh_test_traversal
            use SFC_edge_traversal
            use Samoa_swe

        type num_traversal_data
        end type


#define _GT_NAME t_swe_nh_test_traversal

#define _GT_NODES
#define _GT_EDGES
#define _GT_EDGES_TEMP
#define _GT_ELEMENT_OP element_op
#define _GT_NODE_FIRST_TOUCH_OP node_first_touch_op
#define _GT_NODE_LAST_TOUCH_OP node_last_touch_op

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
end subroutine

elemental subroutine node_last_touch_op(traversal,section, node)
    type(t_swe_NH_test_traversal), intent(in) :: traversal
    type(t_grid_section), intent(in) :: section
    type(t_node_data), intent(inout) :: node


end subroutine

subroutine element_op(traversal, section, element)
    type(t_swe_NH_test_traversal), intent(inout) :: traversal
    type(t_grid_section), intent(inout) :: section
    type(t_element_base), intent(inout), target	:: element

    real (kind=GRID_SR) :: c
    real(kind=GRID_SR):: hu,hv,w
    real (kind=GRID_SR),dimension(2):: normal_x, normal_y
    real (kind = GRID_SR):: nxvn, nxvp, nxhn, nxhp,nyvp,nyvn, nyhn, nyhp, cont1_p,cont1_w, cont2_p,cont2_w, cont3_p,cont3_w

    !assert that q values are ok

    normal_x(1)=1
    normal_x(2)=0

    normal_y(1)=0
    normal_y(2)=1

    c=0.5_GRID_SR * element%cell%geometry%get_leg_size()

    hu= element%cell%data_pers%Q(1)%p(1)
    hv= element%cell%data_pers%Q(1)%p(2)
    w= element%cell%data_pers%Q(1)%w

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
    cont1_w=0.5_GRID_SR*c*c*w

    cont2_p=c*((nxvp+nxhp)*hu+(nyvp+nyhp)*hv)
    cont2_w=c*c*w

    cont3_p=c*(nxhn*hu+nyhn*hv)
    cont3_w=0.5_GRID_SR*c*c*w

    element%nodes(1)%ptr%data_pers%div(1)=element%nodes(1)%ptr%data_pers%r(1)+ cont1_p + cont1_w
    element%nodes(2)%ptr%data_pers%div(1)=element%nodes(2)%ptr%data_pers%r(1)+ cont2_p + cont2_w
    element%nodes(3)%ptr%data_pers%div(1)=element%nodes(3)%ptr%data_pers%r(1)+ cont3_p + cont3_w

    end subroutine
end MODULE
#endif
