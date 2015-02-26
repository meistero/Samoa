#include "Compilation_control.f90"

#if defined (_SWE)
        MODULE SWE_NH_traversal
            use SFC_edge_traversal
            use Samoa_swe

        type num_traversal_data
        end type

        type(swe_gv_qp)  :: gv_qp
        type(swe_gv_r)   :: gv_r

#define _GT_NAME t_swe_nh_traversal

#define _GT_NODES
#define  _GT_EDGES
#define _GT_EDGES_TEMP
#define _GT_ELEMENT_OP element_op
#define _GT_NODE_FIRST_TOUCH_OP node_first_touch_op


#include "SFC_generic_traversal_ringbuffer.f90"

elemental subroutine node_first_touch_op(traversal, section, node)
    type(t_swe_nh_traversal), intent(in) :: traversal
    type(t_grid_section), intent(in) :: section
    type(t_node_data), intent(inout) :: node

    real (kind=GRID_SR)             ::dt

    dt=section%r_dt

    if(.not. (node%data_pers%is_dirichlet_boundary(1))) then
    node%data_pers%w= node%data_pers%w + 2.0_GRID_SR*dt*node%data_pers%qp
    endif


end subroutine



subroutine element_op(traversal, section, element)
    type(t_swe_nh_traversal), intent(inout) :: traversal
    type(t_grid_section), intent(inout) :: section
    type(t_element_base), intent(inout), target	:: element


real (kind=GRID_SR),dimension(2):: p
real (kind=GRID_SR)             :: hu,hv,dt,h, q1,q2, q3,c,m11,m12,m21,m22,s

m11=element%transform_data%plotter_data%jacobian_inv(1,1)
m12=element%transform_data%plotter_data%jacobian_inv(2,1)
m21=element%transform_data%plotter_data%jacobian_inv(1,2)
m22=element%transform_data%plotter_data%jacobian_inv(2,2)
s=sqrt(abs(element%transform_data%plotter_data%det_jacobian))

c=0.5_GRID_SR * element%cell%geometry%get_leg_size()*cfg%scaling
p=element%cell%data_pers%Q(1)%p

h=element%cell%data_pers%Q(1)%h -element%cell%data_pers%Q(1)%b
hu=p(1)
hv=p(2)


dt=section%r_dt
!is_dirichlet=element%cell%data_pers%is_dirichlet

q1= element%nodes(1)%ptr%data_pers%qp(1)
q2= element%nodes(2)%ptr%data_pers%qp(1)
q3= element%nodes(3)%ptr%data_pers%qp(1)


!write (*,*) 'q1: ', q1, 'q2: ',q2 , 'q3: ', q3
!write (*,*) 'hu,hv,w1,w2,w3 before correction: ', hu,',',hv, ',',w1,',',w2,',',w3

!correction q =q2 + (x/(2*c)) * (q1 - q2) + (y/(2*c))  * (q3 - q2)
if( .not. (element%nodes(1)%ptr%data_pers%is_dirichlet_boundary(1) .or. element%nodes(2)%ptr%data_pers%is_dirichlet_boundary(1) .or. element%nodes(3)%ptr%data_pers%is_dirichlet_boundary(1)) .and. h>cfg%dry_tolerance) then
    hu= hu - 0.5_GRID_SR* dt*(h*h*s* (m11 *(q1-q2)/(2.0_GRID_SR*c) + m12*((q3-q2)/(2.0_GRID_SR*c))))
    hv= hv - 0.5_GRID_SR* dt*(h*h*s* (m21 *(q1-q2)/(2.0_GRID_SR*c) + m22*((q3-q2)/(2.0_GRID_SR*c))))
else if(h>0) then
        write(*,*) 'omitting correction in element: '
        write (*,*) 'node 1: ' ,element%nodes(1)%ptr%position(1)*cfg%scaling ,','  ,element%nodes(1)%ptr%position(2)*cfg%scaling

    write (*,*) 'node 2: ' ,element%nodes(2)%ptr%position(1)*cfg%scaling ,','  ,element%nodes(2)%ptr%position(2)*cfg%scaling

    write (*,*) 'node 3: ' ,element%nodes(3)%ptr%position(1)*cfg%scaling,','  ,element%nodes(3)%ptr%position(2)*cfg%scaling
endif
!write (*,*) 'hu,hv,w1,w2,w3 after correction: ', hu,',',hv, ',' ,w1,',',w2,',',w3

!p_local(1)=hu
!p_local(2)=hv

!p(1:2) = samoa_barycentric_to_world_normal(element%transform_data, p_local(1:2))
!p(1:2) = p_local(1:2) *(element%transform_data%custom_data%scaling * sqrt(abs(element%transform_data%plotter_data%det_jacobian)))

p(1)=hu
p(2)=hv
element%cell%data_pers%Q(1)%p=p



end subroutine
end module
#endif
