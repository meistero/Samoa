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
!#define _GT_NODE_FIRST_TOUCH_OP node_first_touch_op


#include "SFC_generic_traversal_ringbuffer.f90"
subroutine element_op(traversal, section, element)
    type(t_swe_nh_traversal), intent(inout) :: traversal
    type(t_grid_section), intent(inout) :: section
    type(t_element_base), intent(inout), target	:: element


real (kind=GRID_SR),dimension(2):: p, p_local, normal_x, normal_y, normal_x_r, normal_y_r
real (kind=GRID_SR)             :: hu,hv,dt,h, q1,q2, q3,c,w,area, m11,m12,m21,m22,s
!logical, dimension(3)           :: is_dirichlet

!normal_x(1)=1
!normal_x(2)=0

!normal_y(1)=0
!normal_y(2)=1

!normal_x_r(1:2) = samoa_world_to_barycentric_normal(element%transform_data, normal_x(1:2))
!normal_x_r(1:2) = normal_x_r(1:2) / (element%transform_data%custom_data%scaling * sqrt(abs(element%transform_data%plotter_data%det_jacobian)))

!normal_y_r(1:2) = samoa_world_to_barycentric_normal(element%transform_data, normal_y(1:2))
!normal_y_r(1:2) = normal_y_r(1:2) / (element%transform_data%custom_data%scaling * sqrt(abs(element%transform_data%plotter_data%det_jacobian)))


!write (*,*) 'h: ', h
!write (*,*) 'normal_x: ',normal_x, '  normal_x rotated:', normal_x_r
!write (*,*) 'normal_y: ',normal_y, '  normal_y rotated:', normal_y_r

!normal_x(1:2) = samoa_barycentric_to_world_normal(element%transform_data, normal_x_r(1:2))
!normal_x(1:2)= normal_x(1:2)* element%transform_data%custom_data%scaling * sqrt(abs(element%transform_data%plotter_data%det_jacobian))

!normal_y(1:2) = samoa_barycentric_to_world_normal(element%transform_data, normal_y_r(1:2))
!normal_y(1:2)= normal_y(1:2)* element%transform_data%custom_data%scaling * sqrt(abs(element%transform_data%plotter_data%det_jacobian))

!write (*,*) 'normal_x rotated back: ',normal_x
!write (*,*) 'normal_y rotated back: ',normal_y

m11=element%transform_data%plotter_data%jacobian_inv(1,1)
    m12=element%transform_data%plotter_data%jacobian_inv(2,1)
    m21=element%transform_data%plotter_data%jacobian_inv(1,2)
    m22=element%transform_data%plotter_data%jacobian_inv(2,2)
    s=sqrt(abs(element%transform_data%plotter_data%det_jacobian))

c=0.5_GRID_SR * element%cell%geometry%get_leg_size()
p=element%cell%data_pers%Q(1)%p
!p_local(1:2) = samoa_world_to_barycentric_normal(element%transform_data, p(1:2))
!p_local(1:2) = p_local(1:2) / (element%transform_data%custom_data%scaling * sqrt(abs(element%transform_data%plotter_data%det_jacobian)))
!p_local=p

hu=p(1)
hv=p(2)
w= element%cell%data_pers%Q(1)%w
dt=section%r_dt
!is_dirichlet=element%cell%data_pers%is_dirichlet

q1= element%nodes(1)%ptr%data_pers%qp(1)
q2= element%nodes(2)%ptr%data_pers%qp(1)
q3= element%nodes(3)%ptr%data_pers%qp(1)

!write (*,*) 'q1: ', q1, 'q2: ',q2 , 'q3: ', q3
!write (*,*) 'hu,hv,w before correction: ', hu,',',hv, ',',w

!correction q =q2 + (x/(2*c)) * (q1 - q2) + (y/(2*c))  * (q3 - q2)
area= (2.0_GRID_SR*c)*(2.0_GRID_SR*c) *0.5_GRID_SR
hu= hu - 0.5_GRID_SR* dt*(h*h*s* (m11 *(q1-q2)/(2*c) + m12*((q3-q2)/(2*c))))
hv= hv - 0.5_GRID_SR* dt*(h*h*s* (m21 *(q1-q2)/(2*c) + m22*((q3-q2)/(2*c))))
w =w + (1._GRID_SR/area) * 2._GRID_SR* dt* (1._GRID_SR/3._GRID_SR)* (2*c*c*(q1+q2+q3))

!write (*,*) 'hu,hv,w after correction: ', hu,',',hv, ',' ,w

!p_local(1)=hu
!p_local(2)=hv

!p(1:2) = samoa_barycentric_to_world_normal(element%transform_data, p_local(1:2))
!p(1:2) = p_local(1:2) *(element%transform_data%custom_data%scaling * sqrt(abs(element%transform_data%plotter_data%det_jacobian)))

p(1)=hu
p(2)=hv
element%cell%data_pers%Q(1)%p=p
element%cell%data_pers%Q(1)%w=w

end subroutine
end module
#endif
