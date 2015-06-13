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
#define _GT_NODE_LAST_TOUCH_OP node_last_touch_op
#define _GT_NODE_MERGE_OP node_merge_op

#include "SFC_generic_traversal_ringbuffer.f90"


elemental subroutine node_first_touch_op(traversal, section, node)
    type(t_swe_nh_traversal), intent(in) :: traversal
    type(t_grid_section), intent(in) :: section
    type(t_node_data), intent(inout) :: node

    real (kind=GRID_SR)             ::dt

    dt=section%r_dt

    node%data_temp%area_ctrl_vol=0.0
    node%data_temp%h_sum_inv=0.0

end subroutine

elemental subroutine node_last_touch_op(traversal, section, node)
    type(t_swe_nh_traversal), intent(in) :: traversal
    type(t_grid_section), intent(in) :: section
    type(t_node_data), intent(inout) :: node

    real (kind=GRID_SR)             ::dt

    dt=section%r_dt

    if(.not. (node%data_pers%is_dirichlet_boundary(1))) then
        node%data_pers%w= node%data_pers%w + 2.0_GRID_SR*dt*node%data_pers%qp *node%data_temp%h_sum_inv /node%data_temp%area_ctrl_vol
    endif

end subroutine



subroutine element_op(traversal, section, element)
    type(t_swe_nh_traversal), intent(inout) :: traversal
    type(t_grid_section), intent(inout) :: section
    type(t_element_base), intent(inout), target	:: element


real (kind=GRID_SR),dimension(2):: p
real (kind=GRID_SR)             :: hu,hv,dt,h,h_old, q1,q2, q3,c,m11,m12,m21,m22,s, area_tri, area_sq

m11=element%transform_data%plotter_data%jacobian_inv(1,1)
m12=element%transform_data%plotter_data%jacobian_inv(2,1)
m21=element%transform_data%plotter_data%jacobian_inv(1,2)
m22=element%transform_data%plotter_data%jacobian_inv(2,2)
s=sqrt(abs(element%transform_data%plotter_data%det_jacobian))

c=0.5_GRID_SR * element%cell%geometry%get_leg_size()*cfg%scaling
p=element%cell%data_pers%Q(1)%p

h=element%cell%data_pers%Q(1)%h -element%cell%data_pers%Q(1)%b
h_old= element%cell%data_pers%h_old(1) -  element%cell%data_pers%Q(1)%b
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

if( .not. (element%nodes(1)%ptr%data_pers%is_dirichlet_boundary(1) .or. element%nodes(2)%ptr%data_pers%is_dirichlet_boundary(1) .or. element%nodes(3)%ptr%data_pers%is_dirichlet_boundary(1)) .and. h>cfg%dry_tolerance ) then! .and. (.not.(abs(q1)>3 .or. abs(q2)>3 .or. abs(q3)>3))) then
    hu= hu - 0.5_GRID_SR* dt*(h_old*s* (m11 *(q1-q2)/(2.0_GRID_SR*c) + m12*((q3-q2)/(2.0_GRID_SR*c))))
    hv= hv - 0.5_GRID_SR* dt*(h_old*s* (m21 *(q1-q2)/(2.0_GRID_SR*c) + m22*((q3-q2)/(2.0_GRID_SR*c))))
end if
!else if(h>0) then
!        write(*,*) 'omitting correction in element: '
!        write (*,*) 'node 1: ' ,element%nodes(1)%ptr%position(1)*cfg%scaling ,','  ,element%nodes(1)%ptr%position(2)*cfg%scaling
!
!    write (*,*) 'node 2: ' ,element%nodes(2)%ptr%position(1)*cfg%scaling ,','  ,element%nodes(2)%ptr%position(2)*cfg%scaling
!
!    write (*,*) 'node 3: ' ,element%nodes(3)%ptr%position(1)*cfg%scaling,','  ,element%nodes(3)%ptr%position(2)*cfg%scaling
!endif
!write (*,*) 'hu,hv,w1,w2,w3 after correction: ', hu,',',hv, ',' ,w1,',',w2,',',w3

!p_local(1)=hu
!p_local(2)=hv

!p(1:2) = samoa_barycentric_to_world_normal(element%transform_data, p_local(1:2))
!p(1:2) = p_local(1:2) *(element%transform_data%custom_data%scaling * sqrt(abs(element%transform_data%plotter_data%det_jacobian)))
area_sq=c*c
area_tri=0.5_GRID_SR*area_sq

if(h>0) then
element%nodes(1)%ptr%data_temp%h_sum_inv=element%nodes(1)%ptr%data_temp%h_sum_inv + (1/h)*area_tri
element%nodes(2)%ptr%data_temp%h_sum_inv=element%nodes(2)%ptr%data_temp%h_sum_inv + (1/h)*area_sq
element%nodes(3)%ptr%data_temp%h_sum_inv=element%nodes(3)%ptr%data_temp%h_sum_inv + (1/h)*area_tri
endif

element%nodes(1)%ptr%data_temp%area_ctrl_vol=element%nodes(1)%ptr%data_temp%area_ctrl_vol + area_tri
element%nodes(2)%ptr%data_temp%area_ctrl_vol=element%nodes(2)%ptr%data_temp%area_ctrl_vol + area_sq
element%nodes(3)%ptr%data_temp%area_ctrl_vol=element%nodes(3)%ptr%data_temp%area_ctrl_vol + area_tri

p(1)=hu
p(2)=hv
element%cell%data_pers%Q(1)%p=p

!if( element%nodes(1)%ptr%position(1) == 0.5_GRID_SR .and. element%nodes(1)%ptr%position(2) == 0.5_GRID_SR ) then

!    open(unit=99999,file='out.txt',status='unknown',access='sequential',position='append',action='write')
!!    write(99999,*) 'rhs: ', element%nodes(1)%ptr%data_pers%rhs
!!    write(99999,*) 'w: ', element%nodes(1)%ptr%data_pers%w
!!    write(99999,*) 'qp: ', element%nodes(1)%ptr%data_pers%qp
!!    write(99999,*) 'r: ', element%nodes(1)%ptr%data_pers%r
!!    !write(99999,*) 'mat_diag: ', element%nodes(1)%ptr%data_pers%mat_diagonal(1)
!!    write(99999,*) 'is_dirichlet: ', element%nodes(1)%ptr%data_pers%is_dirichlet_boundary
!!    write(99999,*) 'A: ', element%cell%data_pers%A
!!    write(99999,*) 'Q: ', element%cell%data_pers%Q
!!    write(99999,'()')
!    write(99999,*) element%nodes(1)%ptr%position, element%nodes(1)%ptr%data_pers%rhs, element%nodes(1)%ptr%data_pers%w, element%nodes(1)%ptr%data_pers%qp, element%nodes(1)%ptr%data_pers%r, element%cell%data_pers%A
!    write(99999,*) element%nodes(2)%ptr%position, element%nodes(2)%ptr%data_pers%rhs, element%nodes(2)%ptr%data_pers%w, element%nodes(2)%ptr%data_pers%qp, element%nodes(2)%ptr%data_pers%r, element%cell%data_pers%A
!    write(99999,*) element%nodes(3)%ptr%position, element%nodes(3)%ptr%data_pers%rhs, element%nodes(3)%ptr%data_pers%w, element%nodes(3)%ptr%data_pers%qp, element%nodes(3)%ptr%data_pers%r, element%cell%data_pers%A
!    close(unit=99999)

!endif
end subroutine


subroutine node_merge_op(local_node, neighbor_node)
    type(t_node_data), intent(inout)			    :: local_node
    type(t_node_data), intent(in)				    :: neighbor_node

    assert_eqv(local_node%data_pers%is_dirichlet_boundary(1),neighbor_node%data_pers%is_dirichlet_boundary(1))
    assert_eq(local_node%data_pers%w(1), neighbor_node%data_pers%w(1))
    assert_eq(local_node%data_pers%qp(1), neighbor_node%data_pers%qp(1))
    assert_eq(local_node%data_pers%rhs(1), neighbor_node%data_pers%rhs(1))
    assert_eq(local_node%data_pers%r(1), neighbor_node%data_pers%r(1))

!    if( local_node%position(1) == 0.5_GRID_SR .and. local_node%position(2) == 0.5_GRID_SR .and. neighbor_node%position(1) == 0.5_GRID_SR .and. neighbor_node%position(2) == 0.5_GRID_SR ) then
!        write(*,*) 'rhs: ', local_node%data_pers%rhs
!        write(*,*) 'w: ', local_node%data_pers%w
!        write(*,*) 'qp: ', local_node%data_pers%qp
!        write(*,*) 'r: ', local_node%data_pers%r
!    endif

    local_node%data_temp%h_sum_inv = local_node%data_temp%h_sum_inv + neighbor_node%data_temp%h_sum_inv
    local_node%data_temp%area_ctrl_vol = local_node%data_temp%area_ctrl_vol + neighbor_node%data_temp%area_ctrl_vol
end subroutine

end module
#endif
