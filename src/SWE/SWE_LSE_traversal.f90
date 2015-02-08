
#include "Compilation_control.f90"

#if defined (_SWE)
        MODULE SWE_lse_traversal
            use SFC_edge_traversal
            use Samoa_swe

        type num_traversal_data
        end type

        type(swe_gv_rhs) :: gv_rhs
        type(swe_gv_qp)  :: gv_qp
        type(swe_gm_a)   :: gm_a
        type(swe_gv_original_lse_orientation):: gv_original_lse_orientation

        type(swe_gv_is_dirichlet_boundary):: gv_is_dirichlet

#define _GT_NAME t_swe_lse_traversal

#define _GT_NODES
#define  _GT_EDGES
#define _GT_EDGES_TEMP
#define _GT_ELEMENT_OP element_op
#define _GT_NODE_FIRST_TOUCH_OP node_first_touch_op


#include "SFC_generic_traversal_ringbuffer.f90"


        !*******************************
!Geometry operators
!*******************************
elemental subroutine node_first_touch_op(traversal, section, node)
    type(t_swe_lse_traversal), intent(in) :: traversal
    type(t_grid_section), intent(in) :: section
    type(t_node_data), intent(inout) :: node
    node%data_pers%rhs = 0.0_GRID_SR
    !this is required, otw. NaN values occur after an Euler traversal, why?
    node%data_pers%qp = 0.0_GRID_SR


    if (node%position(1) > 0.0_GRID_SR .and. node%position(1) < 1.0_GRID_SR .and. node%position(2) > 0.0_GRID_SR .and. node%position(2) < 1.0_GRID_SR ) then
            node%data_pers%is_dirichlet_boundary = .false.
    else
            node%data_pers%is_dirichlet_boundary = .true.
    end if

end subroutine


subroutine element_op(traversal, section, element)
    type(t_swe_lse_traversal), intent(inout) :: traversal
    type(t_grid_section), intent(inout) :: section
    type(t_element_base), intent(inout), target	:: element


    real (kind = GRID_SR) :: qp(3)
    real (kind = GRID_SR) :: rhs(3)
    logical	:: is_dirichlet(3)
    real (kind = GRID_SR)  :: mat(3, 3)
    real (kind=GRID_SR) :: c,dt
    real(kind=GRID_SR):: h,hu,hv,w,q1,q2,q3
    real (kind=GRID_SR),dimension(2):: p, p_local, normal_x, normal_y
    real (kind = GRID_SR):: nxvn, nxvp, nxhn, nxhp,nyvp,nyvn, nyhn, nyhp,s, m11,m12, m21,m22

    !nxvn=-1
    !nxvp=1
    !nyvp=0
    !nxvn=0

    !nxhn=0
    !nxhp=0
    !nyhn=-1
    !nyhp=1
    q1= element%nodes(1)%ptr%data_pers%qp(1)
    q2= element%nodes(2)%ptr%data_pers%qp(1)
    q3= element%nodes(3)%ptr%data_pers%qp(1)
    assert_eq(q1,q1)
    assert_eq(q2,q2)
    assert_eq(q3,q3)



    normal_x(1)=1
    normal_x(2)=0

    normal_y(1)=0
    normal_y(2)=1


    c=0.5_GRID_SR * element%cell%geometry%get_leg_size()
    dt=section%r_dt

    h= element%cell%data_pers%Q(1)%h
    hu= element%cell%data_pers%Q(1)%p(1)
    hv= element%cell%data_pers%Q(1)%p(2)
    w= element%cell%data_pers%Q(1)%w

    normal_x(1:2) = samoa_barycentric_to_world_normal(element%transform_data, normal_x(1:2))
    normal_x(1:2)= normal_x(1:2)* element%transform_data%custom_data%scaling * sqrt(abs(element%transform_data%plotter_data%det_jacobian))

    normal_y(1:2) = samoa_barycentric_to_world_normal(element%transform_data, normal_y(1:2))
    normal_y(1:2)= normal_y(1:2)* element%transform_data%custom_data%scaling * sqrt(abs(element%transform_data%plotter_data%det_jacobian))


    m11=element%transform_data%plotter_data%jacobian_inv(1,1)
    m12=element%transform_data%plotter_data%jacobian_inv(2,1)
    m21=element%transform_data%plotter_data%jacobian_inv(1,2)
    m22=element%transform_data%plotter_data%jacobian_inv(2,2)
    s=sqrt(abs(element%transform_data%plotter_data%det_jacobian))




    !write (*,*) 'h: ', h
    !write (*,*) 'new element'
    !write (*,*) 's: ',s
    !write (*,*)  'dt: ', dt
    !write (*,*) 'node 1: ' ,element%nodes(1)%ptr%position(1) ,','  ,element%nodes(1)%ptr%position(2)

    !write (*,*) 'node 2: ' ,element%nodes(2)%ptr%position(1) ,','  ,element%nodes(2)%ptr%position(2)

    !write (*,*) 'node 3: ' ,element%nodes(3)%ptr%position(1) ,','  ,element%nodes(3)%ptr%position(2)

     !write (*,*) 'Jacobi rotate gradients/normal :',element%transform_data%plotter_data%jacobian_inv

     !write (*,*)  'm11: ', m11
     !write (*,*)  'm12: ', m12
     !write (*,*)  'm21: ', m21
     !write (*,*)  'm22: ', m22

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



    !write (*,*) 'rotation pressure test: x_rotated:' , s*(m11*1), ',' , s*(m21*1)
    !write (*,*) 'rotation pressure test: y_rotated:' , s*(m12*1), ',' , s*(m22*1)

    assert_eq(nxvp, s*(m11))
    assert_eq(nyvp, s*(m21))
    assert_eq(nxhp, s*(m12))
    assert_eq(nyhp, s*(m22))

    !write (*,*) 'orientation: ' , element%transform_data%plotter_data%orientation
    !p=element%cell%data_pers%Q(1)%p
    !rotate p so it points in the right direction (no scaling!)
    !p_local(1:2) = samoa_world_to_barycentric_normal(element%transform_data, p(1:2))
    !p_local(1:2) = p_local(1:2) / (element%transform_data%custom_data%scaling * sqrt(abs(element%transform_data%plotter_data%det_jacobian)))
    !hu=p_local(1)
    !hv=p_local(2)
    !mat= reshape([c*0.5_GRID_SR*dt*h*h+2*dt*c*c*(1._GRID_SR/3._GRID_SR), -c*dt*0.5_GRID_SR*h*h+dt*0.5_GRID_SR*c*c, (1._GRID_SR/6._GRID_SR)*dt*c*c, &
         !   -c*h*h*dt*0.5_GRID_SR+2*dt*c*c*(1._GRID_SR/12._GRID_SR), c*dt*0.5_GRID_SR*h*h+c*dt*0.5_GRID_SR*h*h+dt*c*c, -c*dt*0.5_GRID_SR*h*h+2*dt*c*c*(1._GRID_SR/12._GRID_SR),&
          !  (1._GRID_SR/6._GRID_SR)*dt*c*c,-c*dt*0.5_GRID_SR*h*h+2*dt*c*c*0.25_GRID_SR, c*dt*0.5_GRID_SR*h*h+2*dt*c*c*(1._GRID_SR/3._GRID_SR) &
          !  ],[3,3])
    !write (*,*) 'matrix=',mat
    !write (*,*) ''
    !rhs= [-c*hu+ 0.5_GRID_SR*c*c*w,c*hv+c*hu+c*c*w,-c*hv+0.5_GRID_SR*c*c*w]

    !dt=1
    !h=1
    !c=1

    mat(1,1)= 2*dt*c*c* (1._GRID_SR/3._GRID_SR)- dt*0.25_GRID_SR*s*h*h*(nxvn*m11+nyvn*m21)
    mat(2,1)= 2*dt*c*c* (1._GRID_SR/12._GRID_SR)+dt*0.25_GRID_SR*s*h*h*(nxvn*(m11+m12)+nyvn*(m21+m22))
    mat(3,1)= 2*dt*c*c* (1._GRID_SR/12._GRID_SR)-dt*0.25_GRID_SR*s*h*h*(nxvn*m12+nyvn*m22)

    mat(1,2)= 0.5_GRID_SR*c*c*dt -0.25_GRID_SR*dt*h*h*s*(nxhp*m11+nyhp*m21+nxvp*m11+nyvp*m21)
    mat(2,2)= c*c*dt+ 0.25_GRID_SR *dt*h*h*s*(nxhp*(m11+m12)+nyhp*(m21+m22)+nxvp*(m11+m12)+nyvp*(m21+m22))
    mat(3,2)= 0.5_GRID_SR*c*c*dt - 0.25_GRID_SR*dt*h*h*s*(nxhp*m12+nyhp*m22+nxvp*m12+nyvp*m22)

    mat(1,3)= 2*dt*c*c* (1._GRID_SR/12._GRID_SR) - 0.25_GRID_SR*dt*h*h*s*(nxhn*m11+nyhn*m21)
    mat(2,3)= 2*dt*c*c* (1._GRID_SR/12._GRID_SR) +0.25_GRID_SR*dt*h*h*s*(nxhn*(m11+m12)+nyhn*(m21+m22))
    mat(3,3)= c*c*dt*(2._GRID_SR/3._GRID_SR) - 0.25_GRID_SR*dt*h*h*s*(nxhn*m12+nyhn*m22)

    !write (*,*) 'element diagonal', mat(1,1), ', ', mat(2,2), ' ,' , mat(3,3)
    !write (*,*) 'element matrix', mat
    !assert_gt(mat(1,1),0.01)
    !assert_gt(mat(2,2),0.01)
    !assert_gt(mat(3,3),0.01)

    rhs=[- c*c* 0.5_GRID_SR* w- c*nxvn*hu-c*nyvn*hv, -c*c*w-hu*(nxhp+nxvp)-hv*(c*(nyhp+nyvp)), - c*c* 0.5_GRID_SR* w- c*nxhn*hu-c*nyhn*hv ]

    !write (*,*) 'rhs' , rhs


    if (element%nodes(1)%ptr%data_pers%is_dirichlet_boundary(1)) then
        rhs(1)=0
    endif

    if (element%nodes(2)%ptr%data_pers%is_dirichlet_boundary(1)) then
        rhs(2)=0
    endif

    if (element%nodes(3)%ptr%data_pers%is_dirichlet_boundary(1)) then
        rhs(3)=0
    endif
    !write (*,*) 'rhs=',rhs

    !call gv_original_lse_orientation%write(element, element%transform_data%plotter_data%orientation)
    element%cell%data_pers%original_lse_orientation=element%transform_data%plotter_data%orientation
    call gm_a%write(element, mat)
    call gv_rhs%add_to_element(element, rhs)
end subroutine
end MODULE
#endif
