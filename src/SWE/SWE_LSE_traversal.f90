
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

        type(swe_gv_is_dirichlet_boundary):: gv_is_dirichlet

#define _GT_NAME t_swe_lse_traversal

#define _GT_NODES
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
    node%data_pers%is_dirichlet_boundary = .false.
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
    real(kind=GRID_SR):: h,hu,hv,w
    c=0.5_GRID_SR * element%cell%geometry%get_leg_size() !TODO
    dt=section%r_dt !TODO

    write (*,*) 'element_op'

    h= element%cell%data_pers%Q(1)%h
    hu= element%cell%data_pers%Q(1)%p(1)
    hv= element%cell%data_pers%Q(1)%p(2)
    w= element%cell%data_pers%Q(1)%w


    mat= reshape([c*0.5_GRID_SR*dt*h*h+2*dt*c*c*(1._GRID_SR/3._GRID_SR), -c*dt*0.5_GRID_SR*h*h*dt*0.5_GRID_SR*c*c, (1._GRID_SR/6._GRID_SR)*dt*c*c, &
            -c*h*h*dt*0.5_GRID_SR+2*dt*c*c*(1._GRID_SR/12._GRID_SR), c*dt*0.5_GRID_SR*h*h+c*dt*0.5_GRID_SR*h*h+dt*c*c, -c*dt*0.5_GRID_SR*h*h+2*dt*c*c*(1._GRID_SR/12._GRID_SR),&
            (1._GRID_SR/6._GRID_SR)*dt*c*c,-c*dt*0.5_GRID_SR*h*h+2*dt*c*c*0.25_GRID_SR, c*dt*0.5_GRID_SR*h*h+2*dt*c*c*(1._GRID_SR/3._GRID_SR) &
            ],[3,3])

    rhs= [-c*hu+ 0.5_GRID_SR*c*c*w,c*hv+c*hu+c*c*w,-c*hv+0.5_GRID_SR*c*c*w]



    call gm_a%write(element, mat)
        !call element operator
        !call alpha_volume_op(element, saturation, p, rhs, is_dirichlet, element%cell%data_pers%base_permeability)
    call gv_rhs%add_to_element(element, rhs)
end subroutine
end MODULE
#endif
