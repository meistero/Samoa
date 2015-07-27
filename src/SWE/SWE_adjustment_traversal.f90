
#include "Compilation_control.f90"

#if defined (_SWE)
MODULE SWE_adjustment_traversal
    use SFC_edge_traversal
    use Samoa_swe

type num_traversal_data
end type

type(swe_gv_rhs) :: gv_rhs
type(swe_gv_qp)  :: gv_qp
type(swe_gv_r)  :: gv_r

#define _GT_NAME t_swe_adjustment_traversal

#define _GT_NODES
#define  _GT_EDGES
#define _GT_EDGES_TEMP
!#define	_GT_INNER_ELEMENT_OP	inner_element_op
#define _GT_NODE_MERGE_OP node_merge_op

#include "SFC_generic_traversal_ringbuffer.f90"

subroutine node_merge_op(local_node, neighbor_node)
    type(t_node_data), intent(inout)			    :: local_node
    type(t_node_data), intent(in)				    :: neighbor_node

    real (kind = GRID_SR) :: rhs(1)
    real (kind = GRID_SR) :: qp(1)
    real (kind = GRID_SR) :: r(1)

    if(neighbor_node%owned_globally) then
        call gv_rhs%read(neighbor_node, rhs)
        call gv_rhs%write(local_node, rhs)

        call gv_qp%read(neighbor_node, qp)
        call gv_qp%write(local_node, qp)

        call gv_r%read(neighbor_node, r)
        call gv_r%write(local_node, r)
    end if

end subroutine

end MODULE

#endif
