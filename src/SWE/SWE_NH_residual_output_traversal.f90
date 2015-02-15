

#include "Compilation_control.f90"

#if defined (_SWE)
        MODULE SWE_nh_residual_output_traversal
            use SFC_edge_traversal
            use Samoa_swe

        type num_traversal_data
        end type


#define _GT_NAME t_swe_nh_residual_output_traversal

#define _GT_NODES
#define _GT_EDGES
#define _GT_EDGES_TEMP
#define _GT_ELEMENT_OP element_op

#include "SFC_generic_traversal_ringbuffer.f90"


!*******************************
!Geometry operators
!*******************************

subroutine element_op(traversal, section, element)
    type(t_swe_nh_residual_output_traversal), intent(inout) :: traversal
    type(t_grid_section), intent(inout) :: section
    type(t_element_base), intent(inout), target	:: element



    write (*,*) 'node 1: ' ,element%nodes(1)%ptr%position(1) ,','  ,element%nodes(1)%ptr%position(2)
    write (*,*) 'divergence old:',element%nodes(1)%ptr%data_pers%div_old(1), ' divergence new: ',element%nodes(1)%ptr%data_pers%div(1)
    write (*,*) 'residual:', element%nodes(1)%ptr%data_pers%r(1)

    write (*,*) 'node 2: ' ,element%nodes(2)%ptr%position(1) ,','  ,element%nodes(2)%ptr%position(2)
    write (*,*) 'divergence old:',element%nodes(2)%ptr%data_pers%div_old(1), ' divergence new: ',element%nodes(2)%ptr%data_pers%div(1)
     write (*,*) 'residual:', element%nodes(2)%ptr%data_pers%r(1)

    write (*,*) 'node 3: ' ,element%nodes(3)%ptr%position(1) ,','  ,element%nodes(3)%ptr%position(2)
    write (*,*) 'divergence old:',element%nodes(3)%ptr%data_pers%div_old(1), ' divergence new: ',element%nodes(3)%ptr%data_pers%div(1)
     write (*,*) 'residual:', element%nodes(3)%ptr%data_pers%r(1)

    end subroutine
end MODULE
#endif

