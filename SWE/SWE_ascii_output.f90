#include "Compilation_control.f90"

#if defined(_SWE)
	MODULE SWE_ascii_output
		use ascii_output
        use SFC_edge_traversal
		use Samoa_swe
		use SWE_euler_timestep

		implicit none

        type num_traversal_data
            double precision                :: max_water, min_water, avg_water
            integer                         :: cells
        end type


		type(t_gv_Q)												:: gv_Q
        type(ascy)                                                  :: ascii

#		define _GT_NAME								t_swe_ascii_output_traversal

#		define _GT_EDGES

#		define _GT_PRE_TRAVERSAL_OP					pre_traversal_op
#		define _GT_POST_TRAVERSAL_OP				post_traversal_op
#		define _GT_PRE_TRAVERSAL_GRID_OP				pre_traversal_grid_op
#		define _GT_POST_TRAVERSAL_GRID_OP				post_traversal_grid_op

#		define _GT_ELEMENT_OP						element_op

#		define _GT_CELL_TO_EDGE_OP				    cell_to_edge_op

#		include "SFC_generic_traversal_ringbuffer.f90"

        subroutine pre_traversal_grid_op(traversal, grid)
			type(t_swe_ascii_output_traversal), intent(inout)				:: traversal
			type(t_grid), intent(inout)							        :: grid

			ascii = create_ascy(60, 30, dble(0.001), dble([cfg%scaling, cfg%scaling]), dble(cfg%offset))                                            
            
		end subroutine

        subroutine post_traversal_grid_op(traversal, grid)
			type(t_swe_ascii_output_traversal), intent(inout)				:: traversal
			type(t_grid), intent(inout)							        :: grid

            character (len = 64)							:: s_file_name
			integer(4)										:: i_rank, i_section, e_io
			logical                                         :: l_exists

#           if defined(_MPI)
                call mpi_barrier(MPI_COMM_WORLD)
#           endif
            
            call reduce(traversal%max_water, traversal%children%max_water, MPI_MAX, .true.)
            call reduce(traversal%min_water, traversal%children%min_water, MPI_MIN, .true.)
            call reduce(traversal%avg_water, traversal%children%avg_water, MPI_SUM, .true.) !-----------------
            call reduce(traversal%cells, traversal%children%cells, MPI_SUM, .true.) 
            traversal%avg_water = traversal%avg_water / traversal%cells
            
            ascii%h_min = traversal%min_water
            ascii%h_max = traversal%max_water
            ascii%h_avg = traversal%avg_water
            
            call print_it(ascii)

        end subroutine

		subroutine pre_traversal_op(traversal, section)
			type(t_swe_ascii_output_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)							:: section
            
            type(t_section_info)                                        :: info
            
            info = section%get_capacity()
            
            traversal%min_water = huge(1.0_GRID_SR)
            traversal%max_water = -huge(1.0_GRID_SR)
            traversal%avg_water = 0
            traversal%cells = info%i_cells
            
		end subroutine

		subroutine post_traversal_op(traversal, section)
			type(t_swe_ascii_output_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)							:: section
		end subroutine


		!******************
		!Geometry operators
		!******************

		subroutine element_op(traversal, section, element)
			type(t_swe_ascii_output_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)							:: section
			type(t_element_base), intent(inout)					:: element

			!local variables

			type(t_state), dimension(_SWE_CELL_SIZE)			:: Q
            double precision, dimension(2)                      :: coords1, coords2, coords3
            double precision :: h,b
            
            ! height, bathymetry
            h = dble(t_basis_Q_eval([1.0_GRID_SR/3.0_GRID_SR, 1.0_GRID_SR/3.0_GRID_SR], Q%h))
            b = dble(t_basis_Q_eval([1.0_GRID_SR/3.0_GRID_SR, 1.0_GRID_SR/3.0_GRID_SR], Q%b))
            
            !coordinates of the cell corners
            coords1 = dble(samoa_barycentric_to_world_point(element%transform_data, [1.0_GRID_SR, 0.0_GRID_SR]))
            coords2 = dble(samoa_barycentric_to_world_point(element%transform_data, [0.0_GRID_SR, 0.0_GRID_SR]))
            coords3 = dble(samoa_barycentric_to_world_point(element%transform_data, [0.0_GRID_SR, 1.0_GRID_SR]))
            
			call gv_Q%read(element, Q)

            call fill_sao(ascii, coords1, coords2, coords3, h, b, traversal%min_water, traversal%max_water, traversal%avg_water) 
            traversal%min_water = min(h, traversal%min_water)
            traversal%max_water = max(h, traversal%max_water)
            traversal%avg_water = traversal%avg_water + h
            
		end subroutine
	END MODULE
#endif
