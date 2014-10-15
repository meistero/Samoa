#include "Compilation_control.f90"

#if defined(_SWE)
	MODULE SWE_ascii_output
		use ascii_output
        use SFC_edge_traversal
		use Samoa_swe
		use SWE_euler_timestep

		implicit none

        type num_traversal_data
            real (kind = GRID_SR)           :: max_water, min_water, avg_water
            integer                         :: cells
        end type


	type(t_gv_Q)				:: gv_Q
        type(data_for_ascii)                              :: ascii

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

            try (cfg%i_ascii_width >= 2, "Invalid ascii output width")
            ascii = create_data_for_ascii(cfg%i_ascii_width, (cfg%i_ascii_width/2), 0.001_GRID_SR, [1.0_GRID_SR, -1.0_GRID_SR], [0.0_GRID_SR, 1.0_GRID_SR])

	end subroutine

        subroutine post_traversal_grid_op(traversal, grid)
			type(t_swe_ascii_output_traversal), intent(inout)				:: traversal
			type(t_grid), intent(inout)							        :: grid

        integer                                         :: i_error, alloc_err, i, j, k
        character (len = 64)				:: s_file_name
	integer(4)					:: i_rank, i_section, e_io
	logical                                         :: l_exists
	type(bath_height_data), pointer				:: dummy_ascii(:,:,:) => null() !groesse/dimension auf null setzen
	type(bath_height_data), pointer 				:: big_matrix_array(:,:,:) => null()

#           if defined(_MPI)
                call mpi_barrier(MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
#           endif

            call reduce(traversal%max_water, traversal%children%max_water, MPI_MAX, .true.)
            call reduce(traversal%min_water, traversal%children%min_water, MPI_MIN, .true.)
            call reduce(traversal%avg_water, traversal%children%avg_water, MPI_SUM, .true.) !-----------------
            call reduce(traversal%cells, traversal%children%cells, MPI_SUM, .true.)
            traversal%avg_water = traversal%avg_water / traversal%cells

            ascii%h_min = traversal%min_water
            ascii%h_max = traversal%max_water
            ascii%h_avg = traversal%avg_water

#           if defined(_MPI)

            if (rank_MPI == 0) then
                allocate (big_matrix_array(size(ascii%mat,dim=1),size(ascii%mat,dim=2),size_MPI), stat = alloc_err)
                if (alloc_err > 0) then
                    write(*,'(A)') "Error when trying to allocate rank-0-ascii"
                end if
            else
	      allocate (dummy_ascii(1,1,1), stat = alloc_err)
              if (alloc_err > 0) then
                    write(*,'(A)') "Error when trying to allocate rank-0-ascii"
              end if
	      big_matrix_array => dummy_ascii
            end if

            call MPI_gather(ascii%mat(1,1), sizeof(ascii%mat), MPI_BYTE, &
              big_matrix_array(1,1,1), sizeof(ascii%mat), MPI_BYTE, 0, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)

            if (rank_MPI == 0) then
                do i=2, (size_MPI)
                    do j=1, size(ascii%mat,dim=1)
                        do k=1, size(ascii%mat, dim=2)
                            ascii%mat(j,k)%b = ascii%mat(j,k)%b + big_matrix_array(j,k,i)%b
                            ascii%mat(j,k)%h_sum = ascii%mat(j,k)%h_sum + big_matrix_array(j,k,i)%h_sum
                            ascii%mat(j,k)%h_summands = ascii%mat(j,k)%h_summands + big_matrix_array(j,k,i)%h_summands
                        end do
                    end do
                end do

#           endif

                call print_it(ascii)

#           if defined(_MPI)

            end if

            if (rank_MPI == 0) then
                deallocate (big_matrix_array, stat = alloc_err)
                if (alloc_err > 0) then
                    write(*,'(A)') "Error when trying to deallocate rank-0-ascii"
                end if
            else
                deallocate (dummy_ascii, stat = alloc_err)
                if (alloc_err > 0) then
                    write(*,'(A)') "Error when trying to deallocate rank-X-dummy-ascii"
                end if
            end if

#	    endif

        end subroutine

		subroutine pre_traversal_op(traversal, section)
			type(t_swe_ascii_output_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)							:: section

            type(t_section_info)                                        :: info

            info = section%get_info()

            traversal%min_water = huge(1.0_GRID_SR)
            traversal%max_water = -huge(1.0_GRID_SR)
            traversal%avg_water = 0.0_GRID_SR
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
            real (kind = GRID_SR), dimension(2)                      :: coords1, coords2, coords3
            real (kind = GRID_SR) :: h,b

            call gv_Q%read(element, Q)

            ! height, bathymetry
            h = t_basis_Q_eval([1.0_GRID_SR/3.0_GRID_SR, 1.0_GRID_SR/3.0_GRID_SR], Q%h)
            b = t_basis_Q_eval([1.0_GRID_SR/3.0_GRID_SR, 1.0_GRID_SR/3.0_GRID_SR], Q%b)

            !coordinates of the cell corners
            coords1 = samoa_barycentric_to_world_point(element%transform_data, [1.0_GRID_SR, 0.0_GRID_SR])
            coords2 = samoa_barycentric_to_world_point(element%transform_data, [0.0_GRID_SR, 0.0_GRID_SR])
            coords3 = samoa_barycentric_to_world_point(element%transform_data, [0.0_GRID_SR, 1.0_GRID_SR])

            call fill_ascii(ascii, coords1, coords2, coords3, h, b, traversal%min_water, traversal%max_water, traversal%avg_water)
            traversal%min_water = min(h, traversal%min_water)
            traversal%max_water = max(h, traversal%max_water)
            traversal%avg_water = traversal%avg_water + h

		end subroutine
	END MODULE
#endif
