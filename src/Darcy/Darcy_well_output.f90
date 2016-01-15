
! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_DARCY)
	MODULE Darcy_well_output
        use SFC_edge_traversal
        use Tools_log

        type :: t_darcy_well_output_traversal
            character(len=256)			:: s_file_stamp
            integer (kind = GRID_SI)    :: i_output_iteration = 0
            contains

            procedure, pass :: create
            procedure, pass :: destroy
            procedure, pass :: traverse
        end type

        contains

        subroutine create(traversal)
            class(t_darcy_well_output_traversal), intent(inout)      :: traversal
        end subroutine

        subroutine destroy(traversal)
            class(t_darcy_well_output_traversal), intent(inout)      :: traversal

        end subroutine

		subroutine traverse(traversal, grid)
            class(t_darcy_well_output_traversal), intent(inout)		:: traversal
			type(t_grid), intent(inout)							    :: grid

			character (len = 256)							:: s_file_name
            integer                                         :: i_error
			integer(4)										:: i_rank, i_section, e_io
			logical                                         :: l_exists

            real (kind = GRID_SR)   :: reduction_set(4 * (-_DARCY_PRODUCER_WELLS) + 1 : 4 * (_DARCY_INJECTOR_WELLS) + 4)
            real (kind = GRID_SR)   :: prod_w(-_DARCY_PRODUCER_WELLS : _DARCY_INJECTOR_WELLS)
            real (kind = GRID_SR)   :: prod_n(-_DARCY_PRODUCER_WELLS : _DARCY_INJECTOR_WELLS)
            real (kind = GRID_SR)   :: prod_w_acc(-_DARCY_PRODUCER_WELLS : _DARCY_INJECTOR_WELLS)
            real (kind = GRID_SR)   :: prod_n_acc(-_DARCY_PRODUCER_WELLS : _DARCY_INJECTOR_WELLS)
            real (kind = GRID_SR)   :: p_bh(_DARCY_INJECTOR_WELLS)
            integer                 :: i, f_out

            !$omp single
                do i = -_DARCY_PRODUCER_WELLS, _DARCY_INJECTOR_WELLS
                    reduction_set(4*i + 1) = grid%prod_w(i)
                    reduction_set(4*i + 2) = grid%prod_n(i)
                    reduction_set(4*i + 3) = grid%prod_w_acc(i)
                    reduction_set(4*i + 4) = grid%prod_n_acc(i)
                end do

                call reduce(reduction_set, MPI_SUM)

                if (rank_MPI == 0) then
                    do i = -_DARCY_PRODUCER_WELLS, _DARCY_INJECTOR_WELLS
                        prod_w(i) = reduction_set(4*i + 1)
                        prod_n(i) = reduction_set(4*i + 2)
                        prod_w_acc(i) = reduction_set(4*i + 3)
                        prod_n_acc(i) = reduction_set(4*i + 4)
                    end do

                    prod_w(0) = sum(prod_w(-_DARCY_PRODUCER_WELLS:-1))
                    prod_n(0) = sum(prod_n(-_DARCY_PRODUCER_WELLS:-1))
                    prod_w_acc(0)= sum(prod_w_acc(-_DARCY_PRODUCER_WELLS:-1))
                    prod_n_acc(0) = sum(prod_n_acc(-_DARCY_PRODUCER_WELLS:-1))

                    !we already reduce the bottom hole pressure in each iteration
                    p_bh = grid%p_bh

                    write (s_file_name, "(A, A, I0, A, I0, A)") trim(traversal%s_file_stamp), "_", traversal%i_output_iteration, ".csv"

                    f_out = get_free_file_unit()

                    open(unit=f_out, file=s_file_name, action="write", status="replace")

                    write(f_out, '("time,water rate,oil rate,water cumulative,oil cumulative,water cut,pressure")')

                    do i = 0, -_DARCY_PRODUCER_WELLS, -1
                        write(f_out, '(6(ES18.7E3, ","), ES18.7E3)') grid%r_time / _D, prod_w(i) / (_BBL / _D), prod_n(i) / (_BBL / _D), prod_w_acc(i) / _BBL, prod_n_acc(i) / _BBL, prod_w(i) / (tiny(1.0_SR) + prod_w(i) + prod_n(i)), cfg%r_p_prod / _PPSI
                    end do

                    do i = 1, _DARCY_INJECTOR_WELLS
                        write(f_out, '(6(ES18.7E3, ","), ES18.7E3)') grid%r_time / _D, prod_w(i) / (_BBL / _D), prod_n(i) / (_BBL / _D), prod_w_acc(i) / _BBL, prod_n_acc(i) / _BBL, prod_w(i) / (tiny(1.0_SR) + prod_w(i) + prod_n(i)), p_bh(i) / _PPSI
                    end do

                    close(f_out)
                end if

                traversal%i_output_iteration = traversal%i_output_iteration + 1
            !$omp end single
		end subroutine
	END MODULE
#endif
