MODULE ascii_output
    use SFC_data_types

    implicit none

    type bath_height_data
        real (kind = GRID_SR)            :: b                        ! bathymetry
        real (kind = GRID_SR)            :: h_sum                    ! sum up all water heights for the representing symbol
        integer                     :: h_summands               ! count the summed heights -> average height is computed later
    end type

    type data_for_ascii
        type(bath_height_data), allocatable    :: mat (:,:)                ! matrix to save the above values for alle symbols
        real (kind = GRID_SR)            :: eps, s(2), a(2)          ! epsilon to declare flat water; scaling and offset for coordinate transformation
        real (kind = GRID_SR)            :: h_avg, h_min, h_max      ! average, minimal and maximal water height
    end type

    public

    contains

function create_data_for_ascii(dim_x, dim_y, eps, s, a) result (ascii)
    integer             :: dim_x, dim_y
    real (kind = GRID_SR)    :: eps, s(2), a(2)

    type(data_for_ascii)          :: ascii
    integer             :: data_for_ascii_err, i, j

    ! allocate
    if (allocated(ascii%mat)) then
        deallocate(ascii%mat)
    end if

    allocate (ascii%mat(dim_x, dim_y), stat = data_for_ascii_err)
    if (data_for_ascii_err > 0) then
        write(*,'(A)') "Error when trying to allocate for ascii output"
    end if

    ! setting the base values for the whole matrix
    do i=1, size(ascii%mat,dim=2)
        do j=1, size(ascii%mat,dim=1)
            ascii%mat(j,i)%b = 0.0_GRID_SR
            ascii%mat(j,i)%h_sum = 0.0_GRID_SR
            ascii%mat(j,i)%h_summands = 0
        end do
    end do

    ! setting the other base values
    ascii%eps = eps
    ascii%s = s
    ascii%a = a
    ascii%h_min = huge(1.0_GRID_SR)
    ascii%h_max = -huge(1.0_GRID_SR)
    ascii%h_avg = 0

end function

!when traversing the grid, call the following function from inside a cell to
! transfer its data to the simplified ascii output matrix.
subroutine fill_ascii(ascii_values_mat, coords1, coords2, coords3, h, b, min, max, avg)
    real (kind = GRID_SR), dimension(2)      :: coords1, coords2, coords3    ! world-coordinates of the cells corners, with 1~(0,0)
    real (kind = GRID_SR), intent(in)        :: h, b, min, max, avg          ! water height, bathymetry, min height, max height, avg height
    type(data_for_ascii), intent(inout)           :: ascii_values_mat

    integer, dimension(2)               :: ascii_coords1, ascii_coords2, ascii_coords3           ! coords of the corresponding ascii matrix cell
    integer                             :: ascii_value                   ! associated priority of what is happening
    integer                             :: x_min, x_max, y_min, y_max, x, y
    integer, dimension(2)               :: vs1, vs2, vs3
    integer                             :: s, t, vz

    ascii_values_mat%h_min = min
    ascii_values_mat%h_max = max
    ascii_values_mat%h_avg = avg

    !compute ascii map coordinates from original grid coordinates, after
    ! computing the grid coordinates to 0-1-coordinates
    if (ascii_values_mat%s(1) == 0 .or. ascii_values_mat%s(2) == 0) then
        write(*,'(A)') "Error - scaling is zero"
        stop
    endif

    ascii_coords1(1) = ceiling(((coords1(1)-ascii_values_mat%a(1))/ascii_values_mat%s(1))*(size(ascii_values_mat%mat,dim=1)-1))+1
    ascii_coords1(2) = ceiling(((coords1(2)-ascii_values_mat%a(2))/ascii_values_mat%s(2))*(size(ascii_values_mat%mat,dim=2)-1))+1
    ascii_coords2(1) = ceiling(((coords2(1)-ascii_values_mat%a(1))/ascii_values_mat%s(1))*(size(ascii_values_mat%mat,dim=1)-1))+1
    ascii_coords2(2) = ceiling(((coords2(2)-ascii_values_mat%a(2))/ascii_values_mat%s(2))*(size(ascii_values_mat%mat,dim=2)-1))+1
    ascii_coords3(1) = ceiling(((coords3(1)-ascii_values_mat%a(1))/ascii_values_mat%s(1))*(size(ascii_values_mat%mat,dim=1)-1))+1
    ascii_coords3(2) = ceiling(((coords3(2)-ascii_values_mat%a(2))/ascii_values_mat%s(2))*(size(ascii_values_mat%mat,dim=2)-1))+1

    x_min = minval([ascii_coords1(1), ascii_coords2(1), ascii_coords3(1)])
    x_max = maxval([ascii_coords1(1), ascii_coords2(1), ascii_coords3(1)])
    y_min = minval([ascii_coords1(2), ascii_coords2(2), ascii_coords3(2)])
    y_max = maxval([ascii_coords1(2), ascii_coords2(2), ascii_coords3(2)])

    vs1 = [ascii_coords2(1) - ascii_coords1(1), ascii_coords2(2) - ascii_coords1(2)]
    vs2 = [ascii_coords3(1) - ascii_coords1(1), ascii_coords3(2) - ascii_coords1(2)]

    do x = x_min, x_max
        do y = y_min, y_max
            vs3 = [x-ascii_coords1(1), y-ascii_coords1(2)]
            vz = sign(1, (vs1(1)*vs2(2)-vs2(1)*vs1(2)))
            s = (vs3(1)*vs2(2)-vs2(1)*vs3(2))*vz
            t = (vs1(1)*vs3(2)-vs3(1)*vs1(2))*vz
            if (s>=0 .and. t>=0 .and. s+t<=abs(vs1(1)*vs2(2)-vs2(1)*vs1(2))) then
                !collect height & bathymetry
                ascii_values_mat%mat(x,y)%b = ascii_values_mat%mat(x,y)%b + b
                ascii_values_mat%mat(x,y)%h_sum = ascii_values_mat%mat(x,y)%h_sum + h
                ascii_values_mat%mat(x,y)%h_summands = ascii_values_mat%mat(x,y)%h_summands + 1
            end if
        end do
    end do

end subroutine

! writes the ascii output based on the mat
subroutine print_it(ascii_values_mat)
    type(data_for_ascii)                          :: ascii_values_mat

    integer                             :: i,j

    !printing the ascii output

    do i=1, size(ascii_values_mat%mat,dim=2)                                          ! "dim =" determines which matrix coordinate is meant
        do j=1, (size(ascii_values_mat%mat,dim=1) - 1)
            write(*,'(A,$)') which_ascii(ascii_values_mat, j, i)                       ! "$" forces that no line break is done
        end do
        write(*,'(A)') which_ascii(ascii_values_mat, (size(ascii_values_mat%mat,dim=1)), i)
    end do
    write(*,*)

    ! legend
    write(*,'(A)') ". = flat sea"
    write(*,'(A)') "_ = baseflow"
    write(*,'(A)') "~ = minor flood"
    write(*,'(A)') "^ = flood"
    write(*,'(A)') "| = coast"
    write(*,'(A)') "X = flooded coast"
    write(*,'(A)') "O = land"
    write(*,'(A)') "w = flooded land"
    write(*,'(A)') "# = cell error"

    ! min, avg, max
    write(*,'(A,$)') "Minimum height: "
    write(*,*) ascii_values_mat%h_min
    write(*,'(A,$)') "Average height: "
    write(*,*) ascii_values_mat%h_avg
    write(*,'(A,$)') "Maximum height: "
    write(*,*) ascii_values_mat%h_max


end subroutine

! determines the symbol based on bathymetry & height ---------- for internal use only
function which_ascii(ascii_values_mat, j, i)  result(symb)
    type(data_for_ascii), intent(in)              :: ascii_values_mat
    integer, intent(in)                 :: j, i

    type(bath_height_data)      :: basevalues
    character        :: symb
    real (kind = GRID_SR) :: loc_h_avg, loc_b_avg

    basevalues = ascii_values_mat%mat(j,i)

    if (basevalues%h_summands > 0) then
        !compute local average water height
        loc_h_avg = basevalues%h_sum / real(basevalues%h_summands, GRID_SR)
        loc_b_avg = basevalues%b / real(basevalues%h_summands, GRID_SR)
        if (loc_b_avg<0) then
            if (loc_h_avg>ascii_values_mat%eps) then
                if (loc_h_avg>((ascii_values_mat%h_max) / 2.0_GRID_SR)) then
                    symb = '^'                ! Major wave
                else
                    symb = '~'                ! Minor wave
                end if
            else
                if (loc_h_avg>(-ascii_values_mat%eps)) then
                    symb = '.'                ! flat sea
                else
                    symb = '_'                ! dry sea floor / negative wave
                end if
            end if
        else
            if (loc_b_avg < 0.25_GRID_SR .and. loc_b_avg > -0.25_GRID_SR) then
                if (loc_h_avg>ascii_values_mat%eps) then
                    symb = 'X'                ! Wave hits land = highest priority
                else
                    symb = '|'                ! Coast
                end if
            else
                if (loc_h_avg>ascii_values_mat%eps) then
                    symb = 'w'                ! Land is flooded
                else
                    symb = 'O'                ! Ususal land
                end if
            end if
        end if
    else
        symb = '#'
    end if
end function

subroutine reset(ascii_values_mat)
    type(data_for_ascii)                          :: ascii_values_mat
    integer                             :: i,j

    do i=1, size(ascii_values_mat%mat,dim=2)
        do j=1, size(ascii_values_mat%mat,dim=1)
            ascii_values_mat%mat(j,i)%b = -5.0_GRID_SR
            ascii_values_mat%mat(j,i)%h_sum = 0.0_GRID_SR
            ascii_values_mat%mat(j,i)%h_summands = 0
        end do
    end do

end subroutine

END MODULE
