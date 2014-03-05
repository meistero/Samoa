  ! statt h und b: min und max, wellengroesse trennen bei avg 
MODULE ascii_output

    implicit none
    
    type ascy                                                   
        integer, allocatable            :: mat (:,:)            ! matrix to save the ascii value
        double precision                :: eps, s(2), a(2), h_avg, h_min, h_max      ! epsilon to declare unmoving water, values for coordinate transformation, average, minimal and maximal water height
    end type                                                    !  with s = scaling factor, a = offset

    public
    
    contains
    
function create_ascy(dim_x, dim_y, eps, s, a) result (ascii)
    integer             :: dim_x, dim_y
    double precision    :: eps, s(2), a(2)

    type(ascy)          :: ascii
    integer             :: ascy_err
    
    if (allocated(ascii%mat)) then
        deallocate(ascii%mat)
    end if
    
    allocate (ascii%mat(dim_x, dim_y), stat = ascy_err)
    write(*,'(A)') "Ascii Matrix Allocated"      !--------------------------------------------------------------
    if (ascy_err > 0) then
        write(*,'(A)') "Error when trying to allocate for ascii output"
    end if
    
    ascii%mat = 1                           !setting the base value for the whole matrix
    ascii%eps = eps
    ascii%s = s
    ascii%a = a
    ascii%h_min = huge(1.0)
    ascii%h_max = -huge(1.0)
    ascii%h_avg = 0
    
end function

!when traversing the grid, call the following function from inside a cell to 
! transfer its data to the simplified ascii output matrix.
subroutine fill_sao(sao_values_mat, coords, h, b, min, max, avg)     
    double precision, dimension(2)      :: coords               ! world-coordinates of the cell, value has to be in (0,1)
    double precision, intent(in)        :: h, b, min, max, avg          ! water height, bathymetry, min height, max height, avg height
    type(ascy), intent(inout)           :: sao_values_mat
    
    integer, dimension(2)               :: sao_coords           ! coords of the corresponding ascii matrix cell
    integer                             :: sao_value            ! associated priority of what is happening   
    
    sao_values_mat%h_min = min
    sao_values_mat%h_max = max
    sao_values_mat%h_avg = avg
    
    
    !compute ascii map coordinates from original grid coordinates, after
    ! computing the grid coordinates to 0-1-coordinates
    if (sao_values_mat%s(1) == 0 .or. sao_values_mat%s(2) == 0) then 
        write(*,'(A)') "Error - scaling is zero"
        stop
    endif
    
    sao_coords(1) = ceiling(((coords(1)-sao_values_mat%a(1))/sao_values_mat%s(1))*size(sao_values_mat%mat,dim=1))           
    sao_coords(2) = ceiling(((coords(2)-sao_values_mat%a(2))/sao_values_mat%s(2))*size(sao_values_mat%mat,dim=2))
    
    
    !compute priority value for ascii_output
    if (b<0) then
        if (h>sao_values_mat%eps) then
            if (h>((sao_values_mat%h_max) / dble(2.0))) then
                sao_value = 7                ! Major wave
            else
                sao_value = 6                ! Minor wave
            end if
        else 
            if (h>(-sao_values_mat%eps)) then
                sao_value = 1                ! Ususal sea = lowest priority
            else
                sao_value = 2                ! Either sea floor falls dry or there is just less water
            end if
        end if        
    else
        if (b==0) then   
            if (h>sao_values_mat%eps) then
                sao_value = 8                ! Wave hits land = highest priority
            else
                sao_value = 4                ! Coast
            end if
        else
            if (h>sao_values_mat%eps) then
                sao_value = 5                ! Land is flooded
            else
                sao_value = 3                ! Ususal land
            end if
        end if
    end if
    
    !update the value if required 
     
    if (sao_value .gt. sao_values_mat%mat(sao_coords(1), sao_coords(2))) then
        sao_values_mat%mat(sao_coords(1), sao_coords(2)) = sao_value
    end if
    
end subroutine

! writes the ascii output based on the mat
subroutine print_it(sao_values_mat)                                                  
    type(ascy)                          :: sao_values_mat
    
    integer                             :: i,j
    
    !Max height, min height, avg height
    
    
    !printing the ascii output
    
    do i=1, size(sao_values_mat%mat,dim=2)                                          ! "dim =" determines which matrix coordinate is meant 
        do j=1, (size(sao_values_mat%mat,dim=1) - 1)
            write(*,'(A,$)') which_ascii(sao_values_mat%mat(j,i))                       ! "$" forces that no line break is done
        end do
        write(*,'(A)') which_ascii(sao_values_mat%mat(size(sao_values_mat%mat,dim=1), i))                                                    
    end do
    write(*,*)
    
    !legend
    write(*,'(A,$)') which_ascii(1)
    write(*,'(A)') " = usual sea"
    write(*,'(A,$)') which_ascii(2)
    write(*,'(A)') " = dry sea floor / negative wave"
    write(*,'(A,$)') which_ascii(3)
    write(*,'(A)') " = usual land"
    write(*,'(A,$)') which_ascii(4)
    write(*,'(A)') " = coast"
    write(*,'(A,$)') which_ascii(5)
    write(*,'(A)') " = flooded land"
    write(*,'(A,$)') which_ascii(6)
    write(*,'(A)') " = minor wave"
    write(*,'(A,$)') which_ascii(7)
    write(*,'(A)') " = major wave"
    write(*,'(A,$)') which_ascii(8)
    write(*,'(A)') " = flooded coast"
    
    !min, max, avg
    write(*,'(A,$)') "Maximum height: "
    write(*,*) sao_values_mat%h_max
    write(*,'(A,$)') "Minimum height: "
    write(*,*) sao_values_mat%h_min
    write(*,'(A,$)') "Average height: "
    write(*,*) sao_values_mat%h_avg
    
    sao_values_mat%mat = 1                                                  ! clearing map after printing
                                                                                            
end subroutine

! determines the symbol based on bathymetry & height ---------- for internal use only
function which_ascii(sao_value)  result(symb)                                                   
    integer, intent(in)             :: sao_value
                                                                    
    character        :: symb
    
    select case (sao_value)
        case (1)
            symb = '.'              ! Usual sea ( = lowest priority)
        case(2)
            symb = '_'              ! Either sea floor falls dry or there is just less water
        case (3)
            symb = 'O'              ! Usual land
        case (4)
            symb = '|'              ! Coast
        case (5)
            symb = 'w'              ! Flooded land
        case (6)
            symb = '~'              ! Wave
        case (7)
            symb = '^'              ! "Bigger wave"
        case (8)
            symb = 'X'              ! Wave hits coast ( = highest priority)
    end select
                
end function
        
	END MODULE