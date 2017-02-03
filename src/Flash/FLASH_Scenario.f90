! Artificial scenario selector for the FLASH scenario. 
! To add a new scenario:
! 1) create a new module according to the template below
! 2) add an option for it in the Sconstruct file (FLASH_scenario argument) and define a macro if it is chosen
! 3) "USE" its module in the FLASH_Scenario module (at the end of this file)
! 4) Don't forget to use #if defined(_NEW_MACRO) to avoid conflicts!

#include "Compilation_control.f90"


!*******************
!* MODULE TEMPLATE *
!*******************
#if 0
MODULE FLASH_Scenario_template
    use Samoa_FLASH
    public FLASH_Scenario_get_scaling, FLASH_Scenario_get_offset, FLASH_Scenario_get_bathymetry, FLASH_Scenario_get_initial_Q
    contains

    function FLASH_Scenario_get_scaling() result(scaling)
        real (kind = GRID_SR) :: scaling
        
        scaling = 10.0_GRID_SR
    end function

    function FLASH_Scenario_get_offset() result(offset)
        real (kind = GRID_SR) :: offset(2)
        
        offset = FLASH_Scenario_get_scaling() * [-0.5_GRID_SR, -0.5_GRID_SR]
    end function
    
    function FLASH_Scenario_get_bathymetry(x) result(bathymetry)
        real (kind = GRID_SR), intent(in) :: x(2)
        real (kind = GRID_SR) :: bathymetry
        
        bathymetry = 0.0_GRID_SR
    end function
    
    function FLASH_Scenario_get_initial_Q(x) result(Q)
        real (kind = GRID_SR), intent(in) :: x(2)
        type(t_dof_state) :: Q
        
        Q%p = [0.0_GRID_SR, 0.0_GRID_SR]
        Q%h = 0.0_GRID_SR
    end function

END MODULE FLASH_Scenario_radial_dam_break
#endif



!********************
!* Radial Dam Break *
!********************
#if defined (_FLASH_SCENARIO_RADIAL_DAM_BREAK)
MODULE FLASH_Scenario_radial_dam_break
    use Samoa_FLASH
    public FLASH_Scenario_get_scaling, FLASH_Scenario_get_offset, FLASH_Scenario_get_bathymetry, FLASH_Scenario_get_initial_Q
    contains

    function FLASH_Scenario_get_scaling() result(scaling)
        real (kind = GRID_SR) :: scaling
        
        scaling = 10.0_GRID_SR
    end function

    function FLASH_Scenario_get_offset() result(offset)
        real (kind = GRID_SR) :: offset(2)
        
        offset = FLASH_Scenario_get_scaling() * [-0.5_GRID_SR, -0.5_GRID_SR]
    end function
    
    function FLASH_Scenario_get_bathymetry(x) result(bathymetry)
        real (kind = GRID_SR), intent(in) :: x(2)
        real (kind = GRID_SR) :: bathymetry
        
        if (x(1)*x(1) + x(2)*x(2) < 9.0_GRID_SR) then
            bathymetry = -5.0_GRID_SR
        else 
            bathymetry = -10.0_GRID_SR
        end if
    end function
    
    function FLASH_Scenario_get_initial_Q(x) result(Q)
        real (kind = GRID_SR), intent(in) :: x(2)
        type(t_dof_state) :: Q
        
        Q%p = [0.0_GRID_SR, 0.0_GRID_SR]
        
        if (x(1)*x(1) + x(2)*x(2) < 0.25_GRID_SR) then
            Q%h = 10.0_GRID_SR
        else 
            Q%h = 0.0_GRID_SR
        end if
    end function

END MODULE FLASH_Scenario_radial_dam_break
#endif


!********************
!* Linear Dam Break *
!********************
#if defined (_FLASH_SCENARIO_LINEAR_DAM_BREAK)
MODULE FLASH_Scenario_linear_dam_break
    use Samoa_FLASH
    public FLASH_Scenario_get_scaling, FLASH_Scenario_get_offset, FLASH_Scenario_get_bathymetry, FLASH_Scenario_get_initial_Q
    contains

    function FLASH_Scenario_get_scaling() result(scaling)
        real (kind = GRID_SR) :: scaling
        
        scaling = 100000.0_GRID_SR
    end function

    function FLASH_Scenario_get_offset() result(offset)
        real (kind = GRID_SR) :: offset(2)
        
        offset = FLASH_Scenario_get_scaling() * [-0.5_GRID_SR, -0.5_GRID_SR]
    end function
    
    function FLASH_Scenario_get_bathymetry(x) result(bathymetry)
        real (kind = GRID_SR), intent(in) :: x(2)
        real (kind = GRID_SR) :: bathymetry
        
        bathymetry = -10.0_GRID_SR
    end function
    
    function FLASH_Scenario_get_initial_Q(x) result(Q)
        real (kind = GRID_SR), intent(in) :: x(2)
        type(t_dof_state) :: Q
        
        real (kind = GRID_SR), parameter :: hL = 10.0_GRID_SR, hR = 0.0_GRID_SR
        
        Q%p = [0.0_GRID_SR, 0.0_GRID_SR]
        
        if (x(1) < 0.0_GRID_SR) then
            Q%h = hL
        else 
            Q%h = hR
        end if
    end function

END MODULE FLASH_Scenario_linear_dam_break
#endif

!********************
!* Oscillating Lake *
!********************
! Proposed in: 
! [1] Gallardo et al., 2007. On a well-balanced high-order finite volume scheme for shallow water equations with topography and dry areas. Journal of Computational Physics, volume 227
! [2] Meister & Ortleb, 2014. On unconditionally positive implicit time integration for the dg scheme applied to shallow water flows. International Journal for Numerical Methods in Fluids, volume 76, 69-94.
!
! Domain: [-2, 2]²
! Bathymetry: b(x,y) = 0.1 (x² + y²)
!
! Analytic solution:
! H(x,y,t) = max{ 0.0, 0.05 * (2x cos(wt) + 2y sin(wt) + 0.075 - b(x,y) ) }
! u(x,y,t) = -0.5w sin(wt)
! v(x,y,t) =  0.5w cos(wt)
! --> with w = sqrt( 0.2 g )

#if defined (_FLASH_SCENARIO_OSCILLATING_LAKE)
MODULE FLASH_Scenario_oscillating_lake
    use Samoa_FLASH
    public FLASH_Scenario_get_scaling, FLASH_Scenario_get_offset, FLASH_Scenario_get_bathymetry, FLASH_Scenario_get_initial_Q
    contains

    function FLASH_Scenario_get_scaling() result(scaling)
        real (kind = GRID_SR) :: scaling
        
        scaling = 4.0_GRID_SR
    end function

    function FLASH_Scenario_get_offset() result(offset)
        real (kind = GRID_SR) :: offset(2)
        
        offset = FLASH_Scenario_get_scaling() * [-0.5_GRID_SR, -0.5_GRID_SR]
    end function
    
    function FLASH_Scenario_get_bathymetry(x) result(bathymetry)
        real (kind = GRID_SR), intent(in) :: x(2)
        real (kind = GRID_SR) :: bathymetry
        
        bathymetry = 0.1 * (x(1)*x(1) + x(2)*x(2))
    end function
    
    function FLASH_Scenario_get_initial_Q(x) result(Q)
        real (kind = GRID_SR), intent(in) :: x(2)
        type(t_dof_state) :: Q
        double precision :: w, t, sinwt, coswt, b
        
        b = FLASH_Scenario_get_bathymetry(x)
        
        t = 0.0
        w = sqrt(0.2 * g)
        sinwt = 0.0 ! t = 0
        coswt = 1.0 ! t = 0
        
        Q%h = max( 0.0_GRID_SR, 0.05 * (2*x(1)*coswt + 2*x(2)*sinwt) + 0.075  - b )
        
        Q%p(1) = -0.5*w*sinwt * (Q%h) 
        Q%p(2) =  0.5*w*coswt * (Q%h)

        
        Q%h = Q%h + b
        
    end function

END MODULE FLASH_Scenario_oscillating_lake
#endif


MODULE FLASH_Scenario

#   if defined(_FLASH_SCENARIO_RADIAL_DAM_BREAK)
        USE FLASH_Scenario_radial_dam_break
#   elif defined(_FLASH_SCENARIO_LINEAR_DAM_BREAK)
        USE FLASH_Scenario_linear_dam_break
#   elif defined(_FLASH_SCENARIO_OSCILLATING_LAKE)
        USE FLASH_Scenario_oscillating_lake
#   endif

END MODULE FLASH_Scenario
