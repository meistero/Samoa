! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


! Global compiler control statements and variables

#define COMPILATION_CONTROL

! define debug level if not already defined
#if !defined(_DEBUG_LEVEL)
#	define _DEBUG_LEVEL				1
#endif

! if ASAGI with NUMA support is requested, activate ASAGI
#if defined(_ASAGI_NUMA)
#	define _ASAGI
#endif

!> if true, element coordinates are not computed but stored in the nodes instead (which allows for irregular grids)
!#define _STORE_NODE_COORDS

# define _DARCY_INJ_PRESSURE
# define _DARCY_PROD_PRESSURE

!for the SWE scenario, a skeleton operator is required
#if defined(_SWE) || defined(_NUMA) || defined(_FLASH)
#	define _USE_SKELETON_OP
#endif

#define _SWE_ORDER 					0
#define _SWE_CELL_SIZE				((_SWE_ORDER + 1) * (_SWE_ORDER + 2)) / 2
#define _SWE_EDGE_SIZE 				(_SWE_ORDER + 1)

#define _NUMA_ORDER 				0
#define _NUMA_CELL_SIZE				((_NUMA_ORDER + 1) * (_NUMA_ORDER + 2)) / 2
#define _NUMA_EDGE_SIZE 			(_NUMA_ORDER + 1)

#define _HEAT_EQ_ORDER 				2
#define _HEAT_EQ_SIZE				((_HEAT_EQ_ORDER + 1) * (_HEAT_EQ_ORDER + 2)) / 2
#define _HEAT_EQ_CELL_SIZE			((_HEAT_EQ_ORDER - 1) * (_HEAT_EQ_ORDER - 2)) / 2
#define _HEAT_EQ_EDGE_SIZE 			(_HEAT_EQ_ORDER - 1)
#define _HEAT_EQ_NODE_SIZE			1

#define _FLASH_ORDER 				0
#define _FLASH_CELL_SIZE			((_FLASH_ORDER + 1) * (_FLASH_ORDER + 2)) / 2
#define _FLASH_EDGE_SIZE 			_FLASH_CELL_SIZE
#define _FLASH_EDGE_QUAD_SIZE			_FLASH_CELL_SIZE

#if defined(_DARCY)
#   if defined(_ASAGI)
#	    define _DARCY_INJECTOR_WELLS    1
#	    define _DARCY_PRODUCER_WELLS    4
#   else
#	    define _DARCY_INJECTOR_WELLS    1
#	    define _DARCY_PRODUCER_WELLS    1
#   endif
#endif

!compiler-dependent macros (traditional vs. modern preprocessor)
#	define _id(x) x

#if defined(__GFORTRAN__)
#	define _raise()	                    call abort
#	define _raise_pure()                call raise_error()
#	define _stringify(x)	            "x"
#	define _conc(x, y)	                _id(x)_id(y)
#	define _conc3(x, y, z)	            _id(x)_id(y)_id(z)
#	define _conc4(x, y, z, a)           _id(x)_id(y)_id(z)_id(a)
#	define _conc5(x, y, z, a, b)        _id(x)_id(y)_id(z)_id(a)_id(b)
#	define _conc6(x, y, z, a, b, c)     _id(x)_id(y)_id(z)_id(a)_id(b)_id(c)
#	define _conc7(x, y, z, a, b, c, d)  _id(x)_id(y)_id(z)_id(a)_id(b)_id(c)_id(d)
#else
#	define _raise() 	                print *, 0 / 0
#	define _raise_pure() 	            call raise_error()
#	define _stringify(x)	            #x
#	define _conc(x, y)	                x##y
#	define _conc3(x, y, z)	            x##y##z
#	define _conc4(x, y, z, a)	        x##y##z##a
#	define _conc5(x, y, z, a, b)	    x##y##z##a##b
#	define _conc6(x, y, z, a, b, c)	    x##y##z##a##b##c
#	define _conc7(x, y, z, a, b, c, d)	x##y##z##a##b##c##d
#endif

!> Assertion macros:

!> Checks for a condition to be true and raises an artificial divide-by-zero exception for error handling
#if defined(_ASSERT)
#	define assert(x)				if (.not. (x)) then; PRINT '(a, a, i0, a, a)', __FILE__, "(", __LINE__, "): Assertion failure: ", _stringify(x); flush(6); _raise(); end if
#	define assert_pure(x)			if (.not. (x)) then; _raise_pure(); end if

#	define assert_eq(x, y)			if (.not. (x .eq. y)) then; PRINT '(a, a, i0, a, a, a, a, a, g0, a, g0)', __FILE__, "(", __LINE__, "): Assertion failure: ", _stringify(x), " == ", _stringify(y), ": ", x, " == ", y; flush(6); _raise(); end if
#	define assert_eqv(x, y)			if (.not. (x .eqv. y)) then; PRINT '(a, a, i0, a, a, a, a, a, g0, a, g0)', __FILE__, "(", __LINE__, "): Assertion failure: ", _stringify(x), " == ", _stringify(y), ": ", x, " == ", y; flush(6); _raise(); end if
#	define assert_lt(x, y)			if (.not. (x .lt. y)) then; PRINT '(a, a, i0, a, a, a, a, a, g0, a, g0)', __FILE__, "(", __LINE__, "): Assertion failure: ", _stringify(x), " < ", _stringify(y), ": ", x, " < ", y; flush(6); _raise(); end if
#	define assert_le(x, y)			if (.not. (x .le. y)) then; PRINT '(a, a, i0, a, a, a, a, a, g0, a, g0)', __FILE__, "(", __LINE__, "): Assertion failure: ", _stringify(x), " <= ", _stringify(y), ": ", x, " <= ", y; flush(6); _raise(); end if
#	define assert_gt(x, y)			if (.not. (x .gt. y)) then; PRINT '(a, a, i0, a, a, a, a, a, g0, a, g0)', __FILE__, "(", __LINE__, "): Assertion failure: ", _stringify(x), " > ", _stringify(y), ": ", x, " > ", y; flush(6); _raise(); end if
#	define assert_ge(x, y)			if (.not. (x .ge. y)) then; PRINT '(a, a, i0, a, a, a, a, a, g0, a, g0)', __FILE__, "(", __LINE__, "): Assertion failure: ", _stringify(x), " >= ", _stringify(y), ": ", x, " >= ", y; flush(6); _raise(); end if
#	define assert_ne(x, y)			if (.not. (x .ne. y)) then; PRINT '(a, a, i0, a, a, a, a, a, g0, a, g0)', __FILE__, "(", __LINE__, "): Assertion failure: ", _stringify(x), " /= ", _stringify(y), ": ", x, " /= ", y; flush(6); _raise(); end if

#	define assert_v(x)				if (.not. all(x)) then; PRINT '(a, a, i0, a, a)', __FILE__, "(", __LINE__, "): Assertion failure: ", _stringify(x); flush(6); _raise(); end if
#	define assert_veq(x, y)			if (.not. all(x .eq. y)) then; PRINT '(a, a, i0, a, a, a, a, a, (X, g0), (X, g0))', __FILE__, "(", __LINE__, "): Assertion failure: ", _stringify(x), " == ", _stringify(y), ": ", x, y; flush(6); _raise(); end if
#	define assert_vne(x, y)			if (.not. any(x .ne. y)) then; PRINT '(a, a, i0, a, a, a, a, a, (X, g0), (X, g0))', __FILE__, "(", __LINE__, "): Assertion failure: ", _stringify(x), " != ", _stringify(y), ": ", x, y; flush(6); _raise(); end if
#else
#	define assert(x)
#	define assert_pure(x)

#	define assert_eq(x, y)
#	define assert_eqv(x, y)
#	define assert_lt(x, y)
#	define assert_le(x, y)
#	define assert_gt(x, y)
#	define assert_ge(x, y)
#	define assert_ne(x, y)

#	define assert_v(x)
#	define assert_veq(x, y)
#	define assert_vne(x, y)
#endif

#define try(x, str_exception)	if (.not. (x)) then; PRINT '(a, "(", i0, "): Exception: ", a, ": ", a)', __FILE__, __LINE__, str_exception, _stringify(x) ; flush(6); _raise(); end if

!> Log file macros
#define _log_open_file				call log_open_file
#define _log_close_file				call log_close_file
#define _log_write(dl, f)			if (_DEBUG_LEVEL .ge. dl) write(g_log_file_unit,'(A, A, I0, A, I0, A)',advance='no') term_color(omp_get_thread_num() * size_MPI + rank_MPI), "(r", rank_MPI, ",t", omp_get_thread_num(), ") "; if (_DEBUG_LEVEL .ge. dl) write(g_log_file_unit, f)

!Standard units:

!> Unit meter, defined as the width and height of the full grid
#define _UM     1.0_GRID_SR

!> Second in simulation time
#define _S      1.0_GRID_SR

!> Kilogram (standard mass)
#define _KG     1.0_GRID_SR

!Derived units and their conversion rules:

!> Meter
#define _M      (_UM / cfg%scaling)
!> Minutes (exact)
#define _MIN    (_S * 60.0_GRID_SR)
!> Hours (exact)
#define _H      (_MIN * 60.0_GRID_SR)
!> Days (exact)
#define _D      (_H * 24.0_GRID_SR)
!> Newton (exact)
#define _N      (_KG * _M / (_S * _S))
!> Pascal (exact)
#define _PA     (_N / (_M * _M))
!> Centipoise (exact)
#define _CP     (_PA * 1.0e-3_GRID_SR)
!> Barrel Oil (exact)
#define _BBL    (((_INCH) ** 3) * 9702.0_GRID_SR)
!> Inch (exact)
#define _INCH   (_M * 0.0254_GRID_SR)
!> Feet (exact)
#define _FT     (_INCH * 12.0_GRID_SR)
!> Pound (exact)
#define _LB     (_KG * 0.45359237_GRID_SR)
!> Pound Force (exact)
#define _LBF    (_LB * _G)
!> Pound Force Per Square Inch (exact)
#define _PPSI   (_LBF / (_INCH ** 2))
!> Millidarcy (exact)
#define _MDY    (_DY * 1.0e-3_GRID_SR)
!> Darcy (exact)
#define _DY     ((_M ** 2) / 1.01325e12_GRID_SR)

!Physical constants:

!> Standard gravitational field (exact)
#define _G      9.80665_GRID_SR * _M / (_S ** 2)
