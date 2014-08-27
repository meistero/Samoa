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

!for the SWE scenario, a skeleton operator is required
#if defined(_SWE)
#	define _USE_SKELETON_OP
#endif

#if defined(_NUMA)
#	define _USE_SKELETON_OP
#endif

#if defined(_FLASH)
#	define _USE_SKELETON_OP
#endif

#define _SWE_ORDER 					0
#define _SWE_CELL_SIZE				((_SWE_ORDER + 1) * (_SWE_ORDER + 2)) / 2
#define _SWE_EDGE_SIZE 				(_SWE_ORDER + 1)

#define _NUMA_ORDER 					0
#define _NUMA_CELL_SIZE				((_NUMA_ORDER + 1) * (_NUMA_ORDER + 2)) / 2
#define _NUMA_EDGE_SIZE 			(_NUMA_ORDER + 1)

#define _HEAT_EQ_ORDER 				2
#define _HEAT_EQ_SIZE				((_HEAT_EQ_ORDER + 1) * (_HEAT_EQ_ORDER + 2)) / 2
#define _HEAT_EQ_CELL_SIZE			((_HEAT_EQ_ORDER - 1) * (_HEAT_EQ_ORDER - 2)) / 2
#define _HEAT_EQ_EDGE_SIZE 			(_HEAT_EQ_ORDER - 1)
#define _HEAT_EQ_NODE_SIZE			1

#define _FLASH_ORDER          		        1
#define _FLASH_CELL_SIZE                        ((_FLASH_ORDER + 1) * (_FLASH_ORDER + 2)) / 2
#define _FLASH_EDGE_SIZE                        _FLASH_CELL_SIZE
#define _FLASH_CELL_QUAD_SIZE      		_FLASH_CELL_SIZE
#define _FLASH_EDGE_QUAD_SIZE  			_FLASH_CELL_SIZE

#define _DARCY_P_ORDER				1
#define _DARCY_P_SIZE				((_DARCY_P_ORDER + 1) * (_DARCY_P_ORDER + 2)) / 2
#define _DARCY_P_CELL_SIZE			((_DARCY_P_ORDER - 1) * (_DARCY_P_ORDER - 2)) / 2
#define _DARCY_P_EDGE_SIZE			(_DARCY_P_ORDER - 1)
#define _DARCY_P_NODE_SIZE			1

#define _DARCY_U_ORDER				_DARCY_P_ORDER - 1
#if (_DARCY_U_ORDER > 0)
#	define _DARCY_U_SIZE			((_DARCY_P_ORDER + 1) * (_DARCY_P_ORDER + 2)) / 2
#	define _DARCY_U_CELL_SIZE		((_DARCY_P_ORDER - 1) * (_DARCY_P_ORDER - 2)) / 2
#	define _DARCY_U_EDGE_SIZE		(_DARCY_P_ORDER - 1)
#	define _DARCY_U_NODE_SIZE		1
#else
#	define _DARCY_U_SIZE			1
#	define _DARCY_U_CELL_SIZE		1
#	define _DARCY_U_EDGE_SIZE		0
#	define _DARCY_U_NODE_SIZE		0
#endif

#define _DARCY_FLOW_ORDER			1
#define _DARCY_FLOW_SIZE			3
#define _DARCY_FLOW_CELL_SIZE		0
#define _DARCY_FLOW_EDGE_SIZE		0
#define _DARCY_FLOW_NODE_SIZE		1

!> Assertion macros:

!> Checks for a condition to be true and raises an artificial divide-by-zero exception for error handling
#if defined(_ASSERT)
#	define assert(x)				if (.not. (x)) then; PRINT '(a, a, i0, a, a)', __FILE__, "(", __LINE__, "): Assertion failure: ", #x; flush(6); PRINT *, 0 / 0; end if
#	define assert_pure(x)			if (.not. (x)) then; call raise_error(); end if

#	define assert_eq(x, y)			if (.not. (x .eq. y)) then; PRINT '(a, a, i0, a, a, a, a, a, g0, a, g0)', __FILE__, "(", __LINE__, "): Assertion failure: ", #x, " == ", #y, ": ", x, " == ", y; flush(6); PRINT *, 0 / 0; end if
#	define assert_lt(x, y)			if (.not. (x .lt. y)) then; PRINT '(a, a, i0, a, a, a, a, a, g0, a, g0)', __FILE__, "(", __LINE__, "): Assertion failure: ", #x, " < ", #y, ": ", x, " < ", y; flush(6); PRINT *, 0 / 0; end if
#	define assert_le(x, y)			if (.not. (x .le. y)) then; PRINT '(a, a, i0, a, a, a, a, a, g0, a, g0)', __FILE__, "(", __LINE__, "): Assertion failure: ", #x, " <= ", #y, ": ", x, " <= ", y; flush(6); PRINT *, 0 / 0; end if
#	define assert_gt(x, y)			if (.not. (x .gt. y)) then; PRINT '(a, a, i0, a, a, a, a, a, g0, a, g0)', __FILE__, "(", __LINE__, "): Assertion failure: ", #x, " > ", #y, ": ", x, " > ", y; flush(6); PRINT *, 0 / 0; end if
#	define assert_ge(x, y)			if (.not. (x .ge. y)) then; PRINT '(a, a, i0, a, a, a, a, a, g0, a, g0)', __FILE__, "(", __LINE__, "): Assertion failure: ", #x, " >= ", #y, ": ", x, " >= ", y; flush(6); PRINT *, 0 / 0; end if
#	define assert_ne(x, y)			if (.not. (x .ne. y)) then; PRINT '(a, a, i0, a, a, a, a, a, g0, a, g0)', __FILE__, "(", __LINE__, "): Assertion failure: ", #x, " /= ", #y, ": ", x, " /= ", y; flush(6); PRINT *, 0 / 0; end if

#	define assert_v(x)				if (.not. all(x)) then; PRINT '(a, a, i0, a, a)', __FILE__, "(", __LINE__, "): Assertion failure: ", #x; flush(6); PRINT *, 0 / 0; end if
#	define assert_veq(x, y)			if (.not. all(x .eq. y)) then; PRINT '(a, a, i0, a, a, a, a, a, (X, g0), (X, g0))', __FILE__, "(", __LINE__, "): Assertion failure: ", #x, " == ", #y, ": ", x, y; flush(6); PRINT *, 0 / 0; end if
#	define assert_vne(x, y)			if (.not. any(x .ne. y)) then; PRINT '(a, a, i0, a, a, a, a, a, (X, g0), (X, g0))', __FILE__, "(", __LINE__, "): Assertion failure: ", #x, " != ", #y, ": ", x, y; flush(6); PRINT *, 0 / 0; end if
#	define mpi_isend				mpi_issend
#else
#	define assert(x)
#	define assert_pure(x)

#	define assert_eq(x, y)
#	define assert_lt(x, y)
#	define assert_le(x, y)
#	define assert_gt(x, y)
#	define assert_ge(x, y)
#	define assert_ne(x, y)

#	define assert_v(x)
#	define assert_veq(x, y)
#	define assert_vne(x, y)
#endif

#	define try(x, str_exception)	if (.not. (x)) then; PRINT '(a, "(", i0, "): Exception: ", a, ": ", a)', __FILE__, __LINE__, str_exception, #x; flush(6); PRINT *, 0 / 0; end if

!> Log file macros
#define _log_open_file				call log_open_file
#define _log_close_file				call log_close_file
#define _log_write(dl, f)			if (_DEBUG_LEVEL .ge. dl) write(g_log_file_unit,'(A, A, I0, A, I0, A)',advance='no') term_color(omp_get_thread_num() * size_MPI + rank_MPI), "(r", rank_MPI, ",t", omp_get_thread_num(), ") "; if (_DEBUG_LEVEL .ge. dl) write(g_log_file_unit, f)
