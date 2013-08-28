! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

PROGRAM gridtest
	USE SFC_traversal
	USE M_kracken

	implicit none

	! command line arguments
	logical					                            :: l_help, l_version
	integer          									:: i, i_error, i_threads
	character(256)                                      :: arguments

	!define default command arguments and default values for all scenarios
    write(arguments, '(A, I0)') "-v .false. --version .false. -h .false. --help .false. -asagihints 2 -noprint .false. -sections 4 -threads ", omp_get_max_threads()

	!define additional command arguments and default values depending on the choice of the scenario
#	if defined(_DARCY)
        write(arguments, '(A, A)') trim(arguments), "  -dmin 1 -dmax 14 -tsteps -1 -tmax 2.0e1 -tout -1.0"
#	elif defined(_HEAT_EQ)
        write(arguments, '(A, A)') trim(arguments), "  -dmin 1 -dmax 16 -tsteps -1 -tmax 1.0 -tout -1.0"
#	elif defined(_SWE)
        write(arguments, '(A, A)') trim(arguments), "  -dmin 2 -dmax 14 -tsteps -1 -tmax 2.0e-3 -tout -1.0"
#	elif defined(_NUMA)
        write(arguments, '(A, A)') trim(arguments), "  -dmin 2 -dmax 14 -tsteps -1 -tmax 5 -tout -1.0"
#	elif defined(_TESTS)
        write(arguments, '(A, A)') trim(arguments), "  -dmin 2 -dmax 14 -tsteps -1 -tmax 5 -tout -1.0"
#	elif defined(_PYOP2)
        write(arguments, '(A, A)') trim(arguments), "  -dmin 2 -dmax 14 -tsteps -1 -tmax 5 -tout -1.0"
#	else
#       error No scenario selected!
#	endif

	call kracken('samoa', arguments)

	!  get values
	l_help = lget('samoa_-help') .or. lget('samoa_h')
	l_version = lget('samoa_-version') .or. lget('samoa_v')
	i_min_depth = iget('samoa_dmin')
	i_max_depth = iget('samoa_dmax')
	i_max_time_steps = iget('samoa_tsteps')
	r_max_time = rget('samoa_tmax')
	r_output_time_step = rget('samoa_tout')
	l_log = lget('samoa_noprint')
	i_threads = iget('samoa_threads')
	i_sections_per_thread = iget('samoa_sections')
	i_asagi_mode = iget('samoa_asagihints')

	!if the version option was set was called, display program version
	if (l_version) then
		PRINT '(A)', " sam(oa)²: Space filling curves and Adaptive Meshes for Oceanic and Other Applications"
		PRINT '(A, I0, A, I0, A, I0)', " version ", 0, ".", 5, ".", 0
	end if

	!if the help option was set, display the list of arguments
	if (l_help) then
		PRINT '(A)',            " Usage: samoa [--help | -h] | [--version | -v] | [-asagihints <value>] [-dstart <value>] [-dmin <value>] [-dmax <value>] [-tsteps <value>] [-tmax <value>] [-tout <value>] [-threads <value>] [-sections <value>] [-noprint]"
		PRINT '(A)',            ""
		PRINT '(A)',            " Arguments:"
		PRINT '(A, I0, A)',     " 	-asagihints <value>     ASAGI mode (0: default, 1: pass thorugh, 2: nompi, 3: no mpi + small cache, 4: large grid) (value: ", i_asagi_mode, ")"
		PRINT '(A, I0, A)',     " 	-dmin <value>           minimum grid depth (value: ", i_min_depth, ")"
		PRINT '(A, I0, A)',     "	-dmax <value>           maximum grid depth (value: ", i_max_depth, ")"
		PRINT '(A, I0, A)',     "	-tsteps <value>         maximum number of time steps, less than 0: not defined (value: ", i_max_time_steps, ")"
		PRINT '(A, ES8.1, A)',  "	-tmax <value>           maximum simulation time in seconds, less than 0: not defined (value: ", r_max_time, ")"
		PRINT '(A, ES8.1, A)',  "	-tout <value>           output time step in seconds, less than 0: not defined (value: ", r_output_time_step, ")"
		PRINT '(A, I0, A)',     "	-threads <value>        number of OpenMP threads (value: ", i_threads, ")"
		PRINT '(A, I0, A)',     "	-sections <value>       number of grid sections per OpenMP thread (value: ", i_sections_per_thread, ")"
		PRINT '(A)',            "	-noprint                print log to file instead of console"
		PRINT '(A)',            "	--help, -h              display this help and exit"
		PRINT '(A)',            "	--version, -v           output version information and exit"
    end if

    !stop if the version or help command were called
	if (l_help .or. l_version) then
        stop
    end if

    call omp_set_num_threads(i_threads)

    !init mpi
    call init_mpi()

    !init element transformation data
    call init_transform_data()

    call sfc_generic()

    call finalize_mpi()

	stop
end PROGRAM gridtest
