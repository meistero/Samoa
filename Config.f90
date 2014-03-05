
! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE

#include "Compilation_control.f90"

#if defined(_ASAGI)
#	include "asagi.f90"
#endif

module config
    use omp_lib
	use M_kracken
	use Tools_log

#	if defined(_ASAGI)
		use asagi
#	endif

    implicit none

    private
    public cfg, omp_tasks, i_thread, rank_MPI, size_MPI

    type t_config
        integer			                        :: i_threads			                            !< number of OpenMP threads
        integer			                        :: i_ranks			                                !< number of MPI ranks
        integer (kind = 1)			            :: i_sections_per_thread						    !< number of sections per thread
        integer (kind = selected_int_kind(8))   :: i_max_time_steps				                    !< number of simulation time steps
        double precision                        :: r_max_time, r_output_time_step					!< maximum simulation time, outpout time step
        logical			                        :: l_log                                            !< if true, a log file is used
        integer (kind = 1)                      :: i_min_depth, i_max_depth			                !< minimum and maximum scenario depth
        integer			        	            :: i_asagi_mode			                		    !< ASAGI mode

        double precision                        :: scaling, offset(2)                               !< grid scaling and offset

#    	if defined(_DARCY)
            character(256)                      :: s_permeability_file                              !< permeability file
 			integer					 		    :: afh_permeability			                        !< asagi file handle to permeability data
            integer			        	        :: i_lsolver			                		    !< linear solver

			double precision			        :: r_epsilon				                        !< linear solver error bound
			double precision				    :: r_rel_permeability		                        !< relative permeability of the entering fluid
			double precision				    :: r_rho					                        !< fluid density
			double precision				    :: r_p0			                                    !< initial boundary pressure difference
#    	elif defined(_SWE)
            character(256)                      :: s_bathymetry_file                                !< bathymetry file
            character(256)                      :: s_displacement_file                              !< displacement file
 			integer					 		    :: afh_displacement			                        !< asagi file handle to displacement data
 			integer					 		    :: afh_bathymetry			                        !< asagi file handle to bathymetry da
#    	elif defined(_FLASH)
            character(256)                      :: s_bathymetry_file                                !< bathymetry file
            character(256)                      :: s_displacement_file                              !< displacement file
 			integer					 		    :: afh_displacement			                        !< asagi file handle to displacement data
 			integer					 		    :: afh_bathymetry			                        !< asagi file handle to bathymetry data
#       endif

        contains

        procedure, pass :: read_from_arguments
        procedure, pass :: print => config_print
    end type

    type(t_config) :: cfg
    !$omp threadprivate(cfg)

#   if defined(_OPENMP_TASKS)
        logical, parameter                      :: omp_tasks = .true.
#   else
        logical, parameter                      :: omp_tasks = .false.
#   endif

	integer                     :: i_thread = 1     !start counting at 1
	!$omp threadprivate(i_thread)

	integer 					:: rank_MPI = 0
	integer 					:: size_MPI = 1

    contains

    subroutine read_from_arguments(config)
        class(t_config), intent(inout)          :: config

        logical					                :: l_help, l_version
        integer          					    :: i, i_error
        character(256)                          :: arguments
        character(64)                           :: lsolver_to_char(0:2) = ["Jacobi", "CG", "Pipelined CG"]
        character(64)                           :: asagi_mode_to_char(0:4) = ["default", "pass through", "no mpi", "no mpi + small cache", "large grid"]

        !define default command arguments and default values for all scenarios
        write(arguments, '(A, I0)') "-v .false. --version .false. -h .false. --help .false. -asagihints 2 -noprint .false. -sections 4 -threads ", omp_get_max_threads()

        !define additional command arguments and default values depending on the choice of the scenario
#    	if defined(_DARCY)
            write(arguments, '(A, A)') trim(arguments), "  -dmin 1 -dmax 14 -tsteps -1 -tmax 2.0e1 -tout -1.0 -fperm data/darcy_benchmark/perm.nc -p0 1.0e6 -epsilon 1.0e-5 -rho 0.2 -k_rel 1.5 -lsolver 2"
#    	elif defined(_HEAT_EQ)
            write(arguments, '(A, A)') trim(arguments), "  -dmin 1 -dmax 16 -tsteps -1 -tmax 1.0 -tout -1.0"
#    	elif defined(_SWE)
            write(arguments, '(A, A)') trim(arguments), "  -dmin 2 -dmax 14 -tsteps -1 -tmax 3600.0 -tout -1.0 -fbath data/tohoku_static/bath.nc -fdispl data/tohoku_static/displ.nc"
#	elif defined(_FLASH)
            write(arguments, '(A, A)') trim(arguments), "  -dmin 2 -dmax 14 -tsteps -1 -tmax 3600.0 -tout -1.0 -fbath data/tohoku_static/bath.nc -fdispl data/tohoku_static/displ.nc"
#    	elif defined(_NUMA)
            write(arguments, '(A, A)') trim(arguments), "  -dmin 2 -dmax 14 -tsteps -1 -tmax 5 -tout -1.0"
#    	elif defined(_TESTS)
            write(arguments, '(A, A)') trim(arguments), "  -dmin 2 -dmax 14 -tsteps -1 -tmax 5 -tout -1.0"
#    	elif defined(_GENERIC)
            write(arguments, '(A, A)') trim(arguments), "  -dmin 2 -dmax 14 -tsteps -1 -tmax 5 -tout -1.0"
#    	else
#           error No scenario selected!
#    	endif

        call kracken('samoa', arguments)

        !  get values
        l_help = lget('samoa_-help') .or. lget('samoa_h')
        l_version = lget('samoa_-version') .or. lget('samoa_v')

        config%i_min_depth = iget('samoa_dmin')
        config%i_max_depth = iget('samoa_dmax')
        config%i_max_time_steps = iget('samoa_tsteps')
        config%r_max_time = rget('samoa_tmax')
        config%r_output_time_step = rget('samoa_tout')
        config%l_log = lget('samoa_noprint')
        config%i_threads = iget('samoa_threads')
        config%i_sections_per_thread = iget('samoa_sections')
        config%i_asagi_mode = iget('samoa_asagihints')

#    	if defined(_DARCY)
            config%s_permeability_file = sget('samoa_fperm', 256)
            config%r_epsilon = rget('samoa_epsilon')
			config%r_rel_permeability = rget('samoa_k_rel')
			config%r_rho = rget('samoa_rho')
			config%r_p0 = rget('samoa_p0')
            config%i_lsolver = iget('samoa_lsolver')
#    	elif defined(_SWE)
            config%s_bathymetry_file = sget('samoa_fbath', 256)
            config%s_displacement_file = sget('samoa_fdispl', 256)
#    	elif defined(_FLASH)
            config%s_bathymetry_file = sget('samoa_fbath', 256)
            config%s_displacement_file = sget('samoa_fdispl', 256)
#       endif

        if (rank_MPI == 0) then
             _log_write(0, ' (" sam(oa)²: Space filling curves and Adaptive Meshes for Oceanic and Other Applications")')
            !if the version option was set was called, display program version
            if (l_version) then
                _log_write(0, '(" version ", I0, ".", I0, ".", I0)') 0, 5, 2
            end if

            !if the help option was set, display the list of arguments
            if (l_help) then
                _log_write(0, "")
                PRINT '(A)',            " Usage: samoa [--help | -h] | [--version | -v] | [-asagihints <value>] [-dstart <value>] [-dmin <value>] [-dmax <value>] [-tsteps <value>] [-tmax <value>] [-tout <value>] [-threads <value>] [-sections <value>] [-noprint]"
                PRINT '(A)',            ""
                PRINT '(A)',            " Arguments:"
                PRINT '(A, I0, ": ", A, A)',  " 	-asagihints <value>     ASAGI mode (0: default, 1: pass through, 2: no mpi, 3: no mpi + small cache, 4: large grid) (value: ", config%i_asagi_mode, trim(asagi_mode_to_char(config%i_asagi_mode)), ")"
                PRINT '(A, I0, A)',     " 	-dmin <value>           minimum grid depth (value: ", config%i_min_depth, ")"
                PRINT '(A, I0, A)',     "	-dmax <value>           maximum grid depth (value: ", config%i_max_depth, ")"
                PRINT '(A, I0, A)',     "	-tsteps <value>         maximum number of time steps, less than 0: not defined (value: ", config%i_max_time_steps, ")"
                PRINT '(A, ES8.1, A)',  "	-tmax <value>           maximum simulation time in seconds, less than 0: not defined (value: ", config%r_max_time, ")"
                PRINT '(A, ES8.1, A)',  "	-tout <value>           output time step in seconds, less than 0: not defined (value: ", config%r_output_time_step, ")"
                PRINT '(A, I0, A)',     "	-threads <value>        number of OpenMP threads (value: ", config%i_threads, ")"
                PRINT '(A, I0, A)',     "	-sections <value>       number of grid sections per OpenMP thread (value: ", config%i_sections_per_thread, ")"

#       	    if defined(_DARCY)
                    PRINT '(A, A, A)',  "	-fperm <value>          permeability template xyz(_*).nc (value: ", trim(config%s_permeability_file), ")"
                    PRINT '(A, ES8.1, A)',  "	-epsilon			    linear solver error bound (value: ", config%r_epsilon, ")"
                    PRINT '(A, ES8.1, A)',  "	-k_rel	                relative permeability of the entering fluid (value: ", config%r_rel_permeability, ")"
                    PRINT '(A, ES8.1, A)',  "	-rho				    fluid density (value: ", config%r_rho, ")"
                    PRINT '(A, ES8.1, A)',  "	-p0			            initial boundary pressure difference (value: ", config%r_p0, ")"
                    PRINT '(A, I0, ": ", A, A)',  "	-lsolver			    linear solver (0: Jacobi, 1: CG, 2: Pipelined CG) (value: ", config%i_lsolver, trim(lsolver_to_char(config%i_lsolver)), ")"
#         	    elif defined(_SWE)
                    PRINT '(A, A, A)',  "	-fbath <value>          bathymetry file (value: ", trim(config%s_bathymetry_file), ")"
                    PRINT '(A, A, A)',  "	-fdispl <value>         displacement file (value: ", trim(config%s_displacement_file), ")"
#         	    elif defined(_FLASH)
                    PRINT '(A, A, A)',  "	-fbath <value>          bathymetry file (value: ", trim(config%s_bathymetry_file), ")"
                    PRINT '(A, A, A)',  "	-fdispl <value>         displacement file (value: ", trim(config%s_displacement_file), ")"
#               endif

                PRINT '(A)',            "	-noprint                print log to file instead of console"
                PRINT '(A)',            "	--help, -h              display this help and exit"
                PRINT '(A)',            "	--version, -v           output version information and exit"
            end if

            _log_write(0, "")
        end if

        !stop if the version or help command were called
        if (l_help .or. l_version) then
            stop
        end if

        !init openmp threads
        call omp_set_num_threads(config%i_threads)

        !$omp parallel copyin(cfg)
            !nothing happens here, this parallel section is just meant for copying the data to all thread-private variables
            i_thread = 1 + omp_get_thread_num()
        !$omp end parallel
    end subroutine

    subroutine config_print(config)
        class(t_config)                         :: config
        character(64)                           :: lsolver_to_char(0:2) = ["Jacobi", "CG", "Pipelined CG"]
        character(64)                           :: asagi_mode_to_char(0:4) = ["default", "pass through", "no mpi", "no mpi + small cache", "large grid"]

#	    if defined(_TESTS)
            _log_write(0, '(" Scenario: Tests")')
#		elif defined (_HEAT_EQ)
            _log_write(0, '(" Scenario: Heat Equation")')
#		elif defined(_DARCY)
            _log_write(0, '(" Scenario: Darcy")')
#		elif defined(_SWE)
            _log_write(0, '(" Scenario: SWE")')
#		elif defined(_NUMA)
            _log_write(0, '(" Scenario: NUMA")')
#		elif defined(_GENERIC)
            _log_write(0, '(" Scenario: GENERIC")')
#		endif

#		if defined(_OPENMP)
#		    if defined(_OPENMP_TASKS)
                _log_write(0, '(" OpenMP: Yes, with tasks, threads: ", I0, ", procs: ", I0)') config%i_threads, omp_get_num_procs()
#		    else
                _log_write(0, '(" OpenMP: Yes, without tasks, threads: ", I0, ", procs: ", I0)') config%i_threads, omp_get_num_procs()
#		    endif
#		else
            _log_write(0, '(" OpenMP: No")')
#		endif

#		if defined(_MPI)
            _log_write(0, '(" MPI: Yes, ranks: ", I0)') size_MPI
#		else
            _log_write(0, '(" MPI: No")')
#		endif

#		if defined(_ASAGI)
#		    if defined(_ASAGI_NUMA)
                _log_write(0, '(" ASAGI: Yes, with NUMA support, mode: ", I0, ": ", A)') config%i_asagi_mode, trim(asagi_mode_to_char(config%i_asagi_mode))
#	        else
                _log_write(0, '(" ASAGI: Yes, without NUMA support, mode: ", I0, ": ", A)') config%i_asagi_mode, trim(asagi_mode_to_char(config%i_asagi_mode))
#		    endif
#		else
            _log_write(0, '(" ASAGI: No")')
#		endif

#		if defined(_DEBUG_LEVEL)
            _log_write(0, '(" Debug Level: ", I0)') _DEBUG_LEVEL
#		endif

#		if defined(_ASSERT)
            _log_write(0, '(" Assertions: Yes")')
#		else
            _log_write(0, '(" Assertions: No")')
#		endif

        _log_write(0, '(" Sections per thread: ", I0)') config%i_sections_per_thread
        _log_write(0, '(" Adaptivity: min depth: ", I0, ", max depth: ", I0)') config%i_min_depth, config%i_max_depth
        _log_write(0, '(" Simulation: max time steps: ", I0, ", max time: ", ES9.2, ", output step: ", ES9.2)'), config%i_max_time_steps, config%r_max_time, config%r_output_time_step

#		if defined(_DARCY)
            _log_write(0, '(" Scenario: permeability template: ", A)') trim(config%s_permeability_file)
            _log_write(0, '(" Scenario: linear solver error bound: ", ES8.1)') config%r_epsilon
            _log_write(0, '(" Scenario: relative permeability of the entering fluid: ", ES8.1)') config%r_rel_permeability
            _log_write(0, '(" Scenario: fluid density: ", ES8.1)') config%r_rho
            _log_write(0, '(" Scenario: initial boundary pressure difference: ", ES8.1)') config%r_p0
            _log_write(0, '(" Scenario: linear solver: ", I0, ": ", A)') config%i_lsolver, trim(lsolver_to_char(config%i_lsolver))
#		elif defined(_SWE)
            _log_write(0, '(" Scenario: bathymetry file: ", A, ", displacement file: ", A)') trim(config%s_bathymetry_file), trim(config%s_displacement_file)
#		elif defined(_FLASH)
            _log_write(0, '(" Scenario: bathymetry file: ", A, ", displacement file: ", A)') trim(config%s_bathymetry_file), trim(config%s_displacement_file)

#           if defined (_SWE_LF)
               _log_write(0, '(" Flux solver: ", A)') "Lax Friedrichs"
#           elif defined (_SWE_LLF)
                _log_write(0, '(" Flux solver: ", A)')  "Local Lax Friedrichs"
#           elif defined(_SWE_LF_BATH)
                _log_write(0, '(" Flux solver: ", A)')  "Lax Friedrichs + Bathymetry"
#           elif defined(_SWE_LLF_BATH)
                _log_write(0, '(" Flux solver: ", A)') "Local Lax Friedrichs + Bathymetry"
#           elif defined(_SWE_FWAVE)
                _log_write(0, '(" Flux solver: ", A)')  "FWave"
#           elif defined(_SWE_SSQ_FWAVE)
                _log_write(0, '(" Flux solver: ", A)')  "SSQ-FWave"
#           elif defined(_SWE_AUG_RIEMANN)
                _log_write(0, '(" Flux solver: ", A)')  "Augmented Riemann"
#           endif
#		endif

        _log_write(0, "")
    end subroutine
end module
