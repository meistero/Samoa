
! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE

#include "Compilation_control.f90"

#if defined(_ASAGI)
#	include "asagi.f90"
#endif

module config
	use M_kracken
	use Tools_log
    use Tools_openmp
    use Tools_mpi

#	if defined(_ASAGI)
		use asagi
#	endif

    implicit none

    private
    public cfg

    type t_config
        integer			                        :: i_threads			                            !< number of OpenMP threads
        integer (kind = selected_int_kind(1))   :: i_sections_per_thread						    !< number of sections per thread
        integer (kind = selected_int_kind(8))   :: i_max_time_steps				                    !< number of simulation time steps
        double precision                        :: r_max_time, r_output_time_step					!< maximum simulation time, output time step
        logical			                        :: l_log                                            !< if true, a log file is used
        integer (kind = selected_int_kind(1))   :: i_min_depth, i_max_depth			                !< minimum and maximum scenario depth
        integer			        	            :: i_asagi_mode			                		    !< ASAGI mode
        integer                                 :: i_ascii_width                                    !< width of the ascii output
        logical                                 :: l_ascii_output                                   !< ascii output on/off
        logical                                 :: l_timed_load                                     !< if true, load is estimated by timing, if false load is estimated by counting entities
        double precision                        :: r_cell_weight                                    !< cell weight for the count-based load estimate
        double precision                        :: r_boundary_weight                                !< boundary weight for the count-based load estimate
        logical                                 :: l_split_sections                                 !< if true, MPI load balancing may split sections, if false sections are treated as atomic units
        logical                                 :: l_serial_lb                                      !< if true, MPI load balancing is serialized, if false a distributed algorithm is used

        double precision                        :: scaling, offset(2)                               !< grid scaling and offset
        double precision                        :: courant_number                                   !< time step size relative to the CFL condition

#    	if defined(_DARCY)
            character(256)                      :: s_permeability_file                              !< permeability file
            character(256)                      :: s_porosity_file                                  !< porosity file
 			integer					 		    :: afh_permeability			                        !< asagi file handle to permeability data
 			integer					 		    :: afh_porosity			                            !< asagi file handle to porosity data
            integer			        	        :: i_lsolver			                		    !< linear solver
            integer			        	        :: i_CG_restart			                            !< CG restart interval
            logical                             :: l_lse_output                                     !< print out the linear equation system

			double precision			        :: r_epsilon				                        !< linear solver error bound
			double precision				    :: r_nu_w		                                    !< viscosity of the wetting phase
			double precision				    :: r_nu_n		                                    !< viscosity of the non-wetting phase
			double precision				    :: r_rho_w		                                    !< density of the wetting phase
			double precision				    :: r_rho_n		                                    !< density of the non-wetting phase
			double precision				    :: r_phi					                        !< porosity
			double precision				    :: r_p_in			                                !< injection well pressure
			double precision				    :: r_p_prod			                                !< production well pressure
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

    contains

    subroutine read_from_arguments(config)
        class(t_config), intent(inout)          :: config

        logical					                :: l_help, l_version
        integer          					    :: i, i_error
        character(1024)                         :: arguments
        character(64), parameter           		:: lsolver_to_char(0:3) = [character(64) :: "Jacobi", "CG", "Pipelined CG", "Pipelined CG (unstable)"]
        character(64), parameter             	:: asagi_mode_to_char(0:4) = [character(64) :: "default", "pass through", "no mpi", "no mpi + small cache", "large grid"]

        !define default command arguments and default values for all scenarios
        write(arguments, '(A)') "-v .false. --version .false. -h .false. --help .false."
        write(arguments, '(A, A)') trim(arguments),   " -lbtime .false. -lbsplit .false. -lbserial .false. -lbcellweight 1.0d0 -lbbndweight 0.0d0"
        write(arguments, '(A, A)') trim(arguments),  " -asagihints 2 -asciiout_width 60 -asciiout .false. -noprint .false. -sections 4"
        write(arguments, '(A, A, I0)') trim(arguments), " -threads ", omp_get_max_threads()

        !define additional command arguments and default values depending on the choice of the scenario
#    	if defined(_DARCY)
            write(arguments, '(A, A)') trim(arguments), " -dmin 1 -dmax 14 -tsteps -1 -courant 1.0d0 -tmax 2.0d1 -tout -1.0d0 -fperm data/darcy_benchmark/perm_1024.nc -fpor data/darcy_benchmark/perm_1024.nc -p_in 10.0d3 -p_prod 4.0d3 -epsilon 1.0d-5 -phi 0.2d0 -rho_w 1.0d0 -rho_n 0.67d0 -nu_w 0.3d0 -nu_n 1.01d0 -lsolver 2 -cg_restart 256 -lse_output .false."
#    	elif defined(_HEAT_EQ)
            write(arguments, '(A, A)') trim(arguments), " -dmin 1 -dmax 16 -tsteps -1 -tmax 1.0d0 -tout -1.0d0"
#    	elif defined(_SWE)
            write(arguments, '(A, A)') trim(arguments), " -dmin 2 -dmax 14 -tsteps -1 -courant 0.45d0 -tmax 3600.0d0 -tout -1.0d0 -fbath data/tohoku_static/bath.nc -fdispl data/tohoku_static/displ.nc"
#	    elif defined(_FLASH)
            write(arguments, '(A, A)') trim(arguments), " -dmin 2 -dmax 14 -tsteps -1 -courant 0.45d0 -tmax 3600.0d0 -tout -1.0d0 -fbath data/tohoku_static/bath.nc -fdispl data/tohoku_static/displ.nc"
#    	elif defined(_NUMA)
            write(arguments, '(A, A)') trim(arguments), " -dmin 2 -dmax 14 -tsteps -1 -courant 0.45d0 -tmax 5 -tout -1.0d0"
#    	elif defined(_TESTS)
            write(arguments, '(A, A)') trim(arguments), " -dmin 2 -dmax 14 -tsteps -1 -tmax 5 -tout -1.0d0"
#    	elif defined(_GENERIC)
            write(arguments, '(A, A)') trim(arguments), " -dmin 2 -dmax 14 -tsteps -1 -tmax 5 -tout -1.0d0"
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
        config%l_timed_load = lget('samoa_lbtime')
        config%r_cell_weight = rget('samoa_lbcellweight')
        config%r_boundary_weight = rget('samoa_lbbndweight')
        config%l_split_sections = lget('samoa_lbsplit')
        config%l_serial_lb = lget('samoa_lbserial')
        config%i_sections_per_thread = iget('samoa_sections')
        config%i_asagi_mode = iget('samoa_asagihints')
        config%l_ascii_output = lget('samoa_asciiout')
        config%i_ascii_width = iget('samoa_asciiout_width')
        config%courant_number = rget('samoa_courant')

#    	if defined(_DARCY)
            config%s_permeability_file = sget('samoa_fperm', 256)
            config%s_porosity_file = sget('samoa_fpor', 256)
            config%r_epsilon = rget('samoa_epsilon')
			config%r_nu_w = rget('samoa_nu_w')
			config%r_nu_n = rget('samoa_nu_n')
			config%r_rho_w = rget('samoa_rho_w')
			config%r_rho_n = rget('samoa_rho_n')
			config%r_phi = rget('samoa_phi')
			config%r_p_in = rget('samoa_p_in')
			config%r_p_prod = rget('samoa_p_prod')
            config%i_lsolver = iget('samoa_lsolver')
            config%i_CG_restart = iget('samoa_cg_restart')
            config%l_lse_output = lget('samoa_lse_output')
#    	elif defined(_SWE) || defined(_FLASH)
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
                PRINT '()'
                PRINT '(A)',            " Usage: samoa [--help | -h] | [--version | -v] | [OPTION...]"
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
                PRINT '(A, L, A)',      "	-lbtime                 if true, load is estimated by time measurements, if false load is estimated by cell count (value: ", config%l_timed_load, ")"
                PRINT '(A, L, A)',      "	-lbsplit                if true, MPI load balancing may split sections, if false sections are treated as atomic units (value: ", config%l_split_sections, ")"
                PRINT '(A, L, A)',      "	-lbserial               if true, MPI load balancing is serialized, if false a distributed algorithm is used (value: ", config%l_serial_lb, ")"
                PRINT '(A, F0.3, A)',  "	-lbcellweight           cell weight for the count-based load estimate (value: ", config%r_cell_weight, ")"
                PRINT '(A, F0.3, A)',  "	-lbbndweight            boundary weight for the count-based load estimate (value: ", config%r_boundary_weight, ")"
                PRINT '(A, F0.3, A)',  "	-courant                time step size relative to the CFL condition (value: ", config%courant_number, ")"

#       	    if defined(_DARCY)
                    PRINT '(A, A, A)',  "	-fperm <value>          permeability file (value: ", trim(config%s_permeability_file), ")"
                    PRINT '(A, A, A)',  "	-fpor <value>           porosity file (value: ", trim(config%s_porosity_file), ")"
                    PRINT '(A, ES8.1, A)',  "	-epsilon			    linear solver error bound (value: ", config%r_epsilon, ")"
                    PRINT '(A, ES8.1, A)',  "	-nu_w	                viscosity of the wetting phase (value: ", config%r_nu_w, ")"
                    PRINT '(A, ES8.1, A)',  "	-nu_n	                viscosity of the non-wetting phase (value: ", config%r_nu_n, ")"
                    PRINT '(A, ES8.1, A)',  "	-rho_w	                density of the wetting phase (value: ", config%r_rho_w, ")"
                    PRINT '(A, ES8.1, A)',  "	-rho_n	                density of the non-wetting phase (value: ", config%r_rho_n, ")"
                    PRINT '(A, ES8.1, A)',  "	-phi				    porosity (value: ", config%r_phi, ")"
                    PRINT '(A, ES8.1, A)',  "	-p_in			        injection well pressure (value: ", config%r_p_in, ")"
                    PRINT '(A, ES8.1, A)',  "	-p_prod			        production well pressure (value: ", config%r_p_prod, ")"
                    PRINT '(A, I0, ": ", A, A)',  "	-lsolver			    linear solver (0: Jacobi, 1: CG, 2: Pipelined CG) (value: ", config%i_lsolver, trim(lsolver_to_char(config%i_lsolver)), ")"
                    PRINT '(A, I0, A)',     "	-cg_restart			    CG restart interval (value: ", config%i_CG_restart, ")"
                    PRINT '(A, L, A)',     "	-lse_output             enable LSE output (value: ", config%l_lse_output, ")"
#         	    elif defined(_SWE)
                    PRINT '(A, L, A)',  "	-asciiout               turns on ascii output (value: ", config%l_ascii_output, ")"
                    PRINT '(A, I0, A)', "	-asciiout_width <value> width of ascii output (value: ", config%i_ascii_width, ")"
                    PRINT '(A, A, A)',  "	-fbath <value>          bathymetry file (value: ", trim(config%s_bathymetry_file), ")"
                    PRINT '(A, A, A)',  "	-fdispl <value>         displacement file (value: ", trim(config%s_displacement_file), ")"
#         	    elif defined(_FLASH)
                    PRINT '(A, A, A)',  "	-fbath <value>          bathymetry file (value: ", trim(config%s_bathymetry_file), ")"
                    PRINT '(A, A, A)',  "	-fdispl <value>         displacement file (value: ", trim(config%s_displacement_file), ")"
#               endif

                PRINT '(A, L, A)',      "	-noprint                print log to file instead of console (value: ", config%l_log, ")"
                PRINT '(A)',            "	--help, -h              display this help and exit"
                PRINT '(A)',            "	--version, -v           output version information and exit"
            end if
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

    function logical_to_char(l) result(s)
        logical, intent(in) :: l
        character(3)        :: s

        s = merge("Yes", " No", l)
    end function

    subroutine config_print(config)
        class(t_config)                         :: config
        character(64), parameter           		:: lsolver_to_char(0:3) = [character(64) :: "Jacobi", "CG", "Pipelined CG", "Pipelined CG (unstable)"]
        character(64), parameter             	:: asagi_mode_to_char(0:4) = [character(64) :: "default", "pass through", "no mpi", "no mpi + small cache", "large grid"]

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

#       if defined(_SINGLE_PRECISION)
            _log_write(0, '(" Precision: Single")')
#       elif defined(_DOUBLE_PRECISION)
            _log_write(0, '(" Precision: Double")')
#       elif defined(_QUAD_PRECISION)
            _log_write(0, '(" Precision: Quad")')
#       else
#           error "Invalid floating point precision!"
#       endif

#       if defined(__GFORTRAN__)
            _log_write(0, '(" Compiler: GNU")')
#       else
            _log_write(0, '(" Compiler: Intel")')
#       endif

        _log_write(0, '(" Sections per thread: ", I0)') config%i_sections_per_thread
        _log_write(0, '(" Adaptivity: min depth: ", I0, ", max depth: ", I0)') config%i_min_depth, config%i_max_depth

        _log_write(0, '(" Load balancing: timed load estimate: ", A, ", split sections: ", A, ", serial: ", A)') logical_to_char(config%l_timed_load), logical_to_char(config%l_split_sections), logical_to_char(config%l_serial_lb)
        _log_write(0, '(" Load balancing: cell weight: ", F0.2, ", boundary weight: ", F0.2)') config%r_cell_weight, config%r_boundary_weight

        _log_write(0, '(" Scenario: max time steps: ", I0, ", max time: ", ES9.2, ", output step: ", ES9.2)'), config%i_max_time_steps, config%r_max_time, config%r_output_time_step
        _log_write(0, '(" Scenario: courant number: ", F0.3)'), config%courant_number

#		if defined(_DARCY)
            _log_write(0, '(" Darcy: permeability file: ", A)') trim(config%s_permeability_file)
            _log_write(0, '(" Darcy: linear solver error bound: ", ES8.1)') config%r_epsilon
            _log_write(0, '(" Darcy: vicosities: wetting phase: ", ES8.1, ", non-wetting phase: ", ES8.1)') config%r_nu_w, config%r_nu_n
            _log_write(0, '(" Darcy: densities: wetting phase: ", ES8.1, ", non-wetting phase: ", ES8.1)') config%r_rho_w, config%r_rho_n
            _log_write(0, '(" Darcy: porosity: ", ES8.1)') config%r_phi
            _log_write(0, '(" Darcy: injection well pressure: ", ES8.1, ", production well pressure: ", ES8.1)') config%r_p_in, config%r_p_prod
            _log_write(0, '(" Darcy: linear solver: ", I0, ": ", A)') config%i_lsolver, trim(lsolver_to_char(config%i_lsolver))
            _log_write(0, '(" Darcy: CG restart interval: ", I0)') config%i_CG_restart
#		elif defined(_FLASH)
            _log_write(0, '(" Flash: bathymetry file: ", A, ", displacement file: ", A)') trim(config%s_bathymetry_file), trim(config%s_displacement_file)
#		elif defined(_SWE)
            _log_write(0, '(" SWE: bathymetry file: ", A, ", displacement file: ", A)') trim(config%s_bathymetry_file), trim(config%s_displacement_file)

            if (config%l_ascii_output) then
                _log_write(0, '(" SWE: Ascii Output: Yes, width: ", I0)') config%i_ascii_width
            else
                _log_write(0, '(" SWE: Ascii Output: No")')
            end if

#           if defined (_SWE_LF)
               _log_write(0, '(" SWE: Flux solver: ", A)') "Lax Friedrichs"
#           elif defined (_SWE_LLF)
                _log_write(0, '(" SWE: Flux solver: ", A)')  "Local Lax Friedrichs"
#           elif defined(_SWE_LF_BATH)
                _log_write(0, '(" SWE: Flux solver: ", A)')  "Lax Friedrichs + Bathymetry"
#           elif defined(_SWE_LLF_BATH)
                _log_write(0, '(" SWE: Flux solver: ", A)') "Local Lax Friedrichs + Bathymetry"
#           elif defined(_SWE_FWAVE)
                _log_write(0, '(" SWE: Flux solver: ", A)')  "FWave"
#           elif defined(_SWE_SSQ_FWAVE)
                _log_write(0, '(" SWE: Flux solver: ", A)')  "SSQ-FWave"
#           elif defined(_SWE_AUG_RIEMANN)
                _log_write(0, '(" SWE: Flux solver: ", A)')  "Augmented Riemann"
#           endif
#		endif

        _log_write(0, "()")
    end subroutine
end module
