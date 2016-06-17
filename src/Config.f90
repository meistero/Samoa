
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
        integer (kind = selected_int_kind(8))   :: i_sections_per_thread						    !< number of sections per thread
        double precision                        :: r_max_time					                    !< maximum simulation time
        integer (kind = selected_int_kind(8))   :: i_max_time_steps				                    !< number of simulation time steps
        double precision                        :: r_output_time_step					            !< grid output time step
        integer                                 :: i_output_time_steps					            !< grid output time step
        integer                                 :: i_stats_phases					                !< number of times intermediate stats should be printed during time steps
        logical			                        :: l_log                                            !< if true, a log file is used
        integer (kind = selected_int_kind(1))   :: i_min_depth, i_max_depth, i_start_depth			!< minimum, maximum and start scenario depth
        integer			        	            :: i_asagi_mode			                		    !< ASAGI mode
        logical                                 :: l_timed_load                                     !< if true, load is estimated by timing, if false load is estimated by counting entities
        double precision                        :: r_cell_weight                                    !< cell weight for the count-based load estimate
        double precision                        :: r_boundary_weight                                !< boundary weight for the count-based load estimate
        logical                                 :: l_split_sections                                 !< if true, MPI load balancing may split sections, if false sections are treated as atomic units
        logical                                 :: l_serial_lb                                      !< if true, MPI load balancing is serialized, if false a distributed algorithm is used
        double precision                        :: r_adapt_time_step					            !< grid output time step
        integer			        	            :: i_adapt_time_steps			                            !< number of time steps between each linear solver solution

	    logical 				                :: l_gridoutput			                            !< grid output on/off
	    character(256)				            :: output_dir			                            !< output directory

        double precision                        :: courant_number                                   !< time step size relative to the CFL condition

#    	if defined(_DARCY)
            double precision                    :: scaling, offset(3)                               !< grid scaling and offset of the computational domain
            double precision                    :: x_min(3), x_max(3)                               !< lower and upper bounds of the data domain
            double precision                    :: dx(3), dz                                        !< voxel size of the source data

#           if defined(_ASAGI)
                character(256)                  :: s_permeability_file                              !< permeability file
                character(256)                  :: s_porosity_file                                  !< porosity file

                integer					 	    :: afh_permeability_X			                    !< asagi file handle to X-axis permeability data
                integer					 	    :: afh_permeability_Y			                    !< asagi file handle to Y-axis permeability data
                integer					 	    :: afh_permeability_Z			                    !< asagi file handle to Z-axis permeability data
                integer					 	    :: afh_porosity			                            !< asagi file handle to porosity data
#           endif

            integer			        	        :: i_lsolver			                		    !< linear solver
            integer			        	        :: i_max_iterations			                		!< maximum iterations of the linear solver
            integer			        	        :: i_CG_restart			                            !< CG restart interval
            double precision                    :: r_solver_time_step					            !< grid output time step
            integer			        	        :: i_solver_time_steps			                            !< number of time steps between each linear solver solution
            logical                             :: l_lse_output                                     !< print out the linear equation system

			double precision			        :: r_epsilon				                        !< linear solver error bound
			double precision				    :: r_nu_w, r_nu_w_SI		                        !< viscosity of the wetting phase [Pa s]
			double precision				    :: r_nu_n, r_nu_n_SI		                        !< viscosity of the non-wetting phase [Pa s]
			double precision				    :: r_rho_w, r_rho_w_SI	                            !< density of the wetting phase [kg/m^3]
			double precision				    :: r_rho_n, r_rho_n_SI		                        !< density of the non-wetting phase [kg/m^3]
			double precision			        :: S_wr		                                        !< residual saturation of the wetting phase
			double precision			        :: S_nr		                                        !< residual saturation of the non-wetting phase
			double precision				    :: r_p_in, r_p_in_AU			                    !< injection well pressure [psi]
			double precision				    :: r_p_prod, r_p_prod_AU			                !< production well pressure [psi]
			double precision				    :: r_well_radius, r_well_radius_AU		            !< well radius for injection and production wells [in]
			double precision				    :: r_inflow, r_inflow_AU	                        !< injection well inflow [bbl/d]

            double precision			        :: r_pos_in(2, _DARCY_INJECTOR_WELLS)				!< injector positions [m]
            double precision			        :: r_pos_prod(2, _DARCY_PRODUCER_WELLS)				!< producer positions [m]

			double precision			        :: g(3), g_SI(3)				                    !< gravity vector [m/s^2]

			double precision			        :: p_refinement_threshold				            !< pressure refinement threshold
			double precision			        :: S_refinement_threshold				            !< saturation refinement threshold

            logical 				            :: l_well_output			                        !< well output on/off
#    	elif defined(_SWE)
            double precision                    :: scaling, offset(2)                               !< grid scaling and offset of the computational domain

#           if defined(_ASAGI)
                character(256)                  :: s_bathymetry_file                                !< bathymetry file
                character(256)                  :: s_displacement_file                              !< displacement file
                integer					 		:: afh_displacement			                        !< asagi file handle to displacement data
                integer					 		:: afh_bathymetry			                        !< asagi file handle to bathymetry da
#           endif

            double precision                    :: t_min_eq, t_max_eq						        !< earthquake start and end time [s]
            double precision                    :: dt_eq                                            !< earthquake time step [s]
            double precision                    :: dry_tolerance                                    !< dry tolerance [m]

            logical                             :: l_ascii_output                                   !< ascii output on/off
            integer                             :: i_ascii_width                                    !< width of the ascii output
            logical					            :: l_pointoutput                                    !< test points output on/off
            character(512)				        :: s_testpoints			                            !< test points input string
            double precision, pointer		    :: r_testpoints(:,:)		                        !< test points array
#    	elif defined(_FLASH)
            double precision                    :: scaling, offset(2)                               !< grid scaling and offset of the computational domain

            character(256)                      :: s_bathymetry_file                                !< bathymetry file
            character(256)                      :: s_displacement_file                              !< displacement file
 			integer					 		    :: afh_displacement			                        !< asagi file handle to displacement data
 			integer					 		    :: afh_bathymetry			                        !< asagi file handle to bathymetry data

            double precision                    :: dry_tolerance                                    !< dry tolerance
#       endif

        contains

        procedure, pass :: read_from_program_arguments
        procedure, pass :: print => config_print

        procedure, private, pass :: define_arguments
        procedure, private, pass :: parse_program_arguments
        procedure, private, pass :: parse_string
        procedure, private, pass :: update
    end type

    type(t_config) :: cfg
    !$omp threadprivate(cfg)

    contains

    subroutine read_from_program_arguments(config)
        class(t_config), intent(inout)          :: config

        character(1024)    :: arguments

        call config%define_arguments(arguments)
        call config%parse_program_arguments(arguments)
    end subroutine

    subroutine define_arguments(config, arguments)
        class(t_config), intent(inout)  :: config
        character(*), intent(inout)     :: arguments

        !define default command arguments and default values for all scenarios

        write(arguments, '(A)') "-v .false. --version .false. -h .false. --help .false."
        write(arguments, '(A, A)') trim(arguments),   " -lbtime .false. -lbsplit .false. -lbserial .false. -lbcellweight 1.0d0 -lbbndweight 0.0d0"
        write(arguments, '(A, A)') trim(arguments),  " -asagihints 2 -phases 1 -tadapt -1.0 -nadapt 1 -asciioutput_width 60 -output_dir output -asciioutput .false. -xmloutput .false. -stestpoints '' -noprint .false. -sections 4"
        write(arguments, '(A, A, I0)') trim(arguments), " -threads ", omp_get_max_threads()

        !define additional command arguments and default values depending on the choice of the scenario
#    	if defined(_DARCY)
            write(arguments, '(A, A)') trim(arguments), " -dmin 0 -dmax 14 -dstart 0 -courant 0.5d0 " // &
            "-nout -1 -tout -1.0 -fperm data/darcy_five_spot/spe_perm_renamed.nc -fpor data/darcy_five_spot/spe_phi_renamed.nc "  // &
            "-p_in 10.0d3 -p_prod 4.0d3 -epsilon 1.0d-4 -rho_w 1025.18d0 -rho_n 848.98d0 -nu_w 0.3d-3 -nu_n 3.0d-3 -lsolver 2 " // &
            "-max_iter -1 -tsolver -1.0 -nsolver 1 -cg_restart 256 -lseoutput .false. -welloutput .true. "

#           if (_DARCY_LAYERS > 0)
                write(arguments, '(A, A)') trim(arguments), " -p_ref_th 1.0d0 -S_ref_th 1.0d0 "
#           else
                write(arguments, '(A, A)') trim(arguments), " -p_ref_th 0.75d0 -S_ref_th 0.3d0"
#           endif

#           if defined(_ASAGI)
                write(arguments, '(A, A)') trim(arguments), " -nmax -1 -tmax 172.8d6 -g_x 0.0d0 -g_y 0.0d0 -g_z -9.80665d0 " // &
                "-inflow 5000.0d0 -well_radius 5.0d0 -S_wr 0.2 -S_nr 0.2"
#           else
                write(arguments, '(A, A)') trim(arguments), " -nmax -1 -tmax 2.0d1 -g_x 9.80665d0 -g_y 0.0d0 -g_z 0.0d0 " // &
                "-inflow 0.5434396505 -well_radius 5.0d0  -S_wr 0.0 -S_nr 0.0"
#           endif
#    	elif defined(_HEAT_EQ)
            write(arguments, '(A, A)') trim(arguments), " -dmin 0 -dmax 16 -dstart 0 -nmax -1 -tmax 1.0d0 -nout -1 -tout -1.0d0"
#    	elif defined(_SWE)
            write(arguments, '(A, A)') trim(arguments), " -dmin 0 -dmax 14 -dstart 0 -courant 0.45d0 -nmax -1 -tmax 3600.0d0 -nout -1 -tout -1.0d0 -drytolerance 0.01d0 -fbath data/tohoku_static/bath.nc -fdispl data/tohoku_static/displ.nc"
#	    elif defined(_FLASH)
            write(arguments, '(A, A)') trim(arguments), " -dmin 0 -dmax 14 -dstart 0 -courant 0.45d0 -nmax -1 -tmax 3600.0d0 -nout -1 -tout -1.0d0 -drytolerance 0.01d0 -fbath data/tohoku_static/bath.nc -fdispl data/tohoku_static/displ.nc"
#    	elif defined(_NUMA)
            write(arguments, '(A, A)') trim(arguments), " -dmin 0 -dmax 14 -dstart 0 -courant 0.45d -nmax -1 -tmax 5 -nout -1 -tout -1.0d0"
#    	elif defined(_TESTS)
            write(arguments, '(A, A)') trim(arguments), " -dmin 0 -dmax 14 -dstart 0 -nmax -1 -tmax 5 -nout -1 -tout -1.0d0"
#    	elif defined(_GENERIC)
            write(arguments, '(A, A)') trim(arguments), " -dmin 0 -dmax 14 -dstart 0 -nmax -1 -tmax 5 -nout -1 -tout -1.0d0"
#    	else
#           error No scenario selected!
#    	endif
    end subroutine

    subroutine parse_string(config, string)
        class(t_config), intent(inout)  :: config
        character(*), intent(inout)     :: string

        call parse('samoa', string, "no_add")
        call config%update()
    end subroutine

    subroutine parse_program_arguments(config, arguments)
        class(t_config), intent(inout)  :: config
        character(*), intent(inout)     :: arguments

        call kracken('samoa', arguments)
        call config%update()
    end subroutine

    subroutine update(config)
        class(t_config), intent(inout)  :: config

        logical					        :: l_help, l_version
        character(64), parameter        :: lsolver_to_char(0:3) = [character(64) :: "Jacobi", "CG", "Pipelined CG", "Pipelined CG (unstable)"]
        character(64), parameter        :: asagi_mode_to_char(0:4) = [character(64) :: "default", "pass through", "no mpi", "no mpi + small cache", "large grid"]

        !  get values
        l_help = lget('samoa_-help') .or. lget('samoa_h')
        l_version = lget('samoa_-version') .or. lget('samoa_v')

        config%i_min_depth = int(iget('samoa_dmin'), selected_int_kind(1))
        config%i_max_depth = int(iget('samoa_dmax'), selected_int_kind(1))
        config%i_start_depth = int(iget('samoa_dstart'), selected_int_kind(1))
        config%i_max_time_steps = iget('samoa_nmax')
        config%r_max_time = rget('samoa_tmax')
        config%i_output_time_steps = iget('samoa_nout')
        config%r_output_time_step = rget('samoa_tout')
        config%i_stats_phases = iget('samoa_phases')
        config%l_log = lget('samoa_noprint')
        config%i_threads = iget('samoa_threads')
        config%l_timed_load = lget('samoa_lbtime')
        config%r_cell_weight = rget('samoa_lbcellweight')
        config%r_boundary_weight = rget('samoa_lbbndweight')
        config%l_split_sections = lget('samoa_lbsplit')
        config%l_serial_lb = lget('samoa_lbserial')
        config%i_sections_per_thread = iget('samoa_sections')
        config%i_asagi_mode = iget('samoa_asagihints')
        config%i_adapt_time_steps = iget('samoa_nadapt')
        config%r_adapt_time_step = iget('samoa_tadapt')
        config%courant_number = rget('samoa_courant')
        config%l_gridoutput = lget('samoa_xmloutput')
        config%output_dir = sget('samoa_output_dir', 256)

#    	if defined(_DARCY)
#		    if defined(_ASAGI)
                config%s_permeability_file = sget('samoa_fperm', 256)
                config%s_porosity_file = sget('samoa_fpor', 256)
#           endif

            config%r_epsilon = rget('samoa_epsilon')
			config%S_wr = rget('samoa_S_wr')
			config%S_nr = rget('samoa_S_nr')
            config%i_max_iterations = iget('samoa_max_iter')
            config%i_solver_time_steps = iget('samoa_nsolver')
            config%r_solver_time_step = iget('samoa_tsolver')
            config%i_lsolver = iget('samoa_lsolver')
            config%i_CG_restart = iget('samoa_cg_restart')

			config%r_nu_w_SI = rget('samoa_nu_w')
			config%r_nu_n_SI = rget('samoa_nu_n')
			config%r_rho_w_SI = rget('samoa_rho_w')
			config%r_rho_n_SI = rget('samoa_rho_n')
            config%g_SI = [rget('samoa_g_x'), rget('samoa_g_y'), rget('samoa_g_z')]

			config%r_p_in_AU = rget('samoa_p_in')
			config%r_p_prod_AU = rget('samoa_p_prod')
            config%r_well_radius_AU = rget('samoa_well_radius')
            config%r_inflow_AU = rget('samoa_inflow')

            config%p_refinement_threshold = rget('samoa_p_ref_th')
            config%S_refinement_threshold = rget('samoa_S_ref_th')

            config%l_lse_output = lget('samoa_lseoutput')
            config%l_well_output = lget('samoa_welloutput')
#    	elif defined(_SWE) || defined(_FLASH)
#		    if defined(_ASAGI)
                config%s_bathymetry_file = sget('samoa_fbath', 256)
                config%s_displacement_file = sget('samoa_fdispl', 256)
#           endif

            config%dry_tolerance = rget('samoa_drytolerance')

            config%l_ascii_output = lget('samoa_asciioutput')
            config%i_ascii_width = iget('samoa_asciioutput_width')
            config%s_testpoints = sget('samoa_stestpoints', 512)

            if (len(trim(config%s_testpoints)) .ne. 2) then
                config%l_pointoutput = .true.
                call parse_testpoints(config)
            else
                config%l_pointoutput = .false.
            end if
#       endif

        if (rank_MPI == 0) then
             _log_write(0, ' (" sam(oa)²: Space filling curves and Adaptive Meshes for Oceanic and Other Applications")')
            !if the version option was set was called, display program version
            if (l_version) then
                _log_write(0, '(" version ", I0, ".", I0, ".", I0)') 0, 8, 0
            end if

            !if the help option was set, display the list of arguments
            if (l_help) then
                PRINT '()'
                PRINT '(A)',            " Usage: samoa [--help | -h] | [--version | -v] | [OPTION...]"
                PRINT '(A)',            ""
                PRINT '(A)',            " General arguments:"
                PRINT '(A, I0, ": ", A, A)',  " 	-asagihints <value>     ASAGI mode (0: default, 1: pass through, 2: no mpi, 3: no mpi + small cache, 4: large grid) (value: ", config%i_asagi_mode, trim(asagi_mode_to_char(config%i_asagi_mode)), ")"
                PRINT '(A, I0, A)',     " 	-dmin <value>           minimum grid depth (value: ", config%i_min_depth, ")"
                PRINT '(A, I0, A)',     "	-dmax <value>           maximum grid depth (value: ", config%i_max_depth, ")"
                PRINT '(A, I0, A)',     "	-dstart <value>         start grid depth (value: ", config%i_start_depth, ")"
                PRINT '(A, I0, A)',     "	-nmax <value>           maximum number of time steps, less than 0: disabled (value: ", config%i_max_time_steps, ")"
                PRINT '(A, ES8.1, A)',  "	-tmax <value>           maximum simulation time in seconds, less than 0: disabled (value: ", config%r_max_time, ")"
                PRINT '(A, I0, A)',     "	-nout <value>           output time step interval, less than 1: disabled (value: ", config%i_output_time_steps, ")"
                PRINT '(A, ES8.1, A)',  "	-tout <value>           output time step in seconds, less than 0: disabled (value: ", config%r_output_time_step, ")"
                PRINT '(A, I0, A)',     "	-phases <value>         number of times intermediate stats should be printed during time steps (value: ", config%i_stats_phases, ")"
                PRINT '(A, I0, A)',     "	-threads <value>        number of OpenMP threads (value: ", config%i_threads, ")"
                PRINT '(A, I0, A)',     "	-sections <value>       number of grid sections per OpenMP thread (value: ", config%i_sections_per_thread, ")"
                PRINT '(A, L, A)',      "	-lbtime                 if true, load is estimated by time measurements, if false load is estimated by cell count (value: ", config%l_timed_load, ")"
                PRINT '(A, L, A)',      "	-lbsplit                if true, MPI load balancing may split sections, if false sections are treated as atomic units (value: ", config%l_split_sections, ")"
                PRINT '(A, L, A)',      "	-lbserial               if true, MPI load balancing is serialized, if false a distributed algorithm is used (value: ", config%l_serial_lb, ")"
                PRINT '(A, F0.3, A)',  "	-lbcellweight           cell weight for the count-based load estimate (value: ", config%r_cell_weight, ")"
                PRINT '(A, F0.3, A)',  "	-lbbndweight            boundary weight for the count-based load estimate (value: ", config%r_boundary_weight, ")"
                PRINT '(A, I0, A)',    "	-nadapt <value>         remeshing time step interval, less than 0: disabled (value: ", config%i_adapt_time_steps, ")"
                !This argument is not supported for now
                !PRINT '(A, I0, A)',   "	-tadapt                 remeshing time step in seconds, less than 0: disabled (value: ", config%r_adapt_time_step, ")"
                PRINT '(A, F0.3, A)',  "	-courant <value>        time step size relative to the CFL condition (value: ", config%courant_number, ")"
        		PRINT '(A, L, A)',     "	-xmloutput              [-tout required] turns on grid output (value: ", config%l_gridoutput, ")"
        		PRINT '(A, A, A)',     "	-output_dir <value>     output directory (value: ", trim(config%output_dir), ")"
                PRINT '(A, L, A)',      "	-noprint                print log to file instead of console (value: ", config%l_log, ")"
                PRINT '(A)',            "	--help, -h              display this help and exit"
                PRINT '(A)',            "	--version, -v           output version information and exit"

                PRINT '(A)',            ""
                PRINT '(A)',            " Scenario specific arguments:"
#       	    if defined(_DARCY)
#    	            if defined(_ASAGI)
                        PRINT '(A, A, A)',  "	-fperm <value>          permeability file (value: ", trim(config%s_permeability_file), ")"
                        PRINT '(A, A, A)',  "	-fpor <value>           porosity file (value: ", trim(config%s_porosity_file), ")"
#                   endif

                    PRINT '(A, ES8.1, A)',  "	-nu_w	                viscosity of the wetting phase (value: ", config%r_nu_w_SI, " Pa s)"
                    PRINT '(A, ES8.1, A)',  "	-nu_n	                viscosity of the non-wetting phase (value: ", config%r_nu_n_SI, " Pa s)"
                    PRINT '(A, ES8.1, A)',  "	-rho_w	                density of the wetting phase (value: ", config%r_rho_w_SI, " kg / m^3)"
                    PRINT '(A, ES8.1, A)',  "	-rho_n	                density of the non-wetting phase (value: ", config%r_rho_n_SI, " kg / m^3)"
                    PRINT '(A, ES8.1, A)',  "	-S_wr	                residual saturation of the wetting phase (value: ", config%S_wr, ")"
                    PRINT '(A, ES8.1, A)',  "	-S_nr	                residual saturation of the non-wetting phase (value: ", config%S_nr, ")"
                    PRINT '(A, ES8.1, A)',  "	-p_in                   injection well pressure (value: ", config%r_p_in_AU, " psi)"
                    PRINT '(A, ES8.1, A)',  "	-p_prod                 production well pressure (value: ", config%r_p_prod_AU, " psi)"
                    PRINT '(A, ES9.2, A)',  "	-inflow                 inflow condition (value: ", config%r_inflow_AU, " BBL/d)"
                    PRINT '(A, 3(ES9.2, X), A)',  "	-g_x -g_y -g_z          gravity vector (value: (", config%g_SI, ") m/s^2)"
                    PRINT '(A, ES8.1, A)',  "	-well_radius            injection and production well radius (value: ", config%r_well_radius_AU, " inch)"
                    PRINT '(A, I0, ": ", A, A)',  "	-lsolver                linear solver (0: Jacobi, 1: CG, 2: Pipelined CG) (value: ", config%i_lsolver, trim(lsolver_to_char(config%i_lsolver)), ")"
                    PRINT '(A, ES8.1, A)',  "	-epsilon                linear solver error bound (value: ", config%r_epsilon, ")"
                    PRINT '(A, I0, A)',        "	-max_iter               maximum iterations of the linear solver, less than 1: disabled (value: ", config%i_max_iterations, ")"
                    PRINT '(A, I0, A)',        "	-nsolver            linear solver time step interval, less than 1: disabled (value: ", config%i_solver_time_steps, ")"
                    !This argument is not supported for now
                    !PRINT '(A, I0, A)',        "	-tsolver            linear solver time step in seconds, less than 0: disabled (value: ", config%r_solver_time_step, ")"
                    PRINT '(A, I0, A)',     "	-cg_restart             CG restart interval (value: ", config%i_CG_restart, ")"
                    PRINT '(A, L, A)',     "	-lseoutput              enable LSE output (value: ", config%l_lse_output, ")"
                    PRINT '(A, L, A)',     "	-welloutput             enable well data output (value: ", config%l_well_output, ")"
                    PRINT '(A, ES9.2, A)',  "	-p_ref_th               pressure refinement threshold (value: ", config%p_refinement_threshold, ")"
                    PRINT '(A, ES9.2, A)',  "	-S_ref_th               saturation refinement threshold (value: ", config%S_refinement_threshold, ")"
#          	    elif defined(_SWE)
#    	            if defined(_ASAGI)
                        PRINT '(A, A, A)',  "	-fbath <value>          bathymetry file (value: ", trim(config%s_bathymetry_file), ")"
                        PRINT '(A, A, A)',  "	-fdispl <value>         displacement file (value: ", trim(config%s_displacement_file), ")"
#                   endif

                    PRINT '(A, A, A)',  "	-stestpoints            probe positions (example: -stestpoints 1.2334 4.0,-7.8 0.12 (value: ", trim(config%s_testpoints), ")"
                    PRINT '(A, L, A)',  "	-asciioutput               [usage of -tout required] turns on ascii output (value: ", config%l_ascii_output, ")"
                    PRINT '(A, I0, A)', "	-asciioutput_width <value> width of ascii output (value: ", config%i_ascii_width, ")"
                    PRINT '(A, ES8.1, A)',  "	-drytolerance           dry tolerance, determines up to which water height a cell is considered dry (value: ", config%dry_tolerance, " m)"
#         	    elif defined(_FLASH)
                    PRINT '(A, A, A)',  "	-fbath <value>          bathymetry file (value: ", trim(config%s_bathymetry_file), ")"
                    PRINT '(A, A, A)',  "	-fdispl <value>         displacement file (value: ", trim(config%s_displacement_file), ")"
#               endif
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
        logical                                 :: b_valid

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
#       elif defined(__INTEL_COMPILER)
            _log_write(0, '(" Compiler: Intel")')
#       else
            _log_write(0, '(" Compiler: Not detected")')
#       endif

        _log_write(0, '(" Sections per thread: ", I0)') config%i_sections_per_thread
        _log_write(0, '(" Adaptivity: min depth: ", I0, ", max depth: ", I0, ", start depth: ", I0)') config%i_min_depth, config%i_max_depth, config%i_start_depth

        _log_write(0, '(" Load balancing: timed load estimate: ", A, ", split sections: ", A, ", serial: ", A)') logical_to_char(config%l_timed_load), logical_to_char(config%l_split_sections), logical_to_char(config%l_serial_lb)
        _log_write(0, '(" Load balancing: cell weight: ", F0.2, ", boundary weight: ", F0.2)') config%r_cell_weight, config%r_boundary_weight

        _log_write(0, '(" Scenario: max time steps: ", I0, ", max time: ", ES9.2, ", output step: ", ES9.2)') config%i_max_time_steps, config%r_max_time, config%r_output_time_step

        !check if the output directory exists

#       if defined(__INTEL_COMPILER)
            !This solution is not standard but works for the Intel compiler.
            inquire(directory=trim(config%output_dir), exist=b_valid)
#       else
            !This solution works for the GNU compiler and hopefully for others, too.
            inquire(file=trim(config%output_dir) // '/.', exist=b_valid)
#       endif

        try(b_valid, "output directory " // trim(config%output_dir) // " could not be found")

        _log_write(0, '(" Scenario: output directory: ", A)') trim(config%output_dir)
        _log_write(0, '(" Scenario: courant number: ", F0.3)') config%courant_number

#       if defined (_UPWIND_FLUX)
           _log_write(0, '(" Scenario: Flux solver: ", A)') "Upwind"
#       elif defined (_LF_FLUX)
           _log_write(0, '(" Scenario: Flux solver: ", A)') "Lax Friedrichs"
#       elif defined (_LLF_FLUX)
            _log_write(0, '(" Scenario: Flux solver: ", A)')  "Local Lax Friedrichs"
#       elif defined(_LF_BATH_FLUX)
            _log_write(0, '(" Scenario: Flux solver: ", A)')  "Lax Friedrichs (Bathymetry)"
#       elif defined(_LLF_BATH_FLUX)
            _log_write(0, '(" Scenario: Flux solver: ", A)') "Local Lax Friedrichs (Bathymetry)"
#       elif defined(_FWAVE_FLUX)
            _log_write(0, '(" Scenario: Flux solver: ", A)')  "F-Wave"
#       elif defined(_AUG_RIEMANN_FLUX)
            _log_write(0, '(" Scenario: Flux solver: ", A)')  "Augmented Riemann"
#       else
#           error Invalid flux solver!
#       endif

#       if defined (_ADAPT_INTEGRATE)
            _log_write(0, '(" Scenario: Input data refinement: ", A)') "Integrate"
#       elif defined (_ADAPT_SAMPLE)
            _log_write(0, '(" Scenario: Input data refinement: ", A)') "Sample"
#       else
#           error Invalid data refinement method!
#       endif

#		if defined(_DARCY)
#		    if (_DARCY_LAYERS > 0)
                _log_write(0, '(" Darcy: 3D with ", I0, " layers")') _DARCY_LAYERS
#           else
                _log_write(0, '(" Darcy: 2D")')
#           endif

#           if defined (_PERM_MEAN_ARITHMETIC)
               _log_write(0, '(" Darcy: Permeability averaging: ", A)') "Arithmetic"
#           elif defined (_PERM_MEAN_GEOMETRIC)
               _log_write(0, '(" Darcy: Permeability averaging: ", A)') "Geometric"

#               if defined(_DEBUG)
#                   warning Geometric averaging causes intended floating point overflows, use arithmetic averaging to catch floating point exceptions
#               endif
#           elif defined (_PERM_MEAN_HARMONIC)
               _log_write(0, '(" Darcy: Permeability averaging: ", A)') "Harmonic"

#               if defined(_DEBUG)
#                   warning Harmonic averaging causes intended floating point overflows, use arithmetic averaging to catch floating point exceptions
#               endif
#           else
#               error Invalid permeability averaging!
#           endif

#           if defined (_DARCY_MOB_LINEAR)
               _log_write(0, '(" Darcy: Mobility term: ", A)') "Linear"
#           elif defined (_DARCY_MOB_QUADRATIC)
               _log_write(0, '(" Darcy: Mobility term: ", A)') "Quadratic"
#           elif defined (_DARCY_MOB_BROOKS_COREY)
               _log_write(0, '(" Darcy: Mobility term: ", A)') "Brooks-Corey"
#           else
#               error Invalid mobility term!
#           endif

            !check if the data files exixst
#    	    if defined(_ASAGI)
                inquire(file=trim(config%s_permeability_file), exist=b_valid)
                try(b_valid, "permeability data file " // trim(config%s_permeability_file) // " could not be found")

                inquire(file=trim(config%s_porosity_file), exist=b_valid)
                try(b_valid, "porosity data file " // trim(config%s_porosity_file) // " could not be found")

                _log_write(0, '(" Darcy: permeability file: ", A, ", porosity file: ", A)') trim(config%s_permeability_file),  trim(config%s_porosity_file)
#           endif

            _log_write(0, '(" Darcy: vicosities: wetting phase: ", ES8.1, ", non-wetting phase: ", ES8.1)') config%r_nu_w_SI, config%r_nu_n_SI
            _log_write(0, '(" Darcy: densities: wetting phase: ", ES8.1, ", non-wetting phase: ", ES8.1)') config%r_rho_w_SI, config%r_rho_n_SI
            _log_write(0, '(" Darcy: injection well pressure: ", ES8.1, ", production well pressure: ", ES8.1)') config%r_p_in_AU, config%r_p_prod_AU
            _log_write(0, '(" Darcy: linear solver: ", I0, ": ", A)') config%i_lsolver, trim(lsolver_to_char(config%i_lsolver))
            _log_write(0, '(" Darcy: linear solver: error bound: ", ES8.1, ", max iterations: ", I0)') config%r_epsilon, config%i_max_iterations
            _log_write(0, '(" Darcy: CG restart interval: ", I0)') config%i_CG_restart
            _log_write(0, '(" Darcy: pressure refinement threshold: ", ES9.2, ", saturation refinement threshold: ", ES9.2)') config%p_refinement_threshold, config%S_refinement_threshold
#		elif defined(_FLASH)
            _log_write(0, '(" Flash: bathymetry file: ", A, ", displacement file: ", A)') trim(config%s_bathymetry_file), trim(config%s_displacement_file)
#		elif defined(_SWE)

            !check if the data files exixst
#    	    if defined(_ASAGI)
                inquire(file=trim(config%s_bathymetry_file), exist=b_valid)
                try(b_valid, "bathymetry data file " // trim(config%s_bathymetry_file) // " could not be found")

                inquire(file=trim(config%s_displacement_file), exist=b_valid)
                try(b_valid, "displacement data file " // trim(config%s_displacement_file) // " could not be found")

                _log_write(0, '(" SWE: bathymetry file: ", A, ", displacement file: ", A)') trim(config%s_bathymetry_file), trim(config%s_displacement_file)
#           endif

            _log_write(0, '(" SWE: dry_tolerance: ", ES8.1)') config%dry_tolerance
#		endif

        _log_write(0, '("")')
    end subroutine

#   if defined(_SWE)
        subroutine parse_testpoints(config)
            class(t_config), intent(inout)          :: config

            !local variables
            logical					:: l_wrong_format, l_point, l_comma, l_space, l_number_pre, l_number_post, l_sign, l_ycoord
                integer          			:: i, i_error, points, j, coordstart(50,2), k, counter
                character(512)                          :: checkstring

            ! check for correctness
            l_wrong_format = .false.
            checkstring = config%s_testpoints

            l_sign = .true.
            l_number_pre = .true.
            l_number_post = .false.
            l_space = .false.
            l_comma = .false.
            l_point = .false.
            l_ycoord = .false.

            counter = len(trim(config%s_testpoints))

            do while ((l_wrong_format .eqv. .false.) .and. (counter > 0))
                if ((l_sign .eqv. .true.) .and. (checkstring(1:1) == "-")) then
                    ! sign ('-' only)
                    checkstring = checkstring (2:)
                    counter = counter - 1
                    l_number_pre = .true.
                    l_number_post = .false.
                    l_sign = .false.
                    l_space = .false.
                    l_comma = .false.
                    l_point = .false.
                else if ((l_number_pre .eqv. .true.) .and. &
                   (checkstring(1:1) == "1" .or. checkstring(1:1) == "2" .or. &
                   checkstring(1:1) == "3" .or. checkstring(1:1) == "4" .or. &
                   checkstring(1:1) == "5" .or. checkstring(1:1) == "6" .or. &
                   checkstring(1:1) == "7" .or. checkstring(1:1) == "8" .or. &
                   checkstring(1:1) == "9" .or. checkstring(1:1) == "0")) then
                    ! number before the point
                    checkstring = checkstring(2:)
                    counter = counter - 1
                    l_number_pre = .true.
                    l_number_post = .false.
                    l_sign = .false.
                    l_space = .false.
                    l_comma = .false.
                    l_point = .true.
                else if ((l_number_post .eqv. .true.) .and. &
                   (checkstring(1:1) == "1" .or. checkstring(1:1) == "2" .or. &
                   checkstring(1:1) == "3" .or. checkstring(1:1) == "4" .or. &
                   checkstring(1:1) == "5" .or. checkstring(1:1) == "6" .or. &
                   checkstring(1:1) == "7" .or. checkstring(1:1) == "8" .or. &
                   checkstring(1:1) == "9" .or. checkstring(1:1) == "0")) then
                    ! number behind the point
                    checkstring = checkstring(2:)
                    counter = counter - 1
                    l_number_pre = .false.
                    l_number_post = .true.
                    l_sign = .false.
                    l_space = .true.
                    l_comma = .true.
                    l_point = .false.
                else if ((l_point .eqv. .true.) .and. (checkstring(1:1) == ".")) then
                    ! point
                    checkstring = checkstring(2:)
                    counter = counter - 1
                    l_number_pre = .false.
                    l_number_post = .true.
                    l_sign = .false.
                    l_space = .false.
                    l_comma = .false.
                    l_point = .false.
                else if ((l_comma .eqv. .true.) .and. (l_ycoord .eqv. .true.) &
                    .and. (checkstring(1:1) == ",")) then
                    ! comma
                    checkstring = checkstring(2:)
                    counter = counter - 1
                    l_number_pre = .true.
                    l_number_post = .false.
                    l_sign = .true.
                    l_space = .false.
                    l_comma = .false.
                    l_point = .false.
                    l_ycoord = .false.
                else if ((l_space .eqv. .true.) .and. (l_ycoord .eqv. .false.) &
                   .and. (checkstring(1:1) == " ")) then
                    ! space
                    checkstring = checkstring(2:)
                    counter = counter - 1
                    l_number_pre = .true.
                    l_number_post = .false.
                    l_sign = .true.
                    l_space = .false.
                    l_comma = .false.
                    l_point = .false.
                    l_ycoord = .true.
                else
                    write (*,'(A,$)') 'Position with wrong symbol (counting from the end):'
                    write (*,*) counter
                    l_wrong_format = .true.
                end if
            end do

            try((.not. l_wrong_format), 'Error in submitted testpoints')

            !convert stestpoints (string) to r_testpoints (real array)
            points = 0
            !k = 1
            if (l_wrong_format .eqv. .false.) then
                points = 1
                coordstart(1,1) = 1
                do j=1, len(trim(config%s_testpoints))
                if (config%s_testpoints(j:j) == " ") then
                    coordstart(points,2) = j+1
                end if
                if (config%s_testpoints(j:j) == ",") then
                    points = points + 1
                    coordstart(points,1) = j+1
                end if
                end do

                allocate (config%r_testpoints(points,2), stat = i_error); assert_eq(i_error,0)

                 do j=1, points
                    read(unit=config%s_testpoints(coordstart(j,1):), fmt=*) config%r_testpoints(j,1)
                    read(unit=config%s_testpoints(coordstart(j,2):), fmt=*) config%r_testpoints(j,2)
                 end do
            end if
        end subroutine
#   endif
end module
