# Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
# Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
# This program is licensed under the GPL, for details see the file LICENSE


#----------------------------------------------------------------#
# Compiler Flags                                                 #
#----------------------------------------------------------------#

# Usage: call "make [<scenario>] [<FLAG=VALUE>]*"
#
# make flags:
#  SCENARIO=DARCY|HEAT_EQ|SWE|TESTS|GENERIC
#  SWE_SOLVER=LAX_FRIEDRICHS|LAX_FRIEDRICHS_BATH|FWAVE|SSQ_FWAVE|AUG_RIEMANN
#  TARGET=DEBUG|PROF|OPT
#  MPI=DEFAULT|MPICH2|OPENMPI|INTEL|NO
#  OPENMP=YES|TASKS|NO
#  STD_FORTRAN=YES|NO
#  ASAGI=STANDARD|NUMA|NO
#  ASAGI_TIMING=YES|NO
#  DEBUG_LEVEL = (0-7)
#  ASSERT = YES|NO
#  VEC_REPORT = (0-3)
#  ASAGI_DIR = <path>

#default compiler and compiler-specific flags

FFLAGS			= -implicitnone -nologo -fpp -I"./" -I"Samoa/"
EXEC 			= samoa

#default values for compilation switches

SCENARIO		?= DARCY
SWE_SOLVER		?= AUG_RIEMANN
TARGET			?= OPT
MPI 			?= DEFAULT
OPENMP			?= TASKS
STD_FORTRAN		?= NO
ASAGI_TIMING	?= NO
VEC_REPORT		?= 0
ASAGI_DIR		?= "./ASAGI"

#check switches, set flags of dependent switches and compiler flags accordingly

ifeq ($(MPI), DEFAULT)
  FC			= MPICH_F90=ifort OMPI_FC=ifort I_MPI_F90=ifort mpif90
  LOADER		= MPICH_F90=ifort OMPI_FC=ifort I_MPI_F90=ifort mpif90
else ifeq ($(MPI), OPENMPI)
  FC			= OMPI_FC=ifort mpif90.openmpi
  LOADER		= OMPI_FC=ifort mpif90.openmpi
else ifeq ($(MPI), MPICH2)
  FC			= MPICH_F90=ifort mpif90.mpich2
  LOADER		= MPICH_F90=ifort mpif90.mpich2
else ifeq ($(MPI), INTEL)
  FC			= I_MPI_F90=ifort mpif90.intel
  LOADER		= I_MPI_F90=ifort mpif90.intel
else ifeq ($(MPI), NO)
  FC			= ifort
  LOADER		= ifort
else
  $(error Invalid value for MPI: $(MPI))
endif

ifeq ($(SCENARIO), DARCY)
  EXEC 			:= $(EXEC)_darcy
  FFLAGS		+= -D_DARCY
  ASAGI 		?= STANDARD
  LIB 			?= NO
else ifeq ($(SCENARIO), GENERIC)
  EXEC 			:= $(EXEC)_generic
  FFLAGS		+= -D_GENERIC
  ASAGI 		?= NO
  LIB 			?= YES
else ifeq ($(SCENARIO), SWE)
  EXEC			:= $(EXEC)_swe
  FFLAGS		+= -D_SWE
  ASAGI			?= STANDARD
  LIB 			?= NO
else ifeq ($(SCENARIO), HEAT_EQ)
  EXEC			:= $(EXEC)_heq
  FFLAGS		+= -D_HEAT_EQ
  ASAGI			?= NO
  LIB 			?= NO
else ifeq ($(SCENARIO), TESTS)
  EXEC			:= $(EXEC)_tests
  FFLAGS		+= -D_TESTS
  ASAGI			?= NO
  LIB 			?= NO
else
  $(error Invalid value for SCENARIO: $(SCENARIO))
endif

ifeq ($(OPENMP), YES)
  EXEC			:= $(EXEC)_notasks
  FFLAGS		+= -openmp
  LDFLAGS		+= -openmp
else ifeq ($(OPENMP), TASKS)
  FFLAGS		+= -openmp -D_OPENMP_TASKS
  LDFLAGS		+= -openmp
else ifeq ($(OPENMP), NO)
  EXEC			:= $(EXEC)_noomp
  FFLAGS		+= -openmp-stubs
  LDFLAGS		+= -openmp-stubs
else
  $(error Invalid value for OPENMP: $(OPENMP))
endif

ifeq ($(MPI), NO)
  EXEC 			:= $(EXEC)_nompi
else
  FFLAGS		+= -D_MPI
endif

ifeq ($(ASAGI), STANDARD)
  FFLAGS 		+= -D_ASAGI -I$(ASAGI_DIR)"/include"
  LDFLAGS 		+= "-Wl,-rpath,"$(ASAGI_DIR) -L$(ASAGI_DIR)

  ifeq ($(OPENMP), NO)
    LDFLAGS		+= -lasagi_nomt
  else
    LDFLAGS		+= -lasagi
  endif
else ifeq ($(ASAGI), NUMA)
  FFLAGS 		+= -D_ASAGI -D_ASAGI_NUMA -I$(ASAGI_DIR)"/include"
  LDFLAGS 		+= "-Wl,-rpath,"$(ASAGI_DIR) -L$(ASAGI_DIR)
  LDFLAGS		+= -lasagi

  ifeq ($(OPENMP), NO)
    $(error ASAGI must not be NUMA if OPENMP is NO)
  endif
else ifeq ($(ASAGI), NO)
  EXEC			:= $(EXEC)_noasagi
else
  $(error Invalid value for ASAGI: $(ASAGI))
endif

ifeq ($(ASAGI_TIMING), YES)
  FFLAGS 		+= -D_ASAGI_TIMING

  ifeq ($(ASAGI), NO)
    $(error ASAGI_TIMING must not be YES if ASAGI is NO)
  endif
else ifeq ($(ASAGI_TIMING), NO)
  #nothing to do
else
  $(error Invalid value for ASAGI_TIMING: $(ASAGI_TIMING))
endif

ifeq ($(SWE_SOLVER), LF)
  FFLAGS 		+= -D_SWE_LF
  EXEC			:= $(EXEC)_lf
else ifeq ($(SWE_SOLVER), LF_BATH)
  FFLAGS 		+= -D_SWE_LF_BATH
  EXEC			:= $(EXEC)_lfbath
else ifeq ($(SWE_SOLVER), LLF)
  FFLAGS 		+= -D_SWE_LLF
  EXEC			:= $(EXEC)_llf
else ifeq ($(SWE_SOLVER), LLF_BATH)
  FFLAGS 		+= -D_SWE_LLF_BATH
  EXEC			:= $(EXEC)_llfbath
else ifeq ($(SWE_SOLVER), FWAVE)
  FFLAGS 		+= -D_SWE_FWAVE
  EXEC			:= $(EXEC)_fwave
else ifeq ($(SWE_SOLVER), AUG_RIEMANN)
  FFLAGS 		+= -D_SWE_AUG_RIEMANN
else
  $(error Invalid value for SWE_SOLVER: $(SWE_SOLVER))
endif

ifeq ($(TARGET), DEBUG)
  EXEC			:= $(EXEC)_debug
  DEBUG_LEVEL	?= 3
  ASSERT 		?= YES
  FFLAGS 		+= -g -O0 -traceback -check all -debug all -fpe0
  LDFLAGS 		+= -g -O0 -traceback -check all -debug all -fpe0
else ifeq ($(TARGET), PROF)
  EXEC 			:= $(EXEC)_prof
  DEBUG_LEVEL 	?= 1
  ASSERT 		?= NO
  FFLAGS 		+= -g -trace -fast -inline-level=0 -funroll-loops -unroll
  LDFLAGS 		+= -g -trace -O3 -ip -ipo
else ifeq ($(TARGET), OPT)
  DEBUG_LEVEL 	?= 1
  ASSERT 		?= NO
  FFLAGS 		+= -fast -align all -inline-level=2 -no-inline-min-size -no-inline-max-size -no-inline-max-total-size -no-inline-max-per-routine -no-inline-max-per-compile -no-inline-factor -funroll-loops -unroll
  LDFLAGS 		+= -O3 -ip -ipo
else
  $(error Invalid value for TARGET: $(TARGET))
endif

ifdef VEC_REPORT
  LDFLAGS 		+= -vec-report$(VEC_REPORT)
endif

ifdef DEBUG_LEVEL
  FFLAGS 		+= -D_DEBUG_LEVEL=$(DEBUG_LEVEL)
endif

ifeq ($(ASSERT), YES)
  FFLAGS 		+= -D_ASSERT
else ifeq ($(ASSERT), NO)

else
  $(error Invalid value for ASSERT: $(ASSERT))
endif

ifeq ($(STD_FORTRAN), YES)
  FFLAGS 		+= -std
else ifeq ($(STD_FORTRAN), NO)
else
  $(error Invalid value for STD_FORTRAN: $(STD_FORTRAN))
endif

ifeq ($(LIB), YES)
  EXEC			:= bin/lib$(EXEC).so
  FFLAGS		+= -fpic
  LDFLAGS		+= -fpic -shared
else ifeq ($(LIB), NO)
  EXEC			:= bin/$(EXEC)
else
  $(error Invalid value for LIB: $(LIB))
endif

#----------------------------------------------------------------#
# BUILD RULES                                                    #
#----------------------------------------------------------------#

F90_SOURCES = \
SFC_main.f90 \
SFC_traversal.f90 \
Tests/Tests.f90 \
Tests/Tests_data_types.f90 \
Tests/Tests_initialize.f90 \
Tests/Tests_node_dummy_traversal.f90 \
Tests/Tests_consistency_traversal.f90 \
Tests/Tests_flops_traversal.f90 \
Tests/Tests_memory_traversal.f90 \
Tests/Tests_basis_functions.f90 \
Generic/Generic.f90 \
Generic/Generic_data_types.f90 \
Generic/Generic_initialize.f90 \
Generic/Generic_template.f90 \
Generic/Generic_adapt_template.f90 \
Darcy/Darcy.f90 \
Darcy/Darcy_local_function_spaces.f90 \
Darcy/Darcy_data_types.f90 \
Darcy/Darcy_basis.f90 \
Darcy/Darcy_node_dummy.f90 \
Darcy/Darcy_initialize.f90 \
Darcy/Darcy_output.f90 \
Darcy/Darcy_xml_output.f90 \
Darcy/Darcy_laplace_jacobi.f90 \
Darcy/Darcy_laplace_cg.f90 \
Darcy/Darcy_grad_p.f90 \
Darcy/Darcy_transport_eq.f90 \
Darcy/Darcy_permeability.f90 \
Darcy/Darcy_adapt.f90 \
Heat_Equation/Heat_Eq.f90 \
Heat_Equation/Heat_Eq_local_function_spaces.f90 \
Heat_Equation/Heat_Eq_data_types.f90 \
Heat_Equation/Heat_Eq_basis.f90 \
Heat_Equation/Heat_Eq_initialize.f90 \
Heat_Equation/Heat_Eq_output.f90 \
Heat_Equation/Heat_Eq_xml_output.f90 \
Heat_Equation/Heat_Eq_euler_timestep.f90 \
Heat_Equation/Heat_Eq_midpoint_timestep.f90 \
Heat_Equation/Heat_Eq_heun_timestep.f90 \
Heat_Equation/Heat_Eq_adapt.f90 \
SWE/SWE.f90 \
SWE/SWE_local_function_spaces.f90 \
SWE/SWE_data_types.f90 \
SWE/SWE_basis.f90 \
SWE/SWE_displace.f90 \
SWE/SWE_initialize.f90 \
SWE/SWE_output.f90 \
SWE/SWE_xml_output.f90 \
SWE/SWE_euler_timestep.f90 \
SWE/SWE_adapt.f90 \
geoclaw/c_bind_riemannsolvers.f90 \
Samoa/Samoa.f90 \
Samoa/Tools_quadrature_rule_base.f90 \
Solver/LinearSolver.f90 \
SFC_node_traversal.f90 \
SFC_edge_traversal.f90 \
SFC_data_types.f90 \
LIB_VTK_IO.f90\
M_kracken.f90\
Tools_noise.f90 \
Tools_log.f90 \
Conformity/Conformity.f90 \

F77_SOURCES = \
geoclaw/riemannsolvers.f

F90_OBJS = $(F90_SOURCES:.f90=.o)
F77_OBJS = $(F77_SOURCES:.f=.o)

#----------------------------------------------------------------#
# Build targets                                                  #
#----------------------------------------------------------------#

#default target

all: compile

#if a scenario was defined as target, recursively call make with the desired scenario as parameter

darcy:
	@$(MAKE) SCENARIO=DARCY

swe:
	@$(MAKE) SCENARIO=SWE

numa:
	@$(MAKE) SCENARIO=NUMA

heat_eq:
	@$(MAKE) SCENARIO=HEAT_EQ

tests:
	@$(MAKE) SCENARIO=TESTS

generic:
	@$(MAKE) SCENARIO=GENERIC

dirs:
	@mkdir -p bin output

compile: $(EXEC)
	@rm -f $(F90_OBJS) $(F77_OBJS) *.mod

$(EXEC): $(F90_OBJS) $(F77_OBJS) dirs
	@$(LOADER) $(LDFLAGS) -o $@ $(F90_OBJS) $(F77_OBJS)

%.o: %.f90
	@$(FC) $(FFLAGS) -c -o $@ $<

%.o: %.f
	@$(FC) $(FFLAGS) -c -o $@ $<

clean:
	@rm -f bin/* $(F90_OBJS) $(F77_OBJS) *.mod

-include dependency.mk
