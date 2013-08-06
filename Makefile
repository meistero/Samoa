# Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
# Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
# This program is licensed under the GPL, for details see the file LICENSE


#----------------------------------------------------------------#
# Compiler Flags                                                 #
#----------------------------------------------------------------#

# Usage: call "make <scenario> [<FLAG=VALUE>]*" or "make [<FLAG=VALUE>]*"
#
# make flags:
#  SCENARIO=DARCY|HEAT_EQ|SWE|TESTS
#  SWE_SOLVER=LAX_FRIEDRICHS|LAX_FRIEDRICHS_BATH|FWAVE|SSQ_FWAVE|AUG_JCP
#  TARGET=DEBUG|PROF|OPT
#  MPI=YES|NO
#  OPENMP=YES|NO
#  STD=YES|NO
#  ASAGI=YES|NO
#  ASAGI_TIMING=YES|NO
#  DEBUG_LEVEL = (0-7)
#  ASSERT = YES|NO

#default compiler and compiler-specific flags

FFLAGS			= -implicitnone -nologo -fpp -I"./" -I"Samoa/"
EXEC 			= bin/samoa

#default values for compilation switches

SCENARIO		?= DARCY
TARGET			?= OPT
OPENMP			?= YES
MPI 			?= YES
STD 			?= NO
VEC_REPORT		?= 0
SWE_SOLVER		?= AUG_RIEMANN
ASAGI_TIMING 	?= NO

#check switches, set flags of dependent switches and compiler flags accordingly

ifeq ($(MPI), YES)
  FC			= mpif90 -f90=ifort
  LOADER		= mpif90 -f90=ifort
else ifeq ($(MPI), NO)
  FC			= ifort
  LOADER		= ifort
else
  $(error Invalid value for MPI: $(MPI))
endif

ifeq ($(SCENARIO), DARCY)
  EXEC 			:= $(EXEC)_darcy
  FFLAGS		+= -D_DARCY
  ASAGI 		?= YES
else ifeq ($(SCENARIO), SWE)
  EXEC			:= $(EXEC)_swe
  FFLAGS		+= -D_SWE
  ASAGI			?= YES
else ifeq ($(SCENARIO), NUMA)
  EXEC			:= $(EXEC)_numa
  FFLAGS		+= -D_NUMA
  ASAGI			?= NO
else ifeq ($(SCENARIO), HEAT_EQ)
  EXEC			:= $(EXEC)_heq
  FFLAGS		+= -D_HEAT_EQ
  ASAGI			?= NO
else ifeq ($(SCENARIO), TESTS)
  EXEC			:= $(EXEC)_tests
  FFLAGS		+= -D_TESTS
  ASAGI			?= NO
else
  $(error Invalid value for SCENARIO: $(SCENARIO))
endif

ifeq ($(OPENMP), YES)
  FFLAGS		+= -openmp -D_OMP
  LDFLAGS		+= -openmp
else ifeq ($(OPENMP), NO)
  EXEC			:= $(EXEC)_noomp
  FFLAGS		+= -openmp-stubs
  LDFLAGS		+= -openmp-stubs
else
  $(error Invalid value for OPENMP: $(OPENMP))
endif

ifeq ($(MPI), YES)
  FFLAGS		+= -D_MPI
else ifeq ($(MPI), NO)
  EXEC 			:= $(EXEC)_nompi
else
  $(error Invalid value for MPI: $(MPI))
endif

ifeq ($(ASAGI), YES)
  FFLAGS 		+= -D_ASAGI -I"ASAGI/include"
  LDFLAGS 		+= "-Wl,-rpath,ASAGI" -L"ASAGI/"

  ifeq ($(OPENMP), YES)
    LDFLAGS		+= -lasagi
  else
    LDFLAGS		+= -lasagi_nomt
  endif

  ifeq ($(ASAGI_TIMING), YES)
    FFLAGS 		+= -D_ASAGI_TIMING
  else ifeq ($(ASAGI_TIMING), NO)
  else
    $(error Invalid value for ASAGI_TIMING: $(ASAGI_TIMING))
  endif
else ifeq ($(ASAGI), NO)
  EXEC			:= $(EXEC)_noasagi
else
  $(error Invalid value for ASAGI: $(ASAGI))
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
  FFLAGS 		+= -fast -inline-level=2 -funroll-loops -unroll
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

ifeq ($(STD), YES)
  FFLAGS 		+= -std
else ifeq ($(STD), NO)
else
  $(error Invalid value for STD: $(STD))
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
SWE/SWE_initialize.f90 \
SWE/SWE_output.f90 \
SWE/SWE_xml_output.f90 \
SWE/SWE_euler_timestep.f90 \
SWE/SWE_adapt.f90 \
Samoa/Samoa.f90 \
Samoa/Tools_quadrature_rule_base.f90 \
SFC_grid.f90 \
SFC_fine_grid.f90 \
SFC_node_traversal.f90 \
SFC_edge_traversal.f90 \
SFC_data_types.f90 \
LIB_VTK_IO.f90\
M_kracken.f90\
Tools_noise.f90 \
Tools_log.f90 \
Conformity/Conformity.f90 \
NUMA/NUMA.f90 \
NUMA/NUMA_local_function_spaces.f90 \
NUMA/NUMA_data_types.f90 \
NUMA/NUMA_basis.f90 \
NUMA/NUMA_initialize.f90 \
NUMA/NUMA_output.f90 \
NUMA/NUMA_xml_output.f90 \
NUMA/NUMA_euler_timestep.f90 \
NUMA/NUMA_adapt.f90 \
NUMA/NUMA_constants.f90 \
geoclaw/c_bind_riemannsolvers.f90 \

F77_SOURCES = geoclaw/riemannsolvers.f

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

dirs:
	@mkdir -p bin output

compile: $(EXEC)
	@rm -f $(F90_OBJS) $(F77_OBJS) *.mod

$(EXEC): $(F90_OBJS) $(F77_OBJS) dirs
	$(LOADER) -o $@ $(F90_OBJS) $(F77_OBJS) $(LDFLAGS)

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

%.o: %.f
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	@rm -f bin/* $(F90_OBJS) $(F77_OBJS) *.mod

-include dependency.mk
