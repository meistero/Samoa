# Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
# Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
# This program is licensed under the GPL, for details see the file LICENSE

#default compiler and compiler-specific flags

CXX				= mpicxx.mpich2 -cxx=g++
CXXFLAGS		= -std=c++0x
LDFLAGS 		= -L"../bin/" "-Wl,-rpath,../bin"

#----------------------------------------------------------------#
# Build targets                                                  #
#----------------------------------------------------------------#

#default target

all: CXXFLAGS += -O3
all: LDFLAGS += -lsamoa_generic_noasagi

debug: CXXFLAGS += -g
debug: LDFLAGS += -lsamoa_generic_noasagi_debug

all:
	$(CXX) -o print_kernel_test print_kernel_test.cpp $(CXXFLAGS) $(LDFLAGS)

debug:
	$(CXX) -o print_kernel_test print_kernel_test.cpp $(CXXFLAGS) $(LDFLAGS)

clean:
	@rm -f print_kernel_test
