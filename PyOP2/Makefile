# Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
# Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
# This program is licensed under the GPL, for details see the file LICENSE

#default compiler and compiler-specific flags

CXX				= mpicxx.mpich2 -cxx=g++
CXXFLAGS		= -O3 -std=c++0x
LDFLAGS 		= -L"../bin/" -lsamoa_pyop2_noasagi "-Wl,-rpath,../bin"

#----------------------------------------------------------------#
# Build targets                                                  #
#----------------------------------------------------------------#

#default target

all:
	$(CXX) -o print_kernel_test print_kernel_test.cpp $(CXXFLAGS) $(LDFLAGS)

clean:
	@rm -f print_kernel_test
