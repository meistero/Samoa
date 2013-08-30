# Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
# Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
# This program is licensed under the GPL, for details see the file LICENSE

#default compiler and compiler-specific flags

CXX				= mpicxx.mpich2 -cxx=icc
CXXFLAGS		= -O3 -std=c++0x
LDFLAGS 		= "-Wl,-rpath,../bin" -L"../bin/"

#----------------------------------------------------------------#
# Build targets                                                  #
#----------------------------------------------------------------#

#default target

all:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o print_kernel_test -lsamoa_pyop2_noasagi print_kernel_test.cpp

clean:
	@rm -f print_kernel_test
