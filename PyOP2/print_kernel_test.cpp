// Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
// Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
// This program is licensed under the GPL, for details see the file LICENSE

// Note: this example requires C++0x for lambda expressions

#include <stdio.h>
#include "samoa.h"

int main(int argc, char* args[]) {
    // run a few test cases with the samoa library

	samoa_run_c_kernel([] (const int ic, const int ie[3], const int iv[3], const double coords[6], char* refinement) {
    	printf("cell index: %i edge indices: %i %i %i vertex indices: %i %i %i\n", ic, ie[0], ie[1], ie[2], iv[0], iv[1], iv[2]);
	});

	samoa_run_c_kernel([] (const int ic, const int ie[3], const int iv[3], const double coords[6], char* refinement) {
    	printf("cell index: %i coords: (%.3f %.3f) (%.3f %.3f) (%.3f %.3f)\n", ic, coords[0], coords[1], coords[2], coords[3], coords[4], coords[5]);
	});

    //refine each cell of the grid
	samoa_run_c_kernel([] (const int ic, const int ie[3], const int iv[3], const double coords[6], char* refinement) {
    	*refinement = 1;
	});

	samoa_run_c_kernel([] (const int ic, const int ie[3], const int iv[3], const double coords[6], char* refinement) {
    	printf("cell index: %i edge indices: %i %i %i vertex indices: %i %i %i\n", ic, ie[0], ie[1], ie[2], iv[0], iv[1], iv[2]);
	});

	samoa_run_c_kernel([] (const int ic, const int ie[3], const int iv[3], const double coords[6], char* refinement) {
    	printf("cell index: %i coords: (%.3f %.3f) (%.3f %.3f) (%.3f %.3f)\n", ic, coords[0], coords[1], coords[2], coords[3], coords[4], coords[5]);
	});

    return 0;
}