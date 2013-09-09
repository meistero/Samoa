// Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
// Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
// This program is licensed under the GPL, for details see the file LICENSE

// Note: this example requires C++0x for lambda expressions

#include <stdio.h>
#include "samoa.h"

long long* cells_to_edges, *cells_to_nodes, *edges_to_nodes;
double* coords;

int main(int argc, char* args[]) {
    // run a few test cases with the samoa library
    long long i_cells, i_edges, i_nodes;

    samoa_get_grid(1, i_cells, i_edges, i_nodes, cells_to_edges, cells_to_nodes, edges_to_nodes, coords);
    printf("Grid: cells: %lli edges: %lli vertices: %lli\n", i_cells, i_edges, i_nodes);

	samoa_run_kernel([] (const int is, const long long ic, char& refinement) {
    	printf(" section index: %i, cell index: %lli\n", is, ic);
    	printf("  edge indices: %lli %lli %lli\n", cells_to_edges[3 * ic + 0], cells_to_edges[3 * ic + 1], cells_to_edges[3 * ic + 2]);
    	printf("  node indices: %lli %lli %lli\n", cells_to_nodes[3 * ic + 0], cells_to_nodes[3 * ic + 1], cells_to_nodes[3 * ic + 2]);
	});

    for (int i = 0; i < 2; i++) {
        //refine each cell
        samoa_run_kernel([] (const int is, const long long ic, char& refinement) {
            refinement = 1;
        });

        samoa_get_grid(1, i_cells, i_edges, i_nodes, cells_to_edges, cells_to_nodes, edges_to_nodes, coords);
        printf("Grid: cells: %lli edges: %lli vertices: %lli\n", i_cells, i_edges, i_nodes);
    }

    samoa_run_kernel([] (const int is, const long long ic, char& refinement) {
        printf(" section index: %i, cell index: %lli\n", is, ic);
    	printf("  edge indices: %lli %lli %lli\n", cells_to_edges[3 * ic + 0], cells_to_edges[3 * ic + 1], cells_to_edges[3 * ic + 2]);
    	printf("  node indices: %lli %lli %lli\n", cells_to_nodes[3 * ic + 0], cells_to_nodes[3 * ic + 1], cells_to_nodes[3 * ic + 2]);
    });

    return 0;
}
