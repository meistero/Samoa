// Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
// Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
// This program is licensed under the GPL, for details see the file LICENSE

// Note: this example requires C++0x for lambda expressions

#include <stdio.h>
#include "samoa.h"
#include <assert.h>

int i_sections;
Samoa_grid* grid;

int main(int argc, char* args[]) {
    // run a few test cases with the samoa library

    samoa_get_grid(i_sections, grid);

    for (int is = 0; is < i_sections; is++) {
        printf("section index: %i: cells: %lli edges: %lli nodes: %lli\n", is, grid[is].i_cells, grid[is].i_edges, grid[is].i_nodes);
    }

	samoa_run_kernel([] (const int is, const long long ic, char& refinement) {
    	printf("  section: %i, cell: %lli edges: %lli %lli %lli nodes: %lli %lli %lli \n", is, ic,
            grid[is].cells_to_edges[3 * ic + 0], grid[is].cells_to_edges[3 * ic + 1], grid[is].cells_to_edges[3 * ic + 2],
            grid[is].cells_to_nodes[3 * ic + 0], grid[is].cells_to_nodes[3 * ic + 1], grid[is].cells_to_nodes[3 * ic + 2]);

    	assert(0 <= ic && ic < grid[is].i_cells);
    	assert(0 <= grid[is].cells_to_edges[3 * ic + 0] && grid[is].cells_to_edges[3 * ic + 0] < grid[is].i_edges);
    	assert(0 <= grid[is].cells_to_edges[3 * ic + 1] && grid[is].cells_to_edges[3 * ic + 1] < grid[is].i_edges);
    	assert(0 <= grid[is].cells_to_edges[3 * ic + 2] && grid[is].cells_to_edges[3 * ic + 2] < grid[is].i_edges);
    	assert(0 <= grid[is].cells_to_nodes[3 * ic + 0] && grid[is].cells_to_nodes[3 * ic + 0] < grid[is].i_nodes);
    	assert(0 <= grid[is].cells_to_nodes[3 * ic + 1] && grid[is].cells_to_nodes[3 * ic + 1] < grid[is].i_nodes);
    	assert(0 <= grid[is].cells_to_nodes[3 * ic + 2] && grid[is].cells_to_nodes[3 * ic + 2] < grid[is].i_nodes);
	});

    for (int i = 0; i < 2; i++) {
        //refine each cell
        samoa_run_kernel([] (const int is, const long long ic, char& refinement) {
            refinement = 1;
        });

        samoa_get_grid(i_sections, grid);

        for (int is = 0; is < i_sections; is++) {
            printf("section index: %i: cells: %lli edges: %lli nodes: %lli\n", is, grid[is].i_cells, grid[is].i_edges, grid[is].i_nodes);
        }
    }

    samoa_run_kernel([] (const int is, const long long ic, char& refinement) {
    	printf("  section: %i, cell: %lli edges: %lli %lli %lli nodes: %lli %lli %lli \n", is, ic,
            grid[is].cells_to_edges[3 * ic + 0], grid[is].cells_to_edges[3 * ic + 1], grid[is].cells_to_edges[3 * ic + 2],
            grid[is].cells_to_nodes[3 * ic + 0], grid[is].cells_to_nodes[3 * ic + 1], grid[is].cells_to_nodes[3 * ic + 2]);

    	assert(0 <= ic && ic < grid[is].i_cells);
    	assert(0 <= grid[is].cells_to_edges[3 * ic + 0] && grid[is].cells_to_edges[3 * ic + 0] < grid[is].i_edges);
    	assert(0 <= grid[is].cells_to_edges[3 * ic + 1] && grid[is].cells_to_edges[3 * ic + 1] < grid[is].i_edges);
    	assert(0 <= grid[is].cells_to_edges[3 * ic + 2] && grid[is].cells_to_edges[3 * ic + 2] < grid[is].i_edges);
    	assert(0 <= grid[is].cells_to_nodes[3 * ic + 0] && grid[is].cells_to_nodes[3 * ic + 0] < grid[is].i_nodes);
    	assert(0 <= grid[is].cells_to_nodes[3 * ic + 1] && grid[is].cells_to_nodes[3 * ic + 1] < grid[is].i_nodes);
    	assert(0 <= grid[is].cells_to_nodes[3 * ic + 2] && grid[is].cells_to_nodes[3 * ic + 2] < grid[is].i_nodes);
    });

    return 0;
}
