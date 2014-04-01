// Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
// Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
// This program is licensed under the GPL, for details see the file LICENSE

// Note: this example requires C++0x for lambda expressions

#include <stdio.h>
#include "samoa.h"
#include <assert.h>

int main(int argc, char* args[]) {
    // run a few test cases with the samoa library
    Samoa_grid grid(3);

    for (int is = 0; is < grid.get_nsections(); is++) {
        printf("section index: %i: cells: %lli edges: %lli nodes: %lli\n", is, grid[is].ncells, grid[is].nedges, grid[is].nnodes);
    }

	grid.run_kernel(&grid, [] (const int is, const long long ic, char& refinement, void* data) {
        Samoa_grid& grid = *reinterpret_cast<Samoa_grid*>(data);

    	printf("  section: %i, cell: %lli edges: %lli %lli %lli nodes: %lli %lli %lli \n", is, ic,
            grid[is].cells_to_edges[3 * ic + 0], grid[is].cells_to_edges[3 * ic + 1], grid[is].cells_to_edges[3 * ic + 2],
            grid[is].cells_to_nodes[3 * ic + 0], grid[is].cells_to_nodes[3 * ic + 1], grid[is].cells_to_nodes[3 * ic + 2]);

    	assert(0 <= ic && ic < grid[is].ncells);
    	assert(0 <= grid[is].cells_to_edges[3 * ic + 0] && grid[is].cells_to_edges[3 * ic + 0] < grid[is].nedges);
    	assert(0 <= grid[is].cells_to_edges[3 * ic + 1] && grid[is].cells_to_edges[3 * ic + 1] < grid[is].nedges);
    	assert(0 <= grid[is].cells_to_edges[3 * ic + 2] && grid[is].cells_to_edges[3 * ic + 2] < grid[is].nedges);
    	assert(0 <= grid[is].cells_to_nodes[3 * ic + 0] && grid[is].cells_to_nodes[3 * ic + 0] < grid[is].nnodes);
    	assert(0 <= grid[is].cells_to_nodes[3 * ic + 1] && grid[is].cells_to_nodes[3 * ic + 1] < grid[is].nnodes);
    	assert(0 <= grid[is].cells_to_nodes[3 * ic + 2] && grid[is].cells_to_nodes[3 * ic + 2] < grid[is].nnodes);
	});

    for (int i = 0; i < 3; i++) {
        //refine each cell
        grid.run_kernel(&grid, [] (const int is, const long long ic, char& refinement, void* data) {
            refinement = 1;
        });

        grid.update();

        for (int is = 0; is < grid.get_nsections(); is++) {
            printf("section index: %i: cells: %lli edges: %lli nodes: %lli\n", is, grid[is].ncells, grid[is].nedges, grid[is].nnodes);
        }
    }

	grid.run_kernel(&grid, [] (const int is, const long long ic, char& refinement, void* data) {
        Samoa_grid& grid = *reinterpret_cast<Samoa_grid*>(data);

    	printf("  section: %i, cell: %lli edges: %lli %lli %lli nodes: %lli %lli %lli \n", is, ic,
            grid[is].cells_to_edges[3 * ic + 0], grid[is].cells_to_edges[3 * ic + 1], grid[is].cells_to_edges[3 * ic + 2],
            grid[is].cells_to_nodes[3 * ic + 0], grid[is].cells_to_nodes[3 * ic + 1], grid[is].cells_to_nodes[3 * ic + 2]);

    	assert(0 <= ic && ic < grid[is].ncells);
    	assert(0 <= grid[is].cells_to_edges[3 * ic + 0] && grid[is].cells_to_edges[3 * ic + 0] < grid[is].nedges);
    	assert(0 <= grid[is].cells_to_edges[3 * ic + 1] && grid[is].cells_to_edges[3 * ic + 1] < grid[is].nedges);
    	assert(0 <= grid[is].cells_to_edges[3 * ic + 2] && grid[is].cells_to_edges[3 * ic + 2] < grid[is].nedges);
    	assert(0 <= grid[is].cells_to_nodes[3 * ic + 0] && grid[is].cells_to_nodes[3 * ic + 0] < grid[is].nnodes);
    	assert(0 <= grid[is].cells_to_nodes[3 * ic + 1] && grid[is].cells_to_nodes[3 * ic + 1] < grid[is].nnodes);
    	assert(0 <= grid[is].cells_to_nodes[3 * ic + 2] && grid[is].cells_to_nodes[3 * ic + 2] < grid[is].nnodes);
	});

    return 0;
}
