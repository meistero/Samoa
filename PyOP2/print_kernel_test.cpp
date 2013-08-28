#include <stdio.h>
#include "samoa.h"

//> A simple test kernel
void print_indices_kernel(int* cell_index, int* edge_indices, int* vertex_indices, double* coords) {
    printf("cell index: %i edge indices: %i %i %i vertex indices: %i %i %i\n", *cell_index, edge_indices[0], edge_indices[1], edge_indices[2], vertex_indices[0], vertex_indices[1], vertex_indices[2]);
    printf("coords: (%.3f %.3f) (%.3f %.3f) (%.3f %.3f)\n", coords[0], coords[1], coords[2], coords[3], coords[4], coords[5]);
}

int main(int argc, char* args[]) {
    // run a few test cases with the samoa library

    run_c_kernel(&print_indices_kernel);

    return 0;
}
