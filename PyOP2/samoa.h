// Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
// Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
// This program is licensed under the GPL, for details see the file LICENSE

//> Kernel execution wrapper interface
extern "C" void samoa_run_c_kernel(void (*)(const int cell_index, const int edge_indices[3], const int vertex_indices[3], const double coords[6], char* refinement));
