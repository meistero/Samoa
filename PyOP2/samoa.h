// Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
// Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
// This program is licensed under the GPL, for details see the file LICENSE

/** \brief Executes a kernel on a triangular mesh in Sierpinski order
 * The following order is chosen for local edge and vertex numbering:
 * 2
 * |\
 * | \
 * 0  1
 * |   \
 * |    \
 * 1--2--0
 *
 *
 * \param cell_index        index of the cell
 * \param edge_indices      map from local to global edge indices
 * \param vertex_indices    map from local to global vbertex indices
 * \param coords            vertex coordinates in local order
 *
 */
extern "C" void samoa_run_c_kernel(void (*)(const int cell_index, const int edge_indices[3], const int vertex_indices[3], const double coords[3][2], char* refinement));
