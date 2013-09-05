// Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
// Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
// This program is licensed under the GPL, for details see the file LICENSE

/** \brief Allocates a samoa data array
 *
 * \param cell_size         (in)section index
 * \param cell_size         (in)number of DoFs on each cell
 * \param edge_size         (in)number of DoFs on each edge
 * \param vertex_size       (in)number of DoFs on each vertex
 * \param dim               (in)DoF dimension (flattened)
 * \return                  (out)handle to the array
 *
 */
extern "C" int samoa_malloc(const int cell_size, const int edge_size, const int vertex_size, const int dim);

/** \brief Deallocates a samoa data array
 *
 * \param handle            (in)handle to the array
 * \param cell_size         (in)section index
 *
 */
extern "C" void samoa_free(const int handle);


/** \brief Accesses a samoa data array on a section
 *
 * \param handle            (in)handle to the array
 * \param cell_size         (in)section index
 * \return                  (out)pointer to the array for the specified section
 *
 */
extern "C" double* samoa_access(const int handle, const int section_index);


/** \brief Executes a kernel on a triangular mesh in Sierpinski order
 * The following order is chosen for local edge and vertex numbering:
 *
 * 2
 * |\
 * | \
 * 0  1
 * |   \
 * |    \
 * 1--2--0
 *
 *
 * \param kernel            (in)kernel function
 * \param section_index     (in)index of the current section
 * \param cell_index        (in)index of the cell
 * \param edge_indices      (in)map from local to global edge indices
 * \param vertex_indices    (in)map from local to global vertex indices
 * \param coords            (in)vertex coordinates in local order
 * \param refinement        (out)refinement flag: set to 1 for cell refinement, -1 for coarsening, 0 for no changes (default: 0)
 *
 */
extern "C" void samoa_run_c_kernel(void (*kernel)(const int section_index, const int cell_index, const int edge_indices[3], const int vertex_indices[3], const double coords[3][2], char* refinement));
