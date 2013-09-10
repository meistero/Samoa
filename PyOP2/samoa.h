// Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
// Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
// This program is licensed under the GPL, for details see the file LICENSE

/** \brief Describes all grid entities and their relations via maps
 * Multi-dimensional data is stored in row-wise layout, indices are zero-based.
 *
 * \param cells             number of cells
 * \param edges             number of edges
 * \param nodes             number of nodes
 * \param cells_to_edges    map from cells to edges
 * \param cells_to_nodes    map from cells to nodes
 * \param edges_to_nodes    map from edges to nodes
 * \param coords            node coordinates
 *
 */
struct Samoa_grid {
    long long i_cells, i_edges, i_nodes;
    long long* cells_to_edges, *cells_to_nodes, *edges_to_nodes;
    double* coords;
};

/** \brief Returns the grid data of a specified section
 *
 * \param i_sections        (out)number of sections in the grid
 * \param grid              (out)array of samoa sections that describe all grid entities and their relations
 *
 */
extern "C" void samoa_get_grid(int& i_sections, Samoa_grid*& sections);


/** \brief Allocates a samoa data array in a specified section
 *
 * \param section_index     (in)index of the current section
 * \param dofs              (in)number of DoFs
 * \param cells_to_dofs     (in)map from cell interiors to DoFs
 * \param edges_to_dofs     (in)map from edge interiors to DoFs
 * \param nodes_to_dofs     (in)map from nodes to DoFs
 * \param dim               (in)DoF dimension
 * \return                  (out)array pointer
 *
 */
extern "C" double* samoa_malloc(const int section_index, const long long dofs, const int dim, const long long* cells_to_dofs, const long long* edges_to_dofs, const long long* nodes_to_dofs);

/** \brief Deallocates a samoa data array
 *
 * \param data            (in)array pointer
 *
 */
extern "C" void samoa_free(double* data);


/** \brief Executes a kernel on a triangular mesh in Sierpinski order
 * The following order is assumed for local DoF numbering:
 *
 * 2
 * |\
 * | \
 * 3  4
 * | 6 \
 * |    \
 * 1--5--0
 *
 *
 * \param kernel            (in)kernel function
 * \param section_index     (in)index of the current section
 * \param cell_index        (in)index of the cell
 * \param refinement        (out)refinement flag: set to 1 for cell refinement, -1 for coarsening, 0 for no changes (default: 0)
 *
 */
extern "C" void samoa_run_kernel(void (*kernel)(const int section_index, const long long cell_index, char& refinement));
