// Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
// Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
// This program is licensed under the GPL, for details see the file LICENSE

/** \brief Returns the grid data in a specified section
 *
 * \param section_index     (in)index of the current section
 * \param cells             (out)number of cells
 * \param edges             (out)number of edges
 * \param nodes             (out)number of nodes
 * \param cells_to_edges    (out)map from cells to edges
 * \param cells_to_nodes    (out)map from cells to nodes
 * \param edges_to_nodes    (out)map from edges to nodes
 * \param coords            (out)node coordinates
 *
 */
extern "C" void samoa_get_grid(const int section_index,
                               long long& i_cells, long long& i_edges, long long& i_nodes,
                               long long*& cells_to_edges, long long*& cells_to_nodes, long long*& edges_to_nodes,
                               double*& coords);


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
extern "C" double* samoa_malloc(const int section_index, const long long dofs, const int dim, const long long cells_to_dofs[], const long long edges_to_dofs[], const long long nodes_to_dofs[]);

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
 * 4  5
 * | 7 \
 * |    \
 * 1--6--0
 *
 *
 * \param kernel            (in)kernel function
 * \param section_index     (in)index of the current section
 * \param cell_index        (in)index of the cell
 * \param refinement        (out)refinement flag: set to 1 for cell refinement, -1 for coarsening, 0 for no changes (default: 0)
 *
 */
extern "C" void samoa_run_kernel(void (*kernel)(const int section_index, const long long cell_index, char& refinement));
