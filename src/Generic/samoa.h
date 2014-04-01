// Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
// Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
// This program is licensed under the GPL, for details see the file LICENSE

#include <assert.h>

/** \brief Describes all grid entities and their relations via maps.
 * Multi-dimensional data is stored in row-wise layout, indices are zero-based.
 *
 * \param ncells            number of cells
 * \param nedges            number of edges
 * \param nnodes            number of nodes
 * \param nbcells           number of boundary cells
 * \param nbedges           number of boundary edges
 * \param nbnodes           number of boundary nodes
 * \param cells_to_edges    map from cells to edges
 * \param cells_to_nodes    map from cells to nodes
 * \param edges_to_nodes    map from edges to nodes
 * \param bcells            map from boundary cells to cells
 * \param bedges            map from boundary edges to edges
 * \param bnodes            map from boundary nodes to nodes
 * \param coords            node coordinates
 *
 */
struct Samoa_section {
    long long ncells, nedges, nnodes;
    long long nbcells, nbedges, nbnodes;
    long long* cells_to_edges, *cells_to_nodes, *edges_to_nodes;
    long long* bcells, *bedges, *bnodes;
    double* coords;
};

/** \brief Creates a samoa grid with an initial depth and a section partitioning per thread.
 * The grid is not balanced and may be empty for the current process. It can be retrieved using `samoa_get_grid`
 *
 * \param i_sections_per_thread     (in)number of sections per thread (must not be < 0, usually a good value would be around 16)
 * \param i_depth                   (in)initial depth
 * \return                          handle to the grid
 *
 */
extern "C" int samoa_create_grid(const int i_sections_per_thread, const char i_depth);

/** \brief Destroys the samoa grid.
 *
 * \param handle            (in)handle to the grid
 */
extern "C" void samoa_destroy_grid(const int handle);

/** \brief Returns the samoa grid of the current process.
 * If no such grid exists, a default grid of depth 0 will be created.
 *
 * \param handle            (in)handle to the grid
 * \param nsections         (out)number of sections in the grid
 * \param sections          (out)array of samoa sections that describe all grid entities and their relations
 *
 */
extern "C" void samoa_get_grid(const int handle, int& nsections, Samoa_section*& sections);


/** \brief Kernel to be executed on a triangular mesh in Sierpinski order.
 *
 * \param section_index     (in)index of the current section
 * \param cell_index        (in)index of the cell
 * \param refinement        (out)refinement flag: set to 1 for cell refinement, -1 for coarsening, 0 for no changes (default: 0)
 * \param data              (in)pointer to custom data
 *
 */
typedef void Samoa_kernel(const int section_index, const long long cell_index, char& refinement, void* data);


/** \brief Executes a kernel on a triangular mesh in Sierpinski order.
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
 * \param handle            (in)handle to the grid
 * \param data              (in)pointer to custom data
 * \param kernel            (in)kernel function
 *
 */
extern "C" void samoa_run_kernel(const int handle, void* data, Samoa_kernel* kernel);

/** \brief Describes a grid which is composed of multiple sections
 *
 * \param handle               Samoa handle
 * \param nsections            number of sections
 * \param sections             section data
 *
 */
struct Samoa_grid {
private:
    int handle;
    int nsections;
    Samoa_section* sections;

public:
    Samoa_grid(int depth) {
        handle = samoa_create_grid(1, depth);

        update();
    }

    ~Samoa_grid() {
        samoa_destroy_grid(handle);
    }

    Samoa_section& operator[](int idx) {
        assert(0 <= idx && idx < nsections);
        return sections[idx];
    }

    int get_nsections() {
        return nsections;
    }

    void update() {
        samoa_get_grid(handle, nsections, sections);
    }

    void run_kernel(void* data, Samoa_kernel* kernel) {
        samoa_run_kernel(handle, data, kernel);
    }
};
