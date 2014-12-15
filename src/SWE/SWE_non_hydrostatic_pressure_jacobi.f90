! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if 0
    module Swe_pressure_solver_jacobi
        !this is a dummy module for automated dependency analysis
        use SFC_edge_traversal
        use Samoa_swe
        use linear_solver
    end module
#endif

#if defined(_SWE)
#   define _solver              swe_pressure_solver_jacobi
#   define _solver_use          Samoa_swe

#   define _gv_node_size        _SWE_P_NODE_SIZE
#   define _gv_edge_size        _SWE_P_EDGE_SIZE
#   define _gv_cell_size        _SWE_P_CELL_SIZE

#   define _gm_A                swe_gm_A
#   define _gv_x                swe_gv_qp

#   define _gv_r                swe_gv_r
#   define _gv_rhs              swe_gv_rhs
#   define _gv_trace_A          swe_gv_mat_diagonal
#   define _gv_dirichlet        swe_gv_is_dirichlet_boundary

#   include "../Solver/Jacobi.f90"
#endif
