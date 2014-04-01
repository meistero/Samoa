! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_DARCY)
	MODULE Darcy_edge_dummy
		use SFC_edge_traversal

        type num_traversal_data
        end type

#		define _GT_NAME							t_darcy_edge_dummy_traversal
#		define _GT_NODES

#		include "SFC_generic_traversal_ringbuffer.f90"
	END MODULE
#endif
