! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_TESTS)
	MODULE Tests_node_dummy
		use SFC_node_traversal

#		define _GT_NAME							tests_node_dummy_traversal
#		define _GT_CELL							.false.
#		define _GT_CELL_TEMP					.false.
#		define _GT_EDGES						.false.
#		define _GT_EDGES_TEMP					.false.
#		define _GT_NODES						.true.
#		define _GT_NODES_TEMP					.false.

#		include "SFC_generic_traversal_ringbuffer.f90"

		!undefine macros to avoid compiler warnings
#		undef	_GT_NAME
#		undef	_GT_CELL
#		undef	_GT_CELL_TEMP
#		undef	_GT_EDGES
#		undef	_GT_EDGES_TEMP
#		undef	_GT_NODES
#		undef	_GT_NODES_TEMP
	END MODULE
#endif
