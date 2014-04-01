! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_TESTS)
	!> this traversal routine just pushes memory.
	MODULE Tests_memory_cen
		use SFC_node_traversal

#		define _GT_NAME							Tests_memory_traversal_cen
#		define _GT_EDGES
#		define _GT_EDGES_TEMP
#		define _GT_NODES
#		define _GT_NODES_TEMP

#		include "SFC_generic_traversal_ringbuffer.f90"
	END MODULE

	!> this traversal routine just pushes memory.
	MODULE Tests_memory_cn
		use SFC_node_traversal

#		define _GT_NAME							Tests_memory_traversal_cn
#		define _GT_NODES
#		define _GT_NODES_TEMP

#		include "SFC_generic_traversal_ringbuffer.f90"
	END MODULE

	!> this traversal routine just pushes memory.
	MODULE Tests_memory_ce
		use SFC_node_traversal

#		define _GT_NAME							Tests_memory_traversal_ce
#		define _GT_EDGES
#		define _GT_EDGES_TEMP

#		include "SFC_generic_traversal_ringbuffer.f90"
	END MODULE

	!> this traversal routine just pushes memory.
	MODULE Tests_memory_c
		use SFC_node_traversal

#		define _GT_NAME							Tests_memory_traversal_c

#		include "SFC_generic_traversal_ringbuffer.f90"
	END MODULE

	MODULE Tests_memory
		use Tests_memory_cen
		use Tests_memory_cn
		use Tests_memory_ce
		use Tests_memory_c
	END MODULE
#endif

