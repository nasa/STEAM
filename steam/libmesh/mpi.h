#ifndef __steam_mpi_h__
#define __steam_mpi_h__

//standard includes
//#include <fstream>

//libmesh includes
#include "libmesh/libmesh.h"
//#include "libmesh/mesh.h"
//#include "libmesh/tecplot_io.h"
//#include "libmesh/elem.h"
//#include "libmesh/mesh_refinement.h"
//#include "libmesh/mesh_modification.h"
//
////local includes
#include "mesh.h"
//#include "steam_libmesh_common.h"

void lib_split( libMesh::LibMeshInit* init, _Mesh* mesh );

#endif
