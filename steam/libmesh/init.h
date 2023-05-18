#ifndef __init_h__
#define __init_h__

#include "libmesh/libmesh.h"

/*

  C++ side of wrapper interface for libMesh::LibMeshInit

*/

//function prototypes

/*
  BEGIN Functions that are wrapped by python Init class
*/

/* Old get_new_pointer prototype:
libMesh::LibMeshInit* get_new_pointer(int& argc, char **argv, 
                                      MPI_Comm comm=MPI_COMM_WORLD);
*/
libMesh::LibMeshInit* get_new_pointer( std::vector<std::string> py_arg, 
                                       MPI_Comm comm=MPI_COMM_WORLD );

void                  delete_pointer(libMesh::LibMeshInit* init);

/*
  END Functions that are wrapped by python Init class
*/

#endif
