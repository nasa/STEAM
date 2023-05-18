/* file: model.C */
#include <iostream>
#include <mpi.h>
#include "mpi.h"

using namespace libMesh;

void lib_split( LibMeshInit* init, _Mesh* mesh ){
  int size, rank;
  size = init -> comm().size();
  rank = init -> comm().rank();

  int local_nodes, total_nodes;
  local_nodes = int(mesh->n_local_nodes());
  total_nodes = int(mesh->n_nodes());
  printf( "Hello from C++ rank %d. n_nodes: %d.\n", rank, total_nodes );
  printf( "Hello from C++ rank %d. n_local_nodes: %d.\n", rank, local_nodes );

}

//void root( _Mesh* mesh, LibMeshInit* init){
//  /*  Main function for root rank to run -- First we'll just receive the mesh
//   *  from Python and split it into the number of pawns.  Then distribute the
//   *  mesh segments and let the pawns work.
//   */
//
//  // Recieve the mesh from Python
//
//  // Split the mesh
//  //   Example in /aerolab/aamar/aamar/codes/char/master/src/init.C
//  //   -- search for mesh.prepare_for_use()
//  mesh-> partition( n_pawns );
//
//}
//
//void pawn(){
//  /*  Main function for pawn ranks to run -- First we'll just report how
//   *  many nodes are present in the local segment of the mesh.
//   */
//
//  // Communicate to the root process that we're ready for work
//
//    return
//}
