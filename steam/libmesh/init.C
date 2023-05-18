#include "init.h"

using namespace libMesh;

/*
  BEGIN Functions that are wrapped by python Init class
*/

/*
  initializes libmesh through creation of LibMeshInit object
  Argument Input:  command line arguments
  Argument Output:
  Return: pointer to LibMeshInit object that python can use
*/
//LibMeshInit* get_new_pointer(int& argc, char **argv, MPI_Comm comm){
LibMeshInit* get_new_pointer(std::vector<std::string> py_arg, MPI_Comm comm){
  // Create new argc, argv from python arguments, py_arg
  std::vector<char*> argv;
  for (const auto &arg : py_arg)
    argv.push_back( (char*) arg.data() );
  argv.push_back( nullptr ); 
  int argc = argv.size() - 1;
  
  return(new LibMeshInit(argc, argv.data(), comm));
}

/*
  deletes LibMeshInit object being pointed to
*/
void delete_pointer(LibMeshInit* init) {
  //std::cerr << "Hello (Init)!" << std::endl;
  delete init;
  //std::cerr << "Bye   (Init)!" << std::endl;
}

/*
  END Functions that are wrapped by python Init class
*/
