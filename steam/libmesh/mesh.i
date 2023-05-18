/*
   SWIG interface file to create mesh.py from mesh.C and mesh.h
   which provides a python interface to libMesh::Mesh object
*/

%module mesh

%include "std_vector.i"
%include "std_string.i"
%include "typemaps.i"

%{
#include "mesh.h"
%}

%template(VectorDouble)          std::vector<double>;
%template(VectorInt)             std::vector<unsigned int>;

%newobject get_new_pointer;
%delobject delete_pointer;

%include "mesh.docs"
%include "mesh.h"
