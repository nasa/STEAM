/*
   SWIG interface file to create equation_systems.py from equation_systems.C and equation_systems.h
   which provides a python interface to libMesh::EquationSystems object
*/

%module equation_systems

%include "std_string.i"
%include "std_vector.i"

%{
#include "equation_systems.h"
#include "steam_libmesh_common.h"
%}

/* Create a function to cast a MeshBase* as a _Mesh* */
%inline %{
    _Mesh* BaseToMesh( libMesh::MeshBase* mb ) {
        return static_cast<_Mesh*>( mb );
    }
%}

%template(VectorDouble)          std::vector<double>;
%template(VectorInt)             std::vector<unsigned int>;
%template(VectorStr)             std::vector<std::string>;

%newobject get_new_pointer;
%delobject delete_pointer;

%include "equation_systems.docs"
%include "equation_systems.h"
%include "steam_libmesh_common.h"
