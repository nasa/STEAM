/*
   SWIG interface file to create system.py from system.C and system.h
   which provides a python interface to libMesh::EquationSystems object
*/

%module system

%include "std_string.i"
%include "std_vector.i"

%{
#include "steam_libmesh_common.h"
#include "system.h"
%}


/* Create a function to cast an ExplicitSystem* as a System* */
%inline %{
    libMesh::System* ExpToSys( libMesh::ExplicitSystem* exp_sys ) {
        return static_cast<libMesh::System*>( exp_sys );
    }
%}

%template(VectorDouble)          std::vector<double>;
%template(VectorInt)             std::vector<unsigned int>;
%template(VectorStr)             std::vector<std::string>;

%include "system.docs"
%include "system.h"
%include "steam_libmesh_common.h"
