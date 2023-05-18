// file: mpi.i
%module mpi

//%include "std_string.i"
//%include "std_vector.i"

%{
#include <stdio.h>
#include <mpi.h>
#include "libmesh/libmesh.h"
#include "mpi.h"
%}

//%include "libmesh/libmesh.h"
%include "mpi.h"
//%include mpi4py/mpi4py.i
//%template(VectorStr)        std::vector<std::string>;

//%mpi4py_typemap(Comm, MPI_Comm);
//void sayhello();
