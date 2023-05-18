/*
   SWIG interface file to create init.py from init.C and init.h
   which provides a python interface to libMesh::LibMeshInit object
*/

%module init
%include typemaps.i
%include "std_string.i"
%include "std_vector.i"

%{
#include "init.h"
#include <mpi.h>
%}

%include mpi4py/mpi4py.i
%template(VectorStr) std::vector<std::string>;
%mpi4py_typemap(Comm, MPI_Comm);
/*

   derived from information at:

   http://www.math.uiuc.edu/~gfrancis/illimath/windows/aszgard_mini/szg/python/src/PyTypemaps.i

   PyString_AsString() is replaced with PyUnicode_AsUTF8() to port from Python 2 to Python 3

*/

/*
Conversion from python list to argc, argv now happens in get_new_pointer
// typemap to convert python list of strings (1 argument) to c++ int, char ** (2 arguments)
%typemap(in) (int& argc, char **argv) (int size) {
  if (PyList_Check($input)) {
    size = PyList_Size($input);
    int i = 0;
    $1 = &size;
    $2 = (char **) malloc((size+1)*sizeof(char *));
    for (i = 0; i < size; ++i)
      $2[i] = PyUnicode_AsUTF8(PyList_GetItem($input,i));
    $2[i] = 0;
  } else {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
 }
%typemap(freearg) (int& argc, char **argv) (int size) {
  if ($2) free($2);
 }
 */

%newobject get_new_pointer;
%delobject delete_pointer;

%include "init.docs"
%include "init.h"
