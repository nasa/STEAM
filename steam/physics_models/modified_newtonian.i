
%module modified_newtonian
%include "std_vector.i"
%include "std_string.i"
%include "typemaps.i"

%{
#include "modified_newtonian.h"
%}

%template(VectorDouble)          std::vector<double>;

%include "modified_newtonian.docs"
%include "modified_newtonian.h"
