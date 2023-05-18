
%module streamlines
%include "std_vector.i"
%include "std_string.i"
%include "typemaps.i"

%{
#include "streamlines.h"
%}

%template(VectorDouble)          std::vector<double>;

%include "streamlines.docs"
%include "streamlines.h"
