#ifndef __steam_libmesh_common_h__
#define __steam_libmesh_common_h__

/*
  This file contains "common" definitions that are used by both 
  equation_systems.* and system.*.  It is required because SWIG doesn't like
  the circular requirements that result from including both system.h and 
  equation_systems.h together.  Instead, both of those files can include this
  one and not each other.
*/

//libmesh includes
#include "libmesh/equation_systems.h"
#include "libmesh/explicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/tecplot_io.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/vtk_io.h"

struct ElemNodeStorage{
    enum Enum{ ELEMENTS = 0, NODES = 1 };
};

struct ModelType{
    enum Enum{ MULTIPLICATIVE = 0, ADDITIVE = 1, REPLACEMENT = 2 };
};

struct OutputType{
    enum Enum { TECPLOT_DAT = 0, TECPLOT_PLT = 1, EXODUS = 2, VTK_OUT = 3};
};

struct SolnSplit{
    enum Enum { LOCAL_SOLN = 0, ENTIRE_SOLN = 1 };
};

#endif
