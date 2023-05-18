#ifndef __modified_newtonian_h__
#define __modified_newtonian_h__

//libmesh includes
#include "libmesh/equation_systems.h"
#include "libmesh/explicit_system.h"
#include "libmesh/tecplot_io.h"
#include "libmesh/numeric_vector.h"

//local includes
#include "../libmesh/mesh.h"

//function prototypes

void get_global_solution(_Mesh* mesh,
			 const double& V_mps,
			 const double& density_kgpm3,
			 const double& T_K,
			 const double& gamma,
			 const double& R_JpkgK,
			 const std::vector<double>& fs_in,
			 const std::string& filename = "none");

#endif
