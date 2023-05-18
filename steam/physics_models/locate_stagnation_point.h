#ifndef __locate_stagnation_point_h__
#define __locate_stagnation_point_h__

//libmesh includes


//local includes
#include "../libmesh/mesh.h"

//function prototypes
enum StagnationPointMethodEnum { MOST_FORWARD_ELEM=0 };

void locate_stagnation_points(const StagnationPointMethodEnum& method,
			      _Mesh* mesh,
			      const std::vector<double>& fs_in,
			      std::vector<unsigned int>& stag_elems,
			      const std::string& filename = "none");

void most_forward_elem(_Mesh* mesh,
		       const std::vector<double>& fs_in,
		       std::vector<unsigned int>& stag_elems,
		       const std::string& filename);

#endif
