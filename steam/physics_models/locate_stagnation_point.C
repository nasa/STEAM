#include "locate_stagnation_point.h"

using namespace libMesh;

void locate_stagnation_points(const StagnationPointMethodEnum& method,
			      _Mesh* mesh,
			      const std::vector<double>& fs_in,
			      std::vector<unsigned int>& stag_elems,
			      const std::string& filename) {

  if (method == MOST_FORWARD_ELEM) {

    most_forward_elem(mesh,fs_in,stag_elems,filename);

  } else {

    std::cerr << "ERROR: Invalid choice for stagnation point method." << std::endl;
    libmesh_error();

  }

}

void most_forward_elem(_Mesh* mesh,
		       const std::vector<double>& fs_in,
		       std::vector<unsigned int>& stag_elems,
		       const std::string& filename) {

  //make sure centroid data is populated
  mesh->populate_elem_centroid();

  // freestream directed unit vector
  Point fs1(fs_in[0],fs_in[1],fs_in[2]);

  //minimum value that we'll keep track of
  double min = std::numeric_limits<Real>::max();

  //node number associated with minium
  unsigned int min_n;

  //scalar projection
  double proj;

  for (unsigned int n=0; n<mesh->n_elem(); n++) {

    proj = fs1*mesh->get_elem_centroid(n);

    if (proj < min) {
      min   = proj;
      min_n = n;
    }

  }

  const double local_min(min);

  //fine the minimum value among all processors
  mesh->comm().min(min);

  if (local_min != min)
    //this processor did not win the search for the stagnation element
    min_n = 0;

  /*
     now all should have min_n = 0 except the processor that found
     the stagntion element.  The processor that found the stagnation element
     will have the stagnation element number stored in min_n, which should be
     greater than zero.  Now a parallel max() should get the stagnation
     element number out to all processors.  More than one processor may have
     found the stagnation element, so we have to do a max() instead of a
     broadcast() to prevent two processors from stepping on each other
  */

  mesh->comm().max(min_n);

  if (mesh->processor_id()==0) {

    if (filename != "none") {

      std::ofstream stream;

      stream.open(filename.c_str());

      if (stream.is_open()) {

	stream << "variables = x, y, z\n";
	stream << mesh->get_elem_centroid(min_n)(0) << " "
	       << mesh->get_elem_centroid(min_n)(1) << " "
	       << mesh->get_elem_centroid(min_n)(2) << std::endl;

      } else {

	std::cerr << "ERROR: Could not open output file "
		  << filename
		  << " in most_forward_elem() stagnation point finding routine."
		  << std::endl;

	libmesh_error();

      }

    }

  }

  stag_elems.resize(1,min_n);

}
