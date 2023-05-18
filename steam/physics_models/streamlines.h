#ifndef __streamlines_h__
#define __streamlines_h__

//libmesh includes
#include "libmesh/equation_systems.h"
#include "libmesh/explicit_system.h"
#include "libmesh/tecplot_io.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/elem.h"

//local includes
#include "../libmesh/mesh.h"
#include "locate_stagnation_point.h"

//class definitions

enum StreamVectorMethodEnum { GEOMETRIC=0,
			      CBAERO=1 };

class StreamlineContainer;//forward declaration

class Streamline {

 public:

  // constructor with reference to mesh
  Streamline(StreamlineContainer& _container) :
    container(_container)
  {}

  // constructor with an element number creates the streamline
  Streamline(StreamlineContainer& _container,
	     unsigned int& seed) :
    container(_container)
  {
    create_from_seed_element(seed);
  }

  // function to create a streamline from a seed element
  void create_from_seed_element(const unsigned int& seed);

  // check to see if an element is contained in elem_set
  inline bool contains_element(const unsigned int& n) {
    return(elem_set.count(n));
  }

  inline const std::vector<unsigned int>& get_elems() {
    return(elems);
  }

  void write(std::ostream& stream, const std::string& label = "none");

 private:

  // private function members
  void find_next_element(unsigned int& el,
			 unsigned int& prev_el,
			 libMesh::Point& p,
			 const int& dir,
			 std::vector<unsigned int>& _elems);

  //private member functions
  void find_intersection_with_element_side(libMesh::Point& p,
					   const libMesh::Point& streamv,
					   const libMesh::Point& c,
					   libMesh::Point& p_to_c,
					   const libMesh::Point& sidev);

  void solve2x2 (const double& A00,
		 const double& A01,
		 const double& A10,
		 const double& A11,
		 const double& b0,
		 const double& b1,
		 double& x0,
		 double& x1);


  // private data members

  // vector to hold elements along a streamline in order from beginning to end
  std::vector<unsigned int> elems;

  /*
     set to hold list of elements in a streamline
     this is somewhat redundant to the elem vector
     but the set is more quickly searchable to see
     if an element is contained in the streamline
  */
  std::set<unsigned int> elem_set;

  // keep a reference to the parent vector containing all streamlines
  StreamlineContainer& container;

};

class StreamlineContainer : public std::vector<Streamline> {

 public:
  // default constructor
  StreamlineContainer(_Mesh& _mesh) :
    std::vector<Streamline>(),
    mesh(_mesh)
  {}

  //get a reference to the mesh
  inline _Mesh& get_mesh() const {
    return(mesh);
  }

  //add an element id to the terminating_elems set
  inline void insert_terminating_elem(const unsigned int& n) {
    terminating_elems.insert(n);
  }

  //test to see if an element id is contained in the terminating element set
  inline bool is_this_a_terminating_elem(const unsigned int& n) const {
    return(static_cast<bool>(terminating_elems.count(n)));
  }

  // clear the elem_stream_vectors vector
  inline void clear_elem_stream_vectors() {
    elem_stream_vectors.clear();
  }

  // check to see if elem_stream_vectors is populated
  bool elem_stream_vectors_is_populated();

  // get reference to an element's stream vector
  inline const libMesh::Point& get_elem_stream_vector(const unsigned int& n) const {
    return(elem_stream_vectors[n]);
  }

  // calculate elem_stream_vectors
  void populate_elem_stream_vectors(const StreamVectorMethodEnum& method,
				    const std::vector<double>& f);

  // creates streamlines
  void create_all_streamlines(const StagnationPointMethodEnum& stag_method,
			      const StreamVectorMethodEnum& stream_method,
			      const std::vector<double>& fs);

  //write streamline data to tecplot file for visualization
  void write(const std::string& filename);

  //clear
  void clear_data() {
    terminating_elems.clear();
    elem_stream_vectors.clear();
    this->clear();
  }

private:

  //private data members

  //set to hold list of elements that would terminate a streamline with either forward or backward integration
  std::set<unsigned int> terminating_elems;

  //store a const reference to the mesh
  _Mesh& mesh;

  //storage for element stream vectors
  std::vector<libMesh::Point> elem_stream_vectors;

};

//function prototypes

void create_stream_vectors(const StreamVectorMethodEnum& method,
			   _Mesh* mesh,
			   const std::vector<double>& fs,
			   std::vector<libMesh::Point>* stream_vectors,
			   const std::string& filename = "none");

void create_stream_vectors_geometric_method(_Mesh* mesh,
					    const std::vector<double>& fs_in,
					    std::vector<libMesh::Point>* stream_vectors,
					    const std::string& filename);

void create_stream_vectors_cbaero_method(_Mesh* mesh,
					 const std::vector<double>& fs_in,
					 std::vector<libMesh::Point>* stream_vectors,
					 const std::string& filename);

void create(StreamlineContainer* streamline_container,
	    //const StagnationPointMethodEnum& stag_method,
	    const StreamVectorMethodEnum& stream_method,
	    const std::vector<double>& fs);


void write(StreamlineContainer* streamline_container,
	   const std::string& filename);

void clear(StreamlineContainer* streamline_container);

StreamlineContainer* get_new_pointer(_Mesh* mesh);

void delete_pointer(StreamlineContainer* streamline_container);

#endif
