#ifndef __mesh_h__
#define __mesh_h__

//standard includes
#include <fstream>

//libmesh includes
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/tecplot_io.h"
#include "libmesh/elem.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/mesh_modification.h"

/** 
* @brief _Mesh C/libMesh class
*
* Derived class of libMesh::Mesh to hold additional data that will be useful for our purposes
*/
class _Mesh: public libMesh::Mesh {

public:

  // constructor
  _Mesh(const libMesh::Parallel::Communicator& comm) : libMesh::Mesh(comm) {}

  // public member functions

  // clear vector that holds node normal data
  inline void clear_node_normal() {
    node_normal.clear();
  }

  // populate node normals
  void populate_node_normal();

  // check to see if node_normal is populated
  bool node_normal_is_populated();

  // write node normal pointcloud data to tecplot ascii file
  void write_node_normals(const std::string & filename);

  // get the normal vector at the specified node
  inline const libMesh::Point& get_node_normal(const unsigned int& n) const {
    return(node_normal[n]);
  }

  // clear vector that holds element normal data
  inline void clear_elem_normal() {
    elem_normal.clear();
  }

  // populate elem normals
  void populate_elem_normal();

  // check to see if elem_normal is populated
  bool elem_normal_is_populated();

  // write elem normal pointcloud data to tecplot ascii file
  void write_elem_normals(const std::string & filename);

  // get the normal vector at the specified element
  inline const libMesh::Point& get_elem_normal(const unsigned int& n) const {
    return(elem_normal[n]);
  }

  // clear vector that holds element centroid
  inline void clear_elem_centroid() {
    elem_centroid.clear();
  }

  // populate elem centroids
  void populate_elem_centroid();

  // check to see if elem_centroid is populated
  bool elem_centroid_is_populated();

  // write elem centroid pointcloud data to tecplot ascii file
  void write_elem_centroids(const std::string & filename);

  // get the centroid at the specified element
  inline const libMesh::Point& get_elem_centroid(const unsigned int& n) const {
    return(elem_centroid[n]);
  }

  void write_domain_decomposition(const std::string& filename);

  void uniformly_refine(const unsigned int& num_times);

  // clear node_to_elems vector
  inline void clear_node_to_elems() {
    node_to_elems.clear();
  }

  // populate node_to_elems vector
  void populate_node_to_elems();

  // check to see if node_to_elems is populated
  bool node_to_elems_is_populated();

private:

  //private data members

  /*
     vector to hold node normal data.
     this is populated with call to populate_node_normal()
     it is accessed on a nodal basis with call to get_node_normal(const unsigned int& n)
  */
  std::vector<libMesh::Point> node_normal;

  /*
     vector to hold node normal data.
     this is populated with call to populate_node_normal()
     it is accessed on a nodal basis with call to get_elem_normal(const unsigned int& n)
  */
  std::vector<libMesh::Point> elem_normal;


  /*
    vector to hold element centroid data
    this is populated with call to populate_elem_centroid()
    it is accessed on a nodal basis with call to get_elem_centroid(const unsigned int& n)
  */
  std::vector<libMesh::Point> elem_centroid;

  /*
    vector to hold list of elements that each node is part of
  */
  std::vector<std::vector<unsigned int> > node_to_elems;

};// end definition of class _Mesh

//function prototypes

/*
  BEGIN Functions that are wrapped by python Mesh class
*/

_Mesh* get_new_pointer(libMesh::LibMeshInit* init);

void delete_pointer (_Mesh* mesh);

void print_info(_Mesh* mesh);

void write(_Mesh* mesh,
	   const std::string& filename);

void read(_Mesh* mesh,
	  const std::string& filename);

void make (_Mesh* mesh,
	   const std::vector<double>& x,
	   const std::vector<double>& y,
	   const std::vector<double>& z,
	   const std::vector<unsigned int>& elems);

void get_data (_Mesh* mesh,
	       std::vector<double>& x,
	       std::vector<double>& y,
	       std::vector<double>& z,
	       std::vector<unsigned int>& elems);

unsigned int n_elem(_Mesh* mesh);

unsigned int n_nodes(_Mesh* mesh);

void populate_node_normal(_Mesh* mesh);

void write_node_normals(_Mesh* mesh, const std::string& filename);

void populate_elem_normal(_Mesh* mesh);

void write_elem_normals(_Mesh* mesh, const std::string& filename);

void populate_elem_centroid(_Mesh* mesh);

void write_elem_centroids(_Mesh* mesh, const std::string& filename);

void write_domain_decomposition(_Mesh* mesh, const std::string& filename);

void uniformly_refine(_Mesh* mesh, const unsigned int& num_times);

/*
  END Functions that are wrapped by python Init class
*/

void read_mesh_and_get_data (libMesh::LibMeshInit* init,
			     const std::string & filename,
			     std::vector<double>& x,
			     std::vector<double>& y,
			     std::vector<double>& z,
			     std::vector<unsigned int>& elems);

_Mesh* project_grid (_Mesh* source,
                     _Mesh* target);
#endif
