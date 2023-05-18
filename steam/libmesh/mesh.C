#include "mesh.h"
#include "libmesh/plane.h"
//#include "libmesh/elem.h"

using namespace libMesh;

/*
   BEGIN Member functions of _Mesh class

   Derived class of libMesh::Mesh to hold additional data that will be useful for our purposes

*/

void _Mesh::populate_node_normal() {

  // check to see if we have already populated the node_normal vector. If so, return
  if (this->node_normal_is_populated()) return;

  // create temporary data container to help calculate normals
  std::vector<std::pair<unsigned int, Point> > _node_normal(this->n_nodes(),std::pair<unsigned int, Point>(0,Point(0,0,0)));

  //loop over elements

  //create objects that we do not want to reallocate each time through the loop
  Point v1;
  Point v2;
  Point normal;

  _Mesh::const_element_iterator       el     = this->active_elements_begin();
  const _Mesh::const_element_iterator end_el = this->active_elements_end();

  for ( ; el != end_el; ++el) {

    const Elem* elem = *el;

    //get vectors along tri's sides
    v1 = (elem->point(1) - elem->point(0));
    v2 = (elem->point(2) - elem->point(0));

    //get the unit normal by taking cross product
    normal = (v1.cross(v2)).unit();

    //loop over nodes on element and contribute this elements normal to the accumulated total for a node
    for (unsigned int n=0; n<elem->n_nodes(); n++) {
      _node_normal[elem->node_id(n)].first  += 1;
      _node_normal[elem->node_id(n)].second += normal;
    }//end loop over nodes

  }//end loop over elements

  //resize node_normal vector.  All processors will hold data for all nodes
  node_normal.resize(this->n_nodes());

  //loop over nodes and calculate average normal from accumulated element data
  for (unsigned int n=0; n<node_normal.size(); n++)
    node_normal[n] = (_node_normal[n].second/static_cast<double>(_node_normal[n].first)).unit();

}

/**
* Check it node normals are populated.
*
* @return bool if the node_normal.size() is non-zero.
*/
bool _Mesh::node_normal_is_populated() {

  if (node_normal.size() == this->n_nodes()) {

    // node normals already determined
    return(true);

  } else if (node_normal.size() != 0) {

    std::cerr << "ERROR: _Mesh::node_normal.size() != Mesh.n_nodes()" << std::endl;
    libmesh_error();

    return(false);
    
  } else {

    return(false);

  }

}

void _Mesh::write_node_normals(const std::string & filename) {

  // make sure that we have populated the node_normal vector
  this->populate_node_normal();
  
  if (this->processor_id() == 0) {

    std::ofstream stream(filename.c_str());

    if (stream.is_open()) {

      stream << "variables = x, y, z, nx, ny, nz\n";

      for (unsigned int n=0; n<node_normal.size(); n++) {

	stream << this->point(n)(0) << " "   // x
	       << this->point(n)(1) << " "   // y
	       << this->point(n)(2) << " "   // z
	       << node_normal[n](0) << " "   // nx
	       << node_normal[n](1) << " "   // nx
	       << node_normal[n](2) << "\n"; // nx

      }//end loop over nodes

      stream.close();

    } else {

      std::cerr << "Could not open output file for writing mesh node normals: " << filename << std::endl;
      libmesh_error();

    }//end if (stream.is_open())

  }//end if (this->processor_id() == 0)

  return;

}

void _Mesh::populate_elem_normal() {

  // check to see if we have already populated the elem_normal vector. If so, return
  if (this->elem_normal_is_populated()) return;

  elem_normal.resize(this->n_elem());

  //loop over elements

  //create objects that we do not want to reallocate each time through the loop
  Point v1;
  Point v2;
  Point normal;

  _Mesh::const_element_iterator       el     = this->active_elements_begin();
  const _Mesh::const_element_iterator end_el = this->active_elements_end();

  for ( ; el != end_el; ++el) {

    const Elem* elem = *el;

    //get vectors along tri's sides
    v1 = (elem->point(1) - elem->point(0));
    v2 = (elem->point(2) - elem->point(0));

    //get the unit normal by taking cross product
    elem_normal[elem->id()] = (v1.cross(v2)).unit();

  }//end loop over elements

}

bool _Mesh::elem_normal_is_populated() {

  if (elem_normal.size() == this->n_elem()) {

    // elem normals already determined
    return(true);

  } else if (elem_normal.size() != 0) {

    std::cerr << "ERROR: _Mesh::elem_normal.size() != Mesh.n_elem()" << std::endl;
    libmesh_error();

    return(false);
    
  } else {

    return(false);

  }

}

void _Mesh::write_elem_normals(const std::string & filename) {

  // make sure that we have populated the elem_normal vector
  this->populate_elem_normal();

  // make sure that we have populated the elem_centroid vector
  this->populate_elem_centroid();
  
  if (this->processor_id() == 0) {

    std::ofstream stream(filename.c_str());

    if (stream.is_open()) {

      stream << "variables = x, y, z, nx, ny, nz\n";

      for (unsigned int n=0; n<this->elem_normal.size(); n++) {

	stream << elem_centroid[n](0) << " "   // x
	       << elem_centroid[n](1) << " "   // y
	       << elem_centroid[n](2) << " "   // z
	       << elem_normal[n](0) << " "   // nx
	       << elem_normal[n](1) << " "   // nx
	       << elem_normal[n](2) << "\n"; // nx

      }//end loop over nodes

      stream.close();

    } else {

      std::cerr << "Could not open output file for writing mesh element normals: " << filename << std::endl;
      libmesh_error();

    }//end if (stream.is_open())

  }//end if (this->processor_id() == 0)

  return;

}

void _Mesh::populate_elem_centroid() {

  if (this->elem_centroid_is_populated()) return;

  elem_centroid.resize(this->n_elem());
  
  _Mesh::const_element_iterator       el     = this->active_elements_begin();
  const _Mesh::const_element_iterator end_el = this->active_elements_end();
  
  for ( ; el != end_el; ++el) {

    const Elem* elem = *el;

    elem_centroid[elem->id()] = (elem->point(0) + elem->point(1) + elem->point(2))/3;
    
  }//end loop over elements
  
}

bool _Mesh::elem_centroid_is_populated() {
/*
""" Check to ensure centroids are populated.

  Returns:
   (boolean) If the element centroids are defined.
"""
*/

  if (elem_centroid.size() == this->n_elem()) {

    // elem centroids already determined
    return(true);

  } else if (elem_centroid.size() != 0) {
    
    std::cerr << "ERROR: _Mesh::elem_centroid.size() != Mesh.n_elem()" << std::endl;
    libmesh_error();

    return(false);
    
  } else {
    
    return(false);
    
  }
  
}

void _Mesh::write_elem_centroids(const std::string & filename) {

  // make sure that we have populated the elem_centroid vector
  this->populate_elem_centroid();
  
  if (this->processor_id() == 0) {
    
    std::ofstream stream(filename.c_str());
    
    if (stream.is_open()) {
      
      stream << "variables = x, y, z\n";
      
      for (unsigned int n=0; n<elem_centroid.size(); n++) {
	
	stream << elem_centroid[n](0) << " "   // nx
	       << elem_centroid[n](1) << " "   // nx
	       << elem_centroid[n](2) << "\n"; // nx
	
      }//end loop over elements
      
      stream.close();
      
    } else {
      
      std::cerr << "Could not open output file for writing mesh element centroids: " << filename << std::endl;
      libmesh_error();
      
    }//end if (stream.is_open())
    
  }//end if (this->processor_id() == 0)
  
  return;
  
}

void _Mesh::write_domain_decomposition(const std::string& filename) {

  if (this->processor_id()==0) {

    std::ofstream stream;

    stream.open(filename.c_str());

    stream << "variables = \"x\" \"y\" \"z\" \"processor id\" \n" << std::endl;

    stream << "ZONE T=\"Domain Decomposition\""
	   << "\n DATAPACKING=BLOCK"
	   << "\n NODES=" << this->n_nodes()
	   << "\n ELEMENTS=" << this->n_elem()
	   << "\n ZONETYPE=FEQUADRILATERAL"
	   << "\n VARLOCATION=([1-3]=NODAL,[4]=CELLCENTERED)" << std::endl;

    //write x nodal coordinates
    for (unsigned int i=0; i<this->n_nodes(); i++ ) {
      stream << this->point(i)(0) << "\n";
    }

    //write y nodal coordinates
    for (unsigned int i=0; i<this->n_nodes(); i++ ) {
      stream << this->point(i)(1) << "\n";
    }

    //write z nodal coordinates
    for (unsigned int i=0; i<this->n_nodes(); i++ ) {
      stream << this->point(i)(2) << "\n";
    }

    for (unsigned int i=0; i<this->n_elem(); i++ ) {
      stream << this->elem(i)->processor_id() << "\n";
    }

    // Connectivity list
    for (unsigned int i=0; i<this->n_elem(); i++ ) {
      const Elem* elem = this->elem(i);
      for (unsigned int j=0; j<elem->n_nodes(); j++) {
	stream << elem->node(j)+1 << " ";
      }
      stream << elem->node(2)+1 << " "; //take care of collapsing one edge for triangle elements
      stream << "\n";
    }

  }

}

void _Mesh::uniformly_refine(const unsigned int& num_times) {
/*
""" Uniformly refine a mesh.

 This just wraps the libMesh method.

 Args:
    num_times (int) The number of times to refine each element.
"""
*/

  if (num_times > 0) {
    MeshRefinement mesh_refinement(*this);
    mesh_refinement.uniformly_refine(num_times);
    libMesh::MeshTools::Modification::flatten(*this);
  }
}

void _Mesh::populate_node_to_elems() {

  /*
    This function generates a vector for each node which contains 
    a list of elements it is a part of.  We'll hold the entire
    mapping on each processor (for now).
  */

  // first check to see if we have already done this
  if (this->node_to_elems_is_populated()) return;
  
  /*
    we know our vector needs to have an entry for each node
  */
  node_to_elems.resize(this->n_nodes());

  //loop over all nodes (not just local)
  for (unsigned int i=0; i<node_to_elems.size(); i++) {

    //loop over all elements (not just local)
    _Mesh::const_element_iterator       el     = this->active_elements_begin();
    const _Mesh::const_element_iterator end_el = this->active_elements_end();

    for ( ; el != end_el; ++el) {
      
      const Elem* elem = *el;

      //loop over element's nodes
      for (unsigned int j=0; j<elem->n_nodes(); j++) {
	// check if this element's node is our current node of interest
	if (elem->node_id(j) == i) {
	  // if so, add it to the list of elements that this node belongs to
	  node_to_elems[i].push_back(elem->id());
	  break;
	}
      }// end loop over element's nodes
      
    }// end loop over all elements
    
  }// end loop over all nodes
  
}

bool _Mesh::node_to_elems_is_populated() {

    if (node_to_elems.size() == this->n_nodes()) {

    // node to element mapping already determined
    return(true);

  } else if (node_to_elems.size() != 0) {

    std::cerr << "ERROR: _Mesh::node_to_elems.size() != Mesh.n_nodes()" << std::endl;
    libmesh_error();

    return(false);
    
  } else {

    return(false);

  }
  
}

/*
   END member functions of _Mesh class
*/

/*
  BEGIN Functions that are wrapped by python Mesh class
*/

_Mesh* get_new_pointer(LibMeshInit* init) {
  return(new _Mesh(init->comm())); //_Mesh is derived class created for this application so that we can carry extra data
}

void delete_pointer(_Mesh* mesh) {
  //std::cerr << "Hello (Mesh)!" << std::endl;
  delete mesh;
  //std::cerr << "Bye   (Mesh)!" << std::endl;
}

void print_info(_Mesh* mesh) {
  mesh->print_info();
}

void write(_Mesh* mesh, const std::string& filename) {
  mesh->write(filename);
}

void read(_Mesh* mesh, const std::string& filename) {
  mesh->read(filename);
}

void make (_Mesh* mesh,
	   const std::vector<double>& x,
	   const std::vector<double>& y,
	   const std::vector<double>& z,
	   const std::vector<unsigned int>& elems) {

  libmesh_assert( (x.size() == y.size()) && (x.size() == z.size()) );

  //add points to the mesh
  Point p;
  for (unsigned int i=0; i<x.size(); i++) {
    p(0) = x[i];
    p(1) = y[i];
    p(2) = z[i];
    mesh->add_point(p,i);
  }

  //add the elements to the mesh
  libmesh_assert( (elems.size()/4) == 0 );
  const unsigned int imax = elems.size()/4;
  for (unsigned int i=0; i<imax; i++) {//loop over elements

    const unsigned int j = i*4;//index shifting since our element data is vectorized

    // Build the required type and release it into our custody.
    Elem* elem = Elem::build(TRI3).release();

    //assign the nodes for the element
    elem->set_node(0) = mesh->node_ptr(elems[j]);
    elem->set_node(1) = mesh->node_ptr(elems[j+1]);
    elem->set_node(2) = mesh->node_ptr(elems[j+2]);

    //set the subdomain id
    elem->subdomain_id() = cast_int<subdomain_id_type>(elems[j+3]);

    // Add the element to the mesh
    elem->set_id(i);
    mesh->add_elem(elem);

  }

  mesh->set_mesh_dimension(2);

  /**
      this will keep libMesh from renumbering nodes in prepare_for_use() so that the node
      numbering is preserved as what was input into this function
  **/
  mesh->allow_renumbering(false);

  //find connectivity and partition
  mesh->prepare_for_use();

}

void get_data (_Mesh* mesh,
	       std::vector<double>& x,
	       std::vector<double>& y,
	       std::vector<double>& z,
	       std::vector<unsigned int>& elems) {

  x.resize (mesh->n_nodes());
  y.resize (x.size());
  z.resize (x.size());

  elems.resize (mesh->n_elem() * 4);

  for (unsigned int i=0; i < x.size(); i++) {
    x[i] = mesh->point (i)(0);
    y[i] = mesh->point (i)(1);
    z[i] = mesh->point (i)(2);
  }

  for (unsigned int i=0; i < mesh->n_elem(); i++) {
    Elem* el=mesh->elem_ptr(i);
    const unsigned int k=i*4;
    for (unsigned int j=0; j < el->n_nodes(); j++) {
      elems[k+j]=el->node_id(j);
    }
    elems[k+3]=el->subdomain_id();
  }
}

unsigned int n_elem(_Mesh* mesh) {
  return(static_cast<unsigned int>(mesh->n_elem()));
}

unsigned int n_nodes(_Mesh* mesh) {
  return(static_cast<unsigned int>(mesh->n_nodes()));
}

void populate_node_normal(_Mesh* mesh) {
  mesh->populate_node_normal();
}

void write_node_normals(_Mesh* mesh, const std::string& filename) {
  mesh->write_node_normals(filename);
}

void populate_elem_normal(_Mesh* mesh) {
  mesh->populate_elem_normal();
}

/**
* Write mesh normals to disk.
*
* @param mesh is pointer to _Mesh to write.
* @param filename is the output filename.
*
* This just wraps the libMesh method.
*/
void write_elem_normals(_Mesh* mesh, const std::string& filename) {
  mesh->write_elem_normals(filename);
}

void populate_elem_centroid(_Mesh* mesh) {
  mesh->populate_elem_centroid();
}
   
void write_elem_centroids(_Mesh* mesh, const std::string& filename) {
  mesh->write_elem_centroids(filename);
}

/**
* Write domain decomposition to disk.
*
* @param mesh is pointer to _Mesh to write.
* @param filename is the output filename.
*
* This just wraps the libMesh method.
*/
void write_domain_decomposition(_Mesh* mesh, const std::string& filename) {
  mesh->write_domain_decomposition(filename);
}

/**
* Uniformly refine a mesh.
*
* @param mesh is pointer to _Mesh to refine.
* @param num_times is the number of times to refine each element.
*
* This just wraps the libMesh method.
*/
void uniformly_refine(_Mesh* mesh, const unsigned int& num_times) {
  mesh->uniformly_refine(num_times);
}

/*
  END Functions that are wrapped by python Init class
*/

/*
  BEGIN functions directly callable from python without being wrapped by a python class
*/

void read_mesh_and_get_data (LibMeshInit* init,
			     const std::string & filename,
			     std::vector<double>& x,
			     std::vector<double>& y,
			     std::vector<double>& z,
			     std::vector<unsigned int>& elems) {
/**
* """Read grid from disk and return data to Python.
*
* Args:
*     init (:obj:`libMesh`,pointer): is pointer to libMesh instance.
*     filename (:obj:`string`,input): mesh filename to be read
*     x (:obj:`list`,output): list of floats containing x-coordinate for each node
*     y (:obj:`list`,output): list of floats containing y-coordinate for each node
*     z (:obj:`list`,output): list of floats containing z-coordinate for each node
*     elems (:obj:`list`,output): list of integers containing 3 node numbers for right-handed triangle element plus a 4th integer for component id.  List will be 4*(Number of Elements) long.
*  """
*/

  _Mesh mesh(init->comm());
  mesh.read(filename);
  get_data(&mesh, x, y, z, elems);

}

_Mesh* project_grid (_Mesh* source, _Mesh* target)  {
/**
*""" Projects source grid onto target grid and returns result.
*
*  Args:
*     source (:obj:`_Mesh`,pointer): _Mesh pointer to what will be projected.
*      target (:obj:`_Mesh`,pointer): _Mesh pointer to the projection target.
*
*    Returns:
*    :obj:`_Mesh`,pointer: _Mesh pointer to a projected copy of source.
*
* """
*/

   _Mesh* output = new _Mesh(*source);
   //output->print_info();

    //source->write("old_grid.dat");
    //target->write("target_grid.dat");

    {// projection scope

      //loop over nodes of mesh that we want to project
      for (unsigned int na=0; na<output->n_nodes(); na++) {

	Node& pta = output->node(na);

	Real  minimum_face_distance = std::numeric_limits<Real>::max();
	Real  minimum_edge_distance = std::numeric_limits<Real>::max();
	Real  minimum_node_distance = std::numeric_limits<Real>::max();
	Point best_face_projection;
	Point best_edge_projection;
	Point best_node_projection;
	bool  found_suitable_face_projection = false;
	bool  found_suitable_edge_projection = false;
	Real  distance;

	{// projection to faces and edges scope

	  //loop over elements on target mesh
	  _Mesh::const_element_iterator       el     = target->elements_begin();
	  const _Mesh::const_element_iterator end_el = target->elements_end();

	  for ( ; el != end_el; ++el) {

	    const Elem* elem = *el;

	    //make plane from this 2D element. use the first the nodes
	    const Plane plane(elem->point(0),elem->point(1),elem->point(2));

	    //project our node to the plane
	    Point closest_pt = plane.closest_point(pta);

	    //get the projection distance that we are trying to minimize
	    distance = (pta - closest_pt).size_sq();

	    if ((distance < minimum_face_distance) &&
		(elem->contains_point(closest_pt))) {

	      minimum_face_distance = distance;
	      best_face_projection  = closest_pt;
	      found_suitable_face_projection = true;

	    }

	    //loop over element sides
	    for (unsigned int s=0; s<elem->n_sides(); s++) {

	      //create a vector from node 0 to node 1 along this side
	      const Point va = elem->side(s)->point(1) - elem->side(s)->point(0);

	      //create a vector from node 0 to our point of interest
	      const Point vb = pta - elem->side(s)->point(0);

	      //get the vector projection of our point to the line that contains our edge
	      const Point vc = va*vb/va.size_sq()*va;

	      //now get our projected point in cartesian space
	      closest_pt = elem->side(s)->point(0) + vc;

	      //get the projection distance that we are trying to minimize
	      distance = (pta - closest_pt).size_sq();

	      if ((distance < minimum_edge_distance) &&
		  (elem->side(s)->contains_point(closest_pt))) {

		minimum_edge_distance = distance;
		best_edge_projection  = closest_pt;
		found_suitable_edge_projection = true;

	      }

	    }//end loop over sides of element

	  }//end loop over elements in target mesh

	}//end projection to faces and edges scope

	{//projection to nodes scope

	  //loop over nodes on target mesh
	  for (unsigned int nb=0; nb<target->n_nodes(); nb++) {

	    Node& ptb = target->node(nb);

	    distance = (ptb - pta).size_sq();

	    if (distance < minimum_node_distance) {

	      minimum_node_distance = distance;
	      best_node_projection  = ptb;

	    }

	  }//end loop over nodes on target mesh

	}//end projection to nodes scope

	//decide which projection is best to take
	if ((found_suitable_face_projection) &&
	    (minimum_face_distance < minimum_edge_distance) &&
	    (minimum_face_distance < minimum_node_distance)) {

	  pta = best_face_projection;

	} else if ((found_suitable_edge_projection) &&
		   (minimum_edge_distance < minimum_face_distance) &&
		   (minimum_edge_distance < minimum_node_distance)) {

	  pta = best_edge_projection;

	} else {

	  pta = best_node_projection;

	}

      }//end loop over nodes to project

    }// end projection scope

    {// output scope

//      std::cout << "\nWriting final projected grid file: " << Filename << std::flush;
//      source.write(project_grid);
//      std::cout << "\nComplete.\n" << std::endl;
//      source->write("new_grid.dat");
//        soln = source;

    }// end output scope

   return output;

  }



/*
  END functions directly callable from python without being wrapped by a python class
*/
