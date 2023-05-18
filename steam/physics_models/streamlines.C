#include "streamlines.h"

using namespace libMesh;

/*

  BEGIN Streamline class functions

*/

void Streamline::create_from_seed_element(const unsigned int& seed) {

  // these are recursive calls to integrate along streamline

  //backward integration
  Point p = container.get_mesh().get_elem_centroid(seed);
  std::vector<unsigned int> elems_b; //elements collected during backward integration
  unsigned int el = seed;
  unsigned int prev_el = seed;
  find_next_element(el,prev_el,p,-1,elems_b);

  //forward integration (reset seed location first)
  p = container.get_mesh().get_elem_centroid(seed);
  std::vector<unsigned int> elems_f; //elements collected during forward integration
  el = seed;
  prev_el = seed;
  find_next_element(el,prev_el,p,1,elems_f);

  // lets merge the element lists for the backward and forward integration
  elems.resize(elems_b.size() + elems_f.size() - 1);//-1 takes away the overlap (seed) element that is in both lists
  unsigned int j = 0;
  for (std::vector<unsigned int>::reverse_iterator i = elems_b.rbegin(); i != elems_b.rend(); ++i ) {
    elems[j] = *i;
    j++;
  }

  for (unsigned int i=1; i<elems_f.size(); i++) {
    elems[j] = elems_f[i];
    j++;
  }

  libmesh_assert(elem_set.size() == elems.size());

}

void Streamline::find_next_element(unsigned int& el,
				   unsigned int& prev_el,
				   Point& p,
				   const int& dir,
				   std::vector<unsigned int>& _elems) {

  if (elem_set.count(el) && (_elems.size() > 0)) return;  // the streamline has doubled back on itself

  // otherwise we can add the current element
  elem_set.insert(el);
  _elems.push_back(el);

  if (container.is_this_a_terminating_elem(el)) return; // any streamline must terminate if it hits one of these elements

  // get a const reference to the 'c'urrent 'el'ement
  const Elem& cel = container.get_mesh().elem_ref(el);

  // get the element's stream vector and account for integration direction
  const Point  sv     = dir*container.get_elem_stream_vector(el);

  //get the element's normal
  const Point& normal = container.get_mesh().get_elem_normal(el);

  //find the two 's'mallest 'c'omponents of the 'n'ormal 'v'ector
  //these will be the largest components of the tangent vectors
  unsigned int scnv0 = 0;//assumes that normal(2) is the largest component
  unsigned int scnv1 = 1;//assumes that normal(2) is the largest component

  if ( (normal(0) >= normal(1)) && (normal(0) >= normal(2)) ) {

    scnv0 = 1;
    scnv1 = 2;

  } else if ( (normal(1) >= normal(0)) && (normal(1) >= normal(2))) {

    scnv0 = 0;
    scnv1 = 2;

  }

  unsigned int side;//this will store the side number that we intersect
  double a, b; //magnitude of stream vector components in non-orthogonal basis
  Point tau0; //first non-orthogonal basis vector defined by diretion from integration point (p) to first corner of interest
  Point tau1; //second non-orthogonal basis vector defined by diretion from integration point (p) to second corner of interest
  unsigned int sp1; //side plus 1. (side, sp1) = (0, 1) OR (1, 2) OR (2, 0) hence the funky logic to define sp1
  for (side=0; side<3; side++) {

    if ((cel.neighbor(side))->id() == prev_el) continue;

    sp1 = (side == 2) ? 0 : side+1;

    tau0 = cel.point(side) - p;
    tau1 = cel.point(sp1)  - p;

    solve2x2(tau0(scnv0),tau1(scnv0),tau0(scnv1),tau1(scnv1),sv(scnv0),sv(scnv1),a,b);

    if ((a > 0.0) && (b > 0.0)) {

      find_intersection_with_element_side(p,sv,cel.point(side),tau0,(cel.point(sp1)-cel.point(side)).unit());

      break;

    }

  }

  if (side == 3) {
    //std::cerr << "ERROR: Streamline integration failed to traverse element: " << el << std::endl;
    return;
  }

  prev_el = el;
  el = (cel.neighbor(side))->id();

  find_next_element(el,prev_el,p,dir,_elems);

  return;

}

void Streamline::find_intersection_with_element_side(Point& p,  //integration initiation point
						     const Point& streamv, //stream vector (unit)
						     const Point& c, //triangle corner on side we know we intersect
						     Point& p_to_c, //vector initiation point to corner (c - p)
						     const Point& sidev) {//vector along side pointed away from c (unit)

  /*
                                                  . \
                                               .   | \
                                            .      |  \
                                         .         |   \
                                      .            |    \E
                                   .               |     \
                                .                  |      \    __
             streamv __     .                      |       \  |\ sidev (comes in as unit)
  (comes in as unit) . |  .                        |        \   \
                  .    .   theta                   |     phi \
                      p---------------------------------------c
                                   A                    B
                      --->
		      p_to_c

                     |_____________________D__________________|
                     |                                        |
   */

  const double D = p_to_c.size();

  p_to_c = p_to_c.unit();

  const double cos_theta =  streamv*p_to_c;
  const double cos_phi   =   -sidev*p_to_c;

  const double sin_theta = std::sqrt(1 - cos_theta*cos_theta);
  const double sin_phi   = std::sqrt(1 - cos_phi*cos_phi);

  const double E = D/((sin_phi*cos_theta/sin_theta) + cos_phi);

  p = c + E*sidev;

}

void Streamline::solve2x2 (const double& A00,
			   const double& A01,
			   const double& A10,
			   const double& A11,
			   const double& b0,
			   const double& b1,
			   double& x0,
			   double& x1) {

  const double det = 1.0/(A00*A11 - A01*A10);

  x0 = det*( A11*b0 - A01*b1);
  x1 = det*(-A10*b0 + A00*b1);

}

void Streamline::write(std::ostream& stream, const std::string& label) {

  if (label != "none") {

    stream << "Zone T=\"SL_" << label << "\"\n";

  } else {

    stream << "Zone T=\"SL\"\n";

  }

  //loop over elements in streamline
  for (unsigned int j=0; j<elems.size(); j++) {

    stream << (container.get_mesh().get_elem_centroid(elems[j]))(0) << " "
	   << (container.get_mesh().get_elem_centroid(elems[j]))(1) << " "
	   << (container.get_mesh().get_elem_centroid(elems[j]))(2) << " "
	   << elems[j] << "\n";
  }//end loop over elements in streamline

  stream << std::flush;

}

/*

  END Streamline class functions

*/

/*

  BEGIN StreamlineContainer class functions

*/

bool StreamlineContainer::elem_stream_vectors_is_populated() {

  if (elem_stream_vectors.size() == mesh.n_elem()) {

    // node normals already determined
    return(true);

  } else if (elem_stream_vectors.size() != 0) {

    std::cerr << "ERROR: _Mesh::elem_stream_vectors.size() != Mesh.n_elem()" << std::endl;
    libmesh_error();

    return(false);

  } else {

    return(false);

  }

}

void StreamlineContainer::populate_elem_stream_vectors(const StreamVectorMethodEnum& method,
						       const std::vector<double>& fs) {
  create_stream_vectors(method,&mesh,fs,&elem_stream_vectors);
}

void StreamlineContainer::create_all_streamlines(const StagnationPointMethodEnum& stag_method,
						 const StreamVectorMethodEnum& stream_method,
						 const std::vector<double>& fs) {

  //clear out any existing data
  clear_data();

  {
    std::vector<unsigned int> stag_elems;
    locate_stagnation_points(stag_method,&mesh,fs,stag_elems);
    for (unsigned int i=0; i<stag_elems.size(); i++)
      terminating_elems.insert(stag_elems[i]);
  }

  if (!elem_stream_vectors_is_populated())
    populate_elem_stream_vectors(stream_method,fs);

  //loop over elements to seed
  for (unsigned int el=0; el<mesh.n_elem(); el++) {

    if (terminating_elems.count(el)) continue;

    bool in_another_streamline = false;
    for (unsigned int s=0; s<this->size(); s++) {
      if ((*this)[s].contains_element(el)) {
	in_another_streamline = true;
	break;
      }
    }
    if (in_another_streamline) continue;

    // this element is an eligible seed because it is not a terminating element
    // and it is not already contained in another streamline

    this->push_back(Streamline(*this,el));

  }// end of loop over elements

}

void StreamlineContainer::write(const std::string& filename) {

  if (this->size() == 0) {
    std::cerr << "WARNING: Cannot write streamlines to file because streamlines have not been generated." << std::endl;
    return;
  }

  std::ofstream stream;

  stream.open(filename.c_str());

  if (!stream.is_open()) {
    std::cerr << "ERROR: Could not open file " << filename << " for wirting streamlines." << std::endl;
    libmesh_error();
  }

  stream << "variables = x, y, z, el\n";

  //loop over streamlines
  for (unsigned int i=0; i<this->size(); i++)
    (*this)[i].write(stream,std::to_string(i));
  
}

/*

  END StreamlineContainer class functions

*/

void create_stream_vectors(const StreamVectorMethodEnum& method,
			   _Mesh* mesh,
			   const std::vector<double>& fs,
			   std::vector<Point>* stream_vectors,
			   const std::string& filename) {

  if (method == GEOMETRIC) {

    create_stream_vectors_geometric_method(mesh,fs,stream_vectors,filename);

  } else if (method == CBAERO) {

    create_stream_vectors_cbaero_method(mesh,fs,stream_vectors,filename);

  } else {

    std::cerr << "ERROR: Invalid choice for stream vector method." << std::endl;
    libmesh_error();

  }

}

void create_stream_vectors_geometric_method(_Mesh* mesh,
					    const std::vector<double>& fs_in,
					    std::vector<libMesh::Point>* stream_vectors,
					    const std::string& filename) {

  //create EquationSystems object
  EquationSystems es(*mesh);

  //add explicit system to hold the streamline solution
  es.add_system<ExplicitSystem>("Stream Vectors");

  //get a reference to the explicit system we just made
  ExplicitSystem& system = es.get_system<ExplicitSystem>(0);

  //add the variables for 3 components of streamline vector
  system.add_variable("nx",CONSTANT,XYZ);
  system.add_variable("ny",CONSTANT,XYZ);
  system.add_variable("nz",CONSTANT,XYZ);

  /*
    initialize the equation systems obect for use
    this is a trivial (in this instance) but necessary call to prepare the system for use
    It will initialize the solution to zero, but it doesn't matter since we will overwrite the solution
  */
  es.init();

  /*
    this method needs surface normals
    this call will skip calculation of element normal data if it has already been done and stored for this mesh object
   */
  mesh->populate_elem_normal();

  //initialize the freestream vector from user input and make it a unit vector
  Point fs(fs_in[0],fs_in[1],fs_in[2]);

  //objects that we do not want to re-allocate each time through the loop
  Point streamline;

  //loop over local nodes
  MeshBase::const_element_iterator       el     = mesh->active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_local_elements_end();

  for ( ; el != end_el; ++el) {

    const Elem* elem = *el;

    const Point& normal = mesh->get_elem_normal(elem->id());

    streamline = (fs - ((fs*normal)*normal)).unit();

    system.solution->set (elem->dof_number(0,0,0), streamline(0));
    system.solution->set (elem->dof_number(0,1,0), streamline(1));
    system.solution->set (elem->dof_number(0,2,0), streamline(2));

  }//end loop over nodes

  //close solution vector to prevent further modification and prepare if for parallel communication
  system.solution->close();

  //sync solution amongst processors
  system.update();

  //write tecplot file if user wants.  currently this is hardcoded for ascii since its intent is for debugging
  if (filename != "none" ) {

    TecplotIO tecplot_io(*mesh,false);

    tecplot_io.write_equation_systems (filename,es);

  }

  //package up a local version of the solution
  if (stream_vectors != NULL) {
    std::vector<double> xfer;
    system.solution->localize(xfer);
    stream_vectors->clear();
    stream_vectors->resize(mesh->n_elem());
    libmesh_assert(stream_vectors->size() == (xfer.size()/3));
    for (unsigned int i=0; i<stream_vectors->size(); i++) {
      (*stream_vectors)[i](0) = xfer[3*i];
      (*stream_vectors)[i](1) = xfer[3*i + 1];
      (*stream_vectors)[i](2) = xfer[3*i + 2];
    }

  }

}


void create_stream_vectors_cbaero_method(_Mesh* mesh,
					 const std::vector<double>& fs_in,
					 std::vector<libMesh::Point>* stream_vectors,
					 const std::string& filename) {

  //create EquationSystems object
  EquationSystems es(*mesh);

  //add explicit system to hold the streamline solution
  es.add_system<ExplicitSystem>("Stream Vectors");

  //get a reference to the explicit system we just made
  ExplicitSystem& system = es.get_system<ExplicitSystem>(0);

  //add the variables for 3 components of streamline vector
  system.add_variable("nx",CONSTANT,XYZ);
  system.add_variable("ny",CONSTANT,XYZ);
  system.add_variable("nz",CONSTANT,XYZ);

  /*
    initialize the equation systems obect for use
    this is a trivial (in this instance) but necessary call to prepare the system for use
    It will initialize the solution to zero, but it doesn't matter since we will overwrite the solution
  */
  es.init();

  /*
    this method needs surface normals
    this call will skip calculation of element normal data if it has already been done and stored for this mesh object
   */
  mesh->populate_elem_normal();

  //initialize the freestream vector from user input and make it a unit vector
  Point fs(fs_in[0],fs_in[1],fs_in[2]);

  //objects that we do not want to re-allocate each time through the loop
  Point streamline;

  //loop over local nodes
  MeshBase::const_element_iterator       el     = mesh->active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_local_elements_end();

  for ( ; el != end_el; ++el) {

    const Elem* elem = *el;

    const Point& normal = mesh->get_elem_normal(elem->id());

    streamline = ((normal.cross(fs)).cross(normal)).unit();

    system.solution->set (elem->dof_number(0,0,0), streamline(0));
    system.solution->set (elem->dof_number(0,1,0), streamline(1));
    system.solution->set (elem->dof_number(0,2,0), streamline(2));

  }//end loop over nodes

  //close solution vector to prevent further modification and prepare if for parallel communication
  system.solution->close();

  //sync solution amongst processors
  system.update();

  //write tecplot file if user wants.  currently this is hardcoded for ascii since its intent is for debugging
  if (filename != "none" ) {

    TecplotIO tecplot_io(*mesh,false);

    tecplot_io.write_equation_systems (filename,es);

  }

  //package up a local version of the solution
  if (stream_vectors != NULL) {
    std::vector<double> xfer;
    system.solution->localize(xfer);
    stream_vectors->clear();
    stream_vectors->resize(mesh->n_elem());
    libmesh_assert(stream_vectors->size() == (xfer.size()/3));
    for (unsigned int i=0; i<stream_vectors->size(); i++) {
      (*stream_vectors)[i](0) = xfer[3*i];
      (*stream_vectors)[i](1) = xfer[3*i + 1];
      (*stream_vectors)[i](2) = xfer[3*i + 2];
    }

  }

}


void create(StreamlineContainer* streamline_container,
	    //const StagnationPointMethodEnum& stag_method,
	    const StreamVectorMethodEnum& stream_method,
	    const std::vector<double>& fs) {

  //streamline_container->create_all_streamlines(stag_method,stream_method,fs);
  streamline_container->create_all_streamlines(MOST_FORWARD_ELEM,stream_method,fs);

}

void write(StreamlineContainer* streamline_container, const std::string& filename) {
  streamline_container->write(filename);
}

void clear(StreamlineContainer* streamline_container) {
  streamline_container->clear_data();
}

StreamlineContainer* get_new_pointer(_Mesh* mesh) {
  return(new StreamlineContainer(*mesh));
}

void delete_pointer(StreamlineContainer* streamline_container) {
  delete streamline_container;
}
