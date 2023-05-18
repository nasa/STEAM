#include "modified_newtonian.h"

using namespace libMesh;

void get_global_solution(_Mesh* mesh,
			 const double& V_mps,
			 const double& density_kgpm3,
			 const double& T_K,
			 const double& gamma,
			 const double& R_JpkgK,
			 const std::vector<double>& fs_in,
			 const std::string& filename) {

  //create EquationSystems object
  EquationSystems es(*mesh);

  //add explicit system to hold the pressure silution
  es.add_system<ExplicitSystem>("Modified Newtonian");

  //get a reference to the explicit system we just made
  ExplicitSystem& system = es.get_system<ExplicitSystem>(0);

  //add the pressure variable to the system
  system.add_variable("P_Pa",CONSTANT,XYZ);

  /*
    initialize the equation systems obect for use
    this is a trivial (in this instance) but necessary call to prepare the system for use
    It will initialize the pressure solution to zero, but it doesn't matter since we will overwrite the solution
  */
  es.init();

  /*
    this method needs surface normals
    this call will skip calculation of elem normal data if it has already been done and stored for this mesh object
   */
  mesh->populate_elem_normal();

  //initialize the freestream vector from user input and make it a unit vector
  Point fs(fs_in[0],fs_in[1],fs_in[2]);
  fs = fs.unit();

  double P;      //surface pressure in Pa.  This is the solution we are seeking
  double sin;    //sine of flow deflection angle

  const double M2     = V_mps*V_mps/(gamma*R_JpkgK*T_K);//freestream Mach^2. Derived from input state data
  const double q_inf  = 0.5*density_kgpm3*V_mps*V_mps;  //freestream dynamic pressure in Pa. Derived from input state data
  const double P_inf  = density_kgpm3*R_JpkgK*T_K;      //freestream pressure in Pa. Derived from input state data
  const double Cp_max = 2.0/(gamma*M2)*(std::pow((std::pow((gamma+1),2)*M2)/(4.0*gamma*M2-2.0*(gamma-1.0)),(gamma/(gamma-1.0)))*((1-gamma+2.0*gamma*M2)/(gamma+1.0))-1); //maximum pressure coefficient

  //loop over local nodes
  _Mesh::const_element_iterator           el = mesh->local_elements_begin();
  const _Mesh::const_element_iterator end_el = mesh->local_elements_end();

  for ( ; el != end_el; ++el) {

    //dereference the iterator to get a pointer to the element
    const Elem* elem = *el;

    //get the sine of the flow deflection angle
    sin = fs*mesh->get_elem_normal(elem->id());

    //check to make sure this node has a view of the flow
    if (sin > 0.0) {

      //no view, set pressure to zero
      P = 0.0;

    } else {

      //this node may have a view, we will assume it does for now
      P = P_inf + q_inf*Cp_max*sin*sin;

    }

    /*
      set the solution
      system.solution->set (node->dof_number(system number,variable number, component number), solution value);
      system number = 0 since we have only added one system
      variable number = 0 since we have only added on variable
      component number = 0 (trivial) since we are not invoking that logic in libmesh
    */
    system.solution->set (elem->dof_number(0,0,0), P);

  }

  //close solution vector to prevent further modification and prepare if for parallel communication
  system.solution->close();

  //sync solution amongst processors
  system.update();

  //write tecplot file if user wants.  currently this is hardcoded for ascii since its intent is for debugging
  if (filename != "none" ) {

    TecplotIO tecplot_io(*mesh,false);

    tecplot_io.write_equation_systems (filename,es);

  }

}
