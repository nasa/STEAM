#include "system.h"

using namespace libMesh;

/*
  BEGIN Functions that are wrapped by python EquationSystems class
*/

//System* get_new_pointer(EquationSystems* equation_systems,
//                        const std::string& name, 
//                        const unsigned int number) {
//  return(new System(*equation_systems, name, number ));
//}

////////////////////////////////////////////////

EquationSystems* get_parent( System* system ){
  return &(system->get_equation_systems());
}

////////////////////////////////////////////////

/*
// ExplicitSystem version:
std::string get_name( ExplicitSystem* system ){
  return system->name();
}

unsigned int get_number( ExplicitSystem* system ){
  return system->number();
}
*/

////////////////////////////////////////////////

std::string get_name( System* system ){
  return system->name();
}

////////////////////////////////////////////////

unsigned int get_number( System* system ){
  return system->number();
}

////////////////////////////////////////////////

void delete_pointer(System* system) {
  //std::cerr << "Hello (EqSy)!" << std::endl;
  delete system;
  //std::cerr << "Bye   (EqSy)!" << std::endl;
}

////////////////////////////////////////////////

void print_info(System* system) {
  system->print_info();
}

////////////////////////////////////////////////

void get_data_by_num( 
                    EquationSystems* es,
                    unsigned int sys_num,
                    std::vector<std::string>& var_names,
                    std::vector<double>& data
                    ) {
/*
"""Get data for translation of a libMesh solution into Python.

Retrieve the system by number.
"""
 */
  // Get the system of interest
  System* system = &(es->get_system(sys_num));

  // Get the solution
  NumericVector<Real>& soln = *system->solution;

  // Get the list of all variable names
  const unsigned int n_vars = system->n_vars();
  var_names.resize( n_vars );
  for( unsigned int i = 0; i < n_vars; i++ ){
    var_names[i] = system->variable_name( i );
  }

  // Get the data
  data.resize( soln.size() );
  system->solution->localize( data );
}

////////////////////////////////////////////////

void get_data_by_name(
                    EquationSystems* es,
                    const std::string sys_name,
                    std::vector<std::string>& var_names,
                    std::vector<double>& data,
                    const ElemNodeStorage::Enum node_elem_flag,
                    const SolnSplit::Enum soln_split_flag,
                    const unsigned int output_rank
                    ) {
/*
""" Get data for translation of a libMesh solution into Python.

   Args:
      system(ExplicitSystem*): System pointer to add varible to.
      name(string)           : Name of variable to add.
      val(double vector)     : Node-based values to populate.
"""
*/
  // Get the system of interest
  System* system = &(es->get_system(sys_name));
  const unsigned int sys_num = system->number();
  
  // Initialize array of indices and number of variables
  std::vector<numeric_index_type> indices;
  const unsigned int n_vars = system->n_vars();

  // Get all variable names
  var_names.resize( n_vars );
  for ( unsigned int j = 0; j != n_vars; ++j ) {
    var_names[j] = system->variable_name(j);
  }

  // Get the mesh
  MeshBase& mesh = system->get_mesh();

  switch( node_elem_flag ){

  case ElemNodeStorage::ELEMENTS:
     {
     // Resize indices vector
     indices.resize( n_vars * mesh.n_local_elem() );

     // Loop over elements using the element iterator
     MeshBase::const_element_iterator           n = mesh.local_elements_begin();
     const MeshBase::const_element_iterator end_n = mesh.local_elements_end();
     for ( unsigned int i = 0; n != end_n; ++n ) {

       //dereference the iterator to get a pointer to the element
       const Elem* element = *n;

       // Loop over variables
       for ( unsigned int j = 0; j != n_vars; ++j ) {
         /*
           get the index, which comes from dof_number
           element->dof_number(system number,variable number, component number)
           variable number = j since we are looping over all variables
           component number = 0 (trivial) since we are not invoking 
                                          that logic in libmesh
         */
         indices[i*n_vars + j] = element->dof_number( sys_num, j, 0 );
       }
         
       i = i + 1;   // Increment i because we're looping over elements, not i

     }

     // Actually retrieve the data
     switch( soln_split_flag ){

     case SolnSplit::LOCAL_SOLN:
        {
        // Return just the data for local elements
        system->solution->get( indices, data );
        }
     case SolnSplit::ENTIRE_SOLN:
        {
        // Return the entire solution, but only on the output_rank
        system->solution->localize_to_one( data, output_rank );
        }
   }


  }
     break;
  case ElemNodeStorage::NODES:
     {
     // Resize indices vector
     indices.resize( n_vars * mesh.n_local_nodes() );

     // Loop over nodes using the node iterator
     MeshBase::const_node_iterator           n = mesh.local_nodes_begin();
     const MeshBase::const_node_iterator end_n = mesh.local_nodes_end();
     for ( unsigned int i = 0; n != end_n; ++n ) {

       //dereference the iterator to get a pointer to the node
       const Node* node = *n;

       // Loop over variables
       for ( unsigned int j = 0; j != n_vars; ++j ) {
         /*
           get the index, which comes from dof_number
           node->dof_number(system number,variable number, component number)
           variable number = j since we are looping over all variables
           component number = 0 (trivial) since we are not invoking 
                                          that logic in libmesh
         */
         indices[i*n_vars + j] = node->dof_number( system->number(), j, 0 );
       }
         
       i = i + 1;   // Increment i because we're looping over nodes, not i

     }

   // Actually retrieve the data
   system->solution->get( indices, data );

   }
     break;
   default:
      break;
   }

}
