#include "equation_systems.h"

using namespace libMesh;

/*
  BEGIN Functions that are wrapped by python EquationSystems class
*/

EquationSystems* get_new_pointer(_Mesh* mesh) {
  return(new EquationSystems(*mesh));
}

void delete_pointer(EquationSystems* equation_systems) {
  //std::cerr << "Hello (EqSy)!" << std::endl;
  delete equation_systems;
  //std::cerr << "Bye   (EqSy)!" << std::endl;
}

void print_info(EquationSystems* equation_systems) {
  equation_systems->print_info();
}

void add_transient_explicit_system(EquationSystems* equation_systems, const std::string& name) {
  equation_systems->add_system<TransientExplicitSystem>(name);
}

void init(EquationSystems* equation_systems) {
  equation_systems->init();
}

void reinit(EquationSystems* equation_systems) {
  equation_systems->reinit();
}

unsigned int n_systems(EquationSystems* equation_systems) {
  return equation_systems->n_systems();
}

void delete_system(EquationSystems* equation_systems, const std::string& name) {
   equation_systems->delete_system(name);
}

//_Mesh* get_mesh_ptr( EquationSystems* es ){
//  return static_cast<_Mesh*> (&es->get_mesh() );
//}
MeshBase* get_mesh_ptr( EquationSystems* es ){
  return &es->get_mesh();
}

// Handle Solution data
ExplicitSystem* add_explicit_system(
            EquationSystems* equation_systems, 
            const std::string& name) {
  ExplicitSystem& system = equation_systems->add_system<ExplicitSystem>(name);
  ExplicitSystem *ptr = &(system);
  return ptr;
}

unsigned int add_system_variable( 
                          ExplicitSystem* system, 
                          const std::string& name,
                          const ElemNodeStorage::Enum node_elem_flag ) {

  // OK, add this variable
  switch( node_elem_flag ){

  case ElemNodeStorage::ELEMENTS:
     {
     const unsigned int num = system->add_variable(name,CONSTANT,XYZ);
     return(num);
     }
     break;
  case ElemNodeStorage::NODES:
     {
     const unsigned int num = system->add_variable(name,FIRST,LAGRANGE);
     return(num);
     }
     break;


  default:
     break;

  //return(num);
  }

}




void fill_system_variable(
                          ExplicitSystem* system, 
                          const std::string& name,
                          const std::vector<double>& val,
                          const ElemNodeStorage::Enum node_elem_flag
                          ) {
/*
""" Populate and equaton system with variable from python.

   Args:
      system(ExplicitSystem*): System pointer to add varible to.
      name(string)           : Name of variable to add.
      val(double vector)     : Node-based values to populate.
"""
*/

  // Get the mesh
  MeshBase& mesh = system->get_mesh();

  // Find the variable number
  const unsigned int num = system->variable_number(name);

  switch( node_elem_flag ){

  case ElemNodeStorage::ELEMENTS:
     {
     //loop over local elements
     MeshBase::const_element_iterator         el = mesh.local_elements_begin();
     const MeshBase::const_element_iterator end_el = mesh.local_elements_end();

        unsigned int i = 0;
        for ( ; el != end_el; ++el) {

          //dereference the iterator to get a pointer to the node
          const Elem* elem = *el;

         /*
            set the solution
            system.solution->set (node->dof_number(system number,
                           variable number, component number), solution value);
            component number = 0 (trivial) since we are not invoking 
                                           that logic in libmesh
          */
          system->solution->set(
                            elem->dof_number(system->number(),num,0), val[i]);

          i = i + 1;

        }
     }
     break;
  case ElemNodeStorage::NODES:
     {
        //loop over local nodes
        MeshBase::const_node_iterator           n = mesh.local_nodes_begin();
        const MeshBase::const_node_iterator end_n = mesh.local_nodes_end();


        unsigned int i = 0;
        for ( ; n != end_n; ++n) {

          //dereference the iterator to get a pointer to the node
          const Node* node = *n;

         /*
            set the solution
            system.solution->set (node->dof_number(system number,
                           variable number, component number), solution value);
            component number = 0 (trivial) since we are not invoking 
                                           that logic in libmesh
          */
          system->solution->set(
                            node->dof_number(system->number(),num,0), val[i]);

          i = i + 1;

        }

     }
     break;
   default:
      break;
   }
  // close solution vector to prevent further modification and prepare if 
  // for parallel communication
  system->solution->close();

  // sync solution amongst processors
  system->update();

  // dump some info
  system->print_info();
}




void node_to_element(
                     ExplicitSystem* nsystem, 
                     ExplicitSystem* esystem, 
                     const std::string& name
                     ) {
/*
""" Convert node-based data to element-based data.

   Args:
      system(ExplicitSystem*): Node    system pointer.
      system(ExplicitSystem*): Element system pointer.
      name(string)           : Name of element-based variable to make.

"""
*/

   // Get the mesh - Should be the same for each
   MeshBase& mesh = esystem->get_mesh();

   // Find the variable number
   unsigned int n_snum = nsystem->number();
   unsigned int e_snum = esystem->number();
   unsigned int n_vnum = nsystem->variable_number(name);
   unsigned int e_vnum = esystem->variable_number(name);
   //std::cerr << "n_snum,e_snum " << n_snum << " " << e_snum << std::endl;
   //std::cerr << "n_vnum,e_vnum " << n_vnum << " " << e_vnum << std::endl;

   NumericVector<Real>& nsoln = *nsystem->solution;
   NumericVector<Real>& esoln = *esystem->solution;

   _Mesh::const_element_iterator       el     = mesh.active_local_elements_begin();
   const _Mesh::const_element_iterator end_el = mesh.active_local_elements_end();
   for ( ; el != end_el; ++el) {
      const Elem* elem = *el;
      
      double       value  = 0.0;
      unsigned int number = elem->n_nodes();

/*
   set the solution
   system.solution->set (node->dof_number(system number,variable number, component number), solution value);
   component number = 0 (trivial) since we are not invoking that logic in libmesh
 */
      unsigned int edof   = elem->dof_number(e_snum,e_vnum,0);
      //std::cerr << "elem_stats " << e_snum << " " << e_snum << " " << edof << std::endl;
      //std::cerr <<elem->get_info() << std::endl;

      // std::cerr << "e_dof " << edof << " " << number << std::endl;

      //loop over nodes on element and contribute to the element
      for (unsigned int n=0; n<number; n++) {
         
         unsigned int ndof = elem->node_ptr(n)->dof_number(n_snum,n_vnum,0);
         //std::cerr << "node_stats " << n_snum << " " << n_snum << " " << ndof << std::endl;
         //std::cerr <<elem->node_ptr(n)->get_info() << std::endl;

         value  += nsoln(ndof);
         //std::cerr << "node " << n << " " <<    nsoln(ndof) << std::endl;
         //std::cerr << "node " << n << " " << value << std::endl;
      }//end loop over nodes

      //esoln.set(edof, value/static_cast<double>(number));
      esoln.set(edof, value/static_cast<double>(number));
      //std::cerr << "e_soln " << esoln(edof) << std::endl;

   }

  //close solution vector to prevent further modification and prepare if for parallel communication
  esystem->solution->close();

  //sync solution amongst processors
  esystem->update();

  //dump some info
  esystem->print_info();


   // Write everything out as a debug step
   //el     = mesh.active_local_elements_begin();
   //for ( ; el != end_el; ++el) {
   //   const Elem* elem = *el;
   //   
   //   unsigned int edof   = elem->dof_number(e_snum,e_vnum,0);
   //   std::cerr << "e_soln " << esoln(edof) << std::endl;

   //}

}

void write( EquationSystems* es,
           const OutputType::Enum type,
           const std::string& filename) {

   switch (type) {
   case OutputType::TECPLOT_DAT:
      {
      TecplotIO tecplot_io(es->get_mesh(),false);
      tecplot_io.write_equation_systems (filename,*es);
      }
      break;
   case OutputType::TECPLOT_PLT:
      {
      TecplotIO tecplot_io(es->get_mesh(),true);
      tecplot_io.write_equation_systems (filename,*es);
      }
      break;
   case OutputType::EXODUS:
      {
      ExodusII_IO exodusii_io(es->get_mesh(),false);
      exodusii_io.write_equation_systems (filename,*es);
      }
      break;
   case OutputType::VTK_OUT:
      {
      VTKIO vtk_io(es->get_mesh() );
      vtk_io.write_equation_systems (filename,*es);
      }
      break;
   default:
      break;
   }

}

/*
  END Functions that are wrapped by python EquationSystems class
*/


void model_at_nodes(
                    ExplicitSystem* system,
                    const unsigned int var_num,
                    const ModelType::Enum model_flag,
                    const std::vector<unsigned int>& node_ids,
                    const std::vector<double>& factors
                    ) {

  /*
    Take a pointer to an equation_systems instance and augment the specified
    solution quantity in the desired way.

    Args:
        system:     Pointer to equation_systems object
        var_num:    0-indexed variable on which to operate
        model_flag: Flag indicating how to use the supplied factors. Values:
                        1       Multiply by factors
                        2       Add factors
                        3       Replace with factors
        factors:    Values to use for augmentation
  */
  
  // Get the mesh and other constants
  MeshBase& mesh = system->get_mesh();
  NumericVector<Real>& soln = *system->solution;
  const unsigned int sys_num = system->number();

   /*
      When setting the solution:
      system.solution->set (node->dof_number(system number,variable number, 
                                            component number), solution value);
      component number = 0 (trivial) since we are not invoking 
                                     that logic in libmesh
    */

  unsigned int dof;

  // Accommodate different operations on the value of the variable of interest
  switch (model_flag) {
  case ModelType::MULTIPLICATIVE:
    {
      // Multiplicative model
      for( unsigned int i = 0; i < node_ids.size(); i++ ){
        
        // Allow this to work in parallel 
        // -- only operate on nodes assigned to the local processor
        if( mesh.node( node_ids[i] ).processor_id() == mesh.comm().rank() ){

            dof = mesh.node( node_ids[i] ).dof_number(sys_num,var_num,0);
            soln.set( dof, soln(dof) * factors[i] );
        }
      }
    }
    break;
  case ModelType::ADDITIVE:
    {
      // Additive model
      for( unsigned int i = 0; i < node_ids.size(); i++ ){
        
        // Allow this to work in parallel 
        // -- only operate on nodes assigned to the local processor
        if( mesh.node( node_ids[i] ).processor_id() == mesh.comm().rank() ){

            dof = mesh.node( node_ids[i] ).dof_number(sys_num,var_num,0);
            soln.set( dof, soln(dof) + factors[i] );
        }
      }
    }
    break;
  case ModelType::REPLACEMENT:
    {
      // Replacement model
      for( unsigned int i = 0; i < node_ids.size(); i++ ){
        
        // Allow this to work in parallel 
        // -- only operate on nodes assigned to the local processor
        if( mesh.node( node_ids[i] ).processor_id() == mesh.comm().rank() ){

            dof = mesh.node( node_ids[i] ).dof_number(sys_num,var_num,0);
            soln.set( dof, factors[i] );
        }
      }
    }
    break;
  default:
    break;
  }

  //close solution vector to prevent further modification and 
  //prepare it for parallel communication
  system->solution->close();

  //sync solution amongst processors
  system->update();

  //dump some info
  system->print_info();
}

void model_at_elems(
                    ExplicitSystem* system,
                    const unsigned int var_num,
                    const ModelType::Enum model_flag,
                    const std::vector<unsigned int>& elem_ids,
                    const std::vector<double>& factors
                    ) {

  /*
    Take a pointer to an equation_systems instance and augment the specified
    solution quantity in the desired way.

    Args:
        system:     Pointer to equation_systems object
        var_num:    0-indexed variable on which to operate
        model_flag: Flag indicating how to use the supplied factors. Values:
                        0       Multiply by factors
                        1       Add factors
                        2       Replace with factors
        factors:    Values to use for augmentation
  */
  
  // Get the mesh
  MeshBase& mesh = system->get_mesh();
  NumericVector<Real>& soln = *system->solution;
  const unsigned int sys_num = system->number();

   /*
      When setting the solution:
      system.solution->set (elem->dof_number(system number,variable number, 
                                            component number), solution value);
      component number = 0 (trivial) since we are not invoking 
                                     that logic in libmesh
    */

  unsigned int dof;

  // Accommodate different operations on the value of the variable of interest
  switch (model_flag) {
  case ModelType::MULTIPLICATIVE:
    {
      // Multiplicative model
      for( unsigned int i = 0; i < elem_ids.size(); i++ ){
        const Elem& elem = mesh.elem_ref( elem_ids[i] );
        
        // Allow this to work in parallel 
        // -- only operate on elems assigned to the local processor
        if( elem.processor_id() == mesh.comm().rank() ){

            dof = elem.dof_number(sys_num,var_num,0);
            soln.set( dof, soln(dof) * factors[i] );
        }
      }
    }
    break;
  case ModelType::ADDITIVE:
    {
      // Additive model
      for( unsigned int i = 0; i < elem_ids.size(); i++ ){
        const Elem& elem = mesh.elem_ref( elem_ids[i] );
        
        // Allow this to work in parallel 
        // -- only operate on elems assigned to the local processor
        if( elem.processor_id() == mesh.comm().rank() ){

            dof = elem.dof_number(sys_num,var_num,0);
            soln.set( dof, soln(dof) + factors[i] );
        }
      }
    }
    break;
  case ModelType::REPLACEMENT:
    {
      // Replacement model
      for( unsigned int i = 0; i < elem_ids.size(); i++ ){
        const Elem& elem = mesh.elem_ref( elem_ids[i] );
        
        // Allow this to work in parallel 
        // -- only operate on elems assigned to the local processor
        if( elem.processor_id() == mesh.comm().rank() ){

            dof = elem.dof_number(sys_num,var_num,0);
            soln.set( dof, factors[i] );
        }
      }
    }
    break;
  default:
    break;
  }

  //close solution vector to prevent further modification and 
  //prepare it for parallel communication
  system->solution->close();

  //sync solution amongst processors
  system->update();

  //dump some info
  system->print_info();
}
