#ifndef __equation_systems_h__
#define __equation_systems_h__

//libmesh includes
#include "libmesh/equation_systems.h"
#include "libmesh/explicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/tecplot_io.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/vtk_io.h"

//local includes
#include "mesh.h"
#include "steam_libmesh_common.h"

//function prototypes

/*
  BEGIN Functions that are wrapped by python EquationSystems class
*/

/*
struct ElemNodeStorage{
    enum Enum{ ELEMENTS = 0, NODES = 1 };
};

struct ModelType{
    enum Enum{ MULTIPLICATIVE = 0, ADDITIVE = 1, REPLACEMENT = 2 };
};

struct OutputType{
    enum Enum { TECPLOT_DAT = 0, TECPLOT_PLT = 1, EXODUS = 2, VTK_OUT = 3};
};
*/

libMesh::EquationSystems* get_new_pointer(_Mesh* mesh);

void delete_pointer(libMesh::EquationSystems* equation_systems);

void print_info(libMesh::EquationSystems* equation_systems);

void add_transient_explicit_system(libMesh::EquationSystems* equation_systems, const std::string& name);

void init(libMesh::EquationSystems* equation_systems);
void reinit(libMesh::EquationSystems* equation_systems);
unsigned int n_systems(libMesh::EquationSystems* equation_systems);


void delete_system(libMesh::EquationSystems* equation_systems, const std::string& name);

libMesh::MeshBase* get_mesh_ptr( libMesh::EquationSystems* es );
//_Mesh* get_mesh_ptr( libMesh::EquationSystems* es );

//libMesh::System* get_sys_by_num(libMesh::EquationSystems* equation_systems, 
//                              const unsigned int& num);

libMesh::ExplicitSystem* add_explicit_system(libMesh::EquationSystems* equation_systems, const std::string& name);
//libMesh::System* add_explicit_system(libMesh::EquationSystems* equation_systems, const std::string& name);

unsigned int add_system_variable(
      libMesh::ExplicitSystem* system, 
      const std::string& name,
      const ElemNodeStorage::Enum node_elem_flag
      );

void fill_system_variable(
      libMesh::ExplicitSystem* system, 
      const std::string& name,
      const std::vector<double>& val,
      const ElemNodeStorage::Enum node_elem_flag
      );

void node_to_element(
                     libMesh::ExplicitSystem* nsystem, 
                     libMesh::ExplicitSystem* esystem, 
                     const std::string& name
                     );

void write(libMesh::EquationSystems* es,
           const OutputType::Enum type,
           const std::string& filename);

/*
  END Functions that are wrapped by python EquationSystems class
*/

//void augment_greater_than_0( libMesh::ExplicitSystem* system);
void model_at_nodes(
                    libMesh::ExplicitSystem* system,
                    const unsigned int var_num,
                    const ModelType::Enum model_flag,
                    const std::vector<unsigned int>& node_ids,
                    const std::vector<double>& factors
                    );

void model_at_elems(
                    libMesh::ExplicitSystem* system,
                    const unsigned int var_num,
                    const ModelType::Enum model_flag,
                    const std::vector<unsigned int>& elem_ids,
                    const std::vector<double>& factors
                    );

//void get_data( libMesh::ExplicitSystem* system,
//               std::vector<std::string>& var_names,
//               std::vector<double>& data
//               );
#endif
