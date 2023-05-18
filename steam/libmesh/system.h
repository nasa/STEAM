#ifndef __system_h__
#define __system_h__

//libmesh includes
#include "libmesh/equation_systems.h"
#include "libmesh/explicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/numeric_vector.h"

//local includes
#include "mesh.h"
#include "steam_libmesh_common.h"

//function prototypes

libMesh::EquationSystems* get_parent( libMesh::System* system );

std::string get_name( libMesh::System* system );
unsigned int get_number( libMesh::System* system );

void delete_pointer(libMesh::System* system);

void print_info(libMesh::System* system);

void get_data_by_num( 
                    libMesh::EquationSystems* es,
                    unsigned int sys_num,
                    std::vector<std::string>& var_names,
                    std::vector<double>& data
                    );

void get_data_by_name(
                    libMesh::EquationSystems* es,
                    const std::string sys_name,
                    std::vector<std::string>& var_names,
                    std::vector<double>& data,
                    const ElemNodeStorage::Enum node_elem_flag,
                    const SolnSplit::Enum soln_split_flag,
                    const unsigned int output_rank = 0
                    );

#endif
