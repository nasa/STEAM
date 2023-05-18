## We do a check here to make sure that SWIG is built
## Otherwise, it totally kills the Sphinx build
import logging
logger = logging.getLogger(__name__)
logger.debug("Top of file")

import os
import errno
from warnings import warn
dir_root = os.path.dirname(os.path.abspath(__file__))
if os.path.isfile(dir_root+'/init.py'):
    from . import init
    from . import mesh
    from . import equation_systems
    from . import system
    from . import mpi
import sys
import atexit
from abc import ABCMeta, abstractmethod

################################################################################

def store_at_str2enum( store_at_flag ):
    """ Convert a string flag indicating data storage at either nodes or 
    elements into an integer corresponding to the appropriate value from
    ElemNodeStorage::Enum (defined in steam_libmesh_common.h).

    Args:
        store_at_flag(:obj:`string`): String flag matching either 'ELEMENTS' or 'NODES' (case insensitive)

    Returns:
        :obj:`int`: Integer consistent with ElemNodeStorage::Enum defined in steam_libmesh_common.h
    """

    if store_at_flag.upper().find( 'ELEM' ) == 0:
        return equation_systems.ElemNodeStorage.ELEMENTS
    elif store_at_flag.upper().find( 'NODE' ) == 0:
        return equation_systems.ElemNodeStorage.NODES
    else:
        raise ValueError( 'store_at_flag must match either "ELEM" or "NODE"' )

################################################################################

def soln_split_str2enum( soln_split_flag ):
    """ Convert a string flag indicating whether the solution transferred
    between Python and LibMesh is just the component local to the rank or the
    entire solution.  This function will return the appropriate value from
    SolnSplit::Enum (defined in steam_libmesh_common.h).

    Args:
        soln_split_flag(:obj:`string`): String flag matching either 'LOCAL_SOLN' or 'ENTIRE_SOLN' (case insensitive)

    Returns:
        :obj:`int`: Integer consistent with SolnSplit::Enum defined in steam_libmesh_common.h
    """

    if soln_split_flag.upper().startswith( 'LOCAL' ):
        return system.SolnSplit.LOCAL_SOLN
    elif soln_split_flag.upper().startswith( 'ENTIRE' ):
        return system.SolnSplit.ENTIRE_SOLN
    else:
        raise ValueError( 
                    'soln_split_flag must match either "LOCAL" or "ENITRE"' )

################################################################################

class LibObj( metaclass = ABCMeta ):
    """
    SuperClass for all Python implementations of libMesh objects.  This
    basically means all other classes defined in this file (except for Init).

    NOTE: This is an abstract superclass and as such instances of LibObj 
          cannot be directly instantiated.  Instead, other classes can inherit
          from this one and then use the functionality provided here.
    """

    pointer_dict = {}

    def __init__(self, *args, **kwargs):
        """ Constructor for the LibObj superclass.  This method should not be
        replaced or augmented by the derived classes.  Instead, it will call
        the derived classes' self.derived_init method and then perform the
        necessary bookkeeping of pointers.
        """

        ### Call the derived class init method
        self.derived_init(*args, **kwargs)

        ### Add the pointer to the storage dictionary if it doesn't exist
        self._ptr_id = id(self.p)
        if self._ptr_id not in LibObj.pointer_dict:
            LibObj.pointer_dict[self._ptr_id] = 0

        ### Increment the pointer count
        LibObj.pointer_dict[self._ptr_id] += 1

    def __del__(self):
        """ Destructor for the LibObj superclass.  This method should not be
        replaced or augmented by the derived classes.  Instead, it will call
        the derived classes' self.delete_pointer methods when appropriate.
        """

        ### Decrease the pointer count
        LibObj.pointer_dict[self._ptr_id] -= 1
        
        if LibObj.pointer_dict[self._ptr_id] <0:
            raise ValueError( 'THIS IS A PROBLEM' )

        if LibObj.pointer_dict[self._ptr_id] == 0:
            ### Delete the pointer from the dictionary and delete the pointer
            del LibObj.pointer_dict[self._ptr_id]
            self.delete_pointer()

    ### delete_pointer and derived_init must be overwritten by derived classes
    @abstractmethod
    def delete_pointer(self):
        pass

    @abstractmethod
    def derived_init(self):
        pass

################################################################################

class Init( LibObj ):
    """
    Class serves as python interface to libMesh::LibMeshInit object.
    It is intended to be a smart pointer that will deallocate memory newed 
    up in C++ when this python object goes out of scope.

    Note that only one Init object should be created per code.  When creating
    an Init object, STEAM checks to see if one exists already and will point
    to a preexisting one if possible.
    """

    def derived_init (self, mpi_comm=None):
        """Constructor calls SWIG function to create libMesh::LibMeshInit object
        """
        ### This init object will store a pointer to itself as 
        ### steam.libmesh._lmp because all other LibMesh-derived objects 
        ### need access to it. 
        global _lmp
        try:
            _lmp
        except NameError:
            ### We want to end up here -- if not then another LibMeshInit
            ### object already exists, which isn't allowed
            _lmp = self

            if mpi_comm:
                self.p = init.get_new_pointer(sys.argv, mpi_comm)
            else:
                self.p = init.get_new_pointer(sys.argv)

        else:
            #raise RuntimeError( 'LibMesh initialized more than once.' )
            warn( RuntimeWarning( 'LibMesh initialized more than once.  ' + 
                                  'Returning first LibMeshInit object.' ) )
            self.p = _lmp.p

    def delete_pointer (self):
        """Destructor calls SWIG function to delete libMesh::LibMeshInit object
        """
        init.delete_pointer(self.p)

    def __call__(self):
        """ overloaded () to return pointer """
        return self.p

## Now, make a pointer that can be used by anything that imported this package.
#if os.path.isfile(dir_root+'/init.py'):
#    try:
#        _lmp
#    except NameError:
#        _lmp = Init()

# Add a hook to delete the libmesh pointer before we leave python
def clean_lmp():
    global _lmp
    try:
        _lmp
    except NameError:
        pass
    else:
        del _lmp

atexit.register(clean_lmp)

#########################
#def init_lib( mpi_comm=None ):
#    """ Allow a user to reinitialize -- create a new LibMeshInit object with
#    a user-supplied MPI communicator.
#
#        Args:
#            mpi_comm (:obj:`mpi4py.MPI.Intracomm`): MPI communicator that will become the new communicator for LibMesh.  If mpi_comm is `None`, MPI_COMM_WORLD will be used by default.
#    """
#
#    global _lmp
#
#    ### Delete the old LibMeshInit object
#    clean_lmp()
#
#    try:
#        _lmp
#    except NameError:
#        ### We want to end up here -- if not then we haven't properly cleaned
#        ###   up the old _lmp.
#        _lmp = Init( mpi_comm=mpi_comm )
#    else:
#        raise RuntimeError( 'Improperly cleaned up old _lmp' )


################################################################################

class Mesh( LibObj ):
    """
    Class serves as python interface to |obj_Mesh| object which is derived 
    from libMesh::Mesh object.  It is intended to be a smart pointer that 
    will deallocate memory newed up in C++ when this python object goes out 
    of scope.
    """

    def derived_init (self,filename=None,pointer=None):
        """
        Constructor for Mesh object

        Args:
            filename(:obj:`string`,optional,input): string containing mesh filename
            pointer(pointer to C++ |obj_Mesh| object,input): pointer to C++ |obj_Mesh| object

        """
        # For garbage collection to work in the right order, save dependant
        try:
            self.parent = _lmp
        except NameError:
            raise RuntimeError( 'LibMesh not yet initialized.  ' +
                  'Must call steam.libmesh.Init() before using any other ' +
                  'LibMesh-derived object.' )

        # throw an error if mis-using the constructor.  This is an example of why Python sucks.
        if ((filename is not None) and (pointer is not None)):
            print('ERROR: Cannot construct Mesh object given a valid pointer and filename.')
            print('       One of the two arguments must be invalid (equal to None).')
            raise RuntimeError('ERROR: Cannot construct Mesh object given a valid pointer and filename.')
            return

        # If we're passing in a pointer, then take it.
        if (pointer is not None):
            self.p = pointer
            return

        # get a pointer to a fresh _Mesh object
        self.p = mesh.get_new_pointer(_lmp())

        # if provided a filename, initialize the mesh from the file
        if (filename is not None):
            if not os.path.isfile( filename ):
                raise FileNotFoundError(
                        errno.ENOENT, os.strerror(errno.ENOENT), filename)
            mesh.read(self.p,filename)

    @classmethod
    def from_python(cls, py_mesh):
        """
        Create mesh object from passed-in database object.

        Args:
            py_mesh(:obj:`steam.mesh.Mesh()`): Python mesh object.
        Returns:
            :obj:`Mesh()`: libMesh mesh object.
        """
        # Package grid elements for libmesh
        v1  = (py_mesh.xyz_pt['X']).tolist()
        v2  = (py_mesh.xyz_pt['Y']).tolist()
        v3  = (py_mesh.xyz_pt['Z']).tolist()
        n1  = (py_mesh.conn[0]    ).tolist()
        n2  = (py_mesh.conn[1]    ).tolist()
        n3  = (py_mesh.conn[2]    ).tolist()
        comp= (py_mesh.comp[0]    ).tolist()

        # Populate Libmesh Data Types
        x     = mesh.VectorDouble ()
        y     = mesh.VectorDouble ()
        z     = mesh.VectorDouble ()
        elems = mesh.VectorInt ()

        for i in range (len (v1)):
            x.append (v1[i])
            y.append (v2[i])
            z.append (v3[i])

        for i in range (len (n1)):
            elems.append (n1[i]  )
            elems.append (n2[i]  )
            elems.append (n3[i]  )
            elems.append (comp[i])

        new_mesh = cls()

        # Send to Libmesh to populate a mesh object
        mesh.make(new_mesh.p, x, y, z, elems)

        return new_mesh

    def delete_pointer (self):
        """ Destructor. Cleans up memory allocated in C++ code.

        """

        # call c++ function to delete _Mesh object that out pointer owns
        mesh.delete_pointer(self.p)

    def __call__(self):
        """ Overloaded () to return pointer to C++ |obj_Mesh| object

        """
        return self.p

    def __str__(self):
        """ Overloaded print() to wrap libMesh::Mesh::print_info()

        Returns:
            (:obj:`string`): Returns " " just for neatness

        """
        mesh.print_info(self.p)
        return " "

    def print_info(self):
        """ Wrapper for libMesh::Mesh::print_info()

        """
        mesh.print_info(self.p)

    def write(self,filename):
        """ Wrapper for libMesh::Mesh::write()

        Args:
            filename(:obj:`string`,input):  Filename for writing mesh file.  The extension dictates the file format.  See libMesh documentation for details.

        """
        mesh.write(self.p,filename)

    def read(self,filename):
        """ Wrapper for libMesh::Mesh::read()

        Args:
            filename(:obj:`string`,input):  Filename for reading mesh file.  The extension dictates the file format.  See libMesh documentation for details.

        """
        mesh.read(self.p,filename)

    def make(self,x,y,z,elems):
        """
        Populates :obj:`steam.mesh.Mesh()` object given nodal coordinates and triangle element definition.

        Args:
            x(:obj:`list`,input): list of floats containing x-coordinate for each node
            y(:obj:`list`,input): list of floats containing y-coordinate for each node
            z(:obj:`list`,input): list of floats containing z-coordinate for each node
            elems(:obj:`list`,input): list of integers containing 3 node numbers for right-handed triangle element plus a 4th integer for component id.  List will be 4*(Number of Elements) long.

        """
        mesh.make(self.p,x,y,z,elems)

    def get_data(self,x,y,z,elems):
        """
        Retrieves :obj:`steam.mesh.Mesh()` object given nodal coordinates and triangle element definition.

        Args:
            x(:obj:`list`,input): x-coordinate for each node
            y(:obj:`list`,input): y-coordinate for each node
            z(:obj:`list`,input): z-coordinate for each node
            elems(:obj:`list`,input): integer list containing 3 node numbers for right-handed triangle element plus a 4th integer for component id.  List will be 4*(Number of Elements) long.

        """
        mesh.get_data(self.p,x,y,z,elems)

    def n_elem(self):
        """
        Wrapper for libMesh::Mesh::n_elem()

        Returns:
            :obj:`int`: number of elements in mesh

        """
        return mesh.n_elem(self.p)

    def n_nodes(self):
        """
        Wrapper for libMesh::Mesh::n_nodes()

        Returns:
            :obj:`int`: number of nodes in mesh

        """
        return mesh.n_nodes(self.p)

    def populate_node_normal(self):
        """
        Calculates and saves an outward unit normal at each node by simply averaging the unit normal components for each element that contains the node.

        """
        mesh.populate_node_normal(self.p)

    def write_node_normals(self,filename):
        """
        Writes node unit normals to tecplot ascii file

        Args:
            filename(:obj:`string`,input): filename to write node unit normals to

        """
        mesh.write_node_normals(self.p,filename)

    def populate_elem_normal(self):
        """
        Calculates and saves an outward unit normal at each element.

        """
        mesh.populate_elem_normal(self.p)

    def write_elem_normals(self,filename):
        """
        Writes element unit normals to tecplot ascii file

        Args:
            filename(:obj:`string`,input): filename to write element unit normals to

        """
        mesh.write_elem_normals(self.p,filename)

    def populate_elem_centroid(self):
        """
        Calculates and saves centroid for each element

        """
        mesh.populate_elem_centroid(self.p)

    def write_elem_centroids(self,filename):
        """
        Writes element centroids to tecplot ascii file

        Args:
            filename(:obj:`string`,input): filename to write element centroids to

        """
        mesh.write_elem_centroids(self.p,filename)

    def write_domain_decomposition(self,filename):
        """
        Writes parallel domain decomposition data to tecplot ascii file

        Args:
            filename(:obj:`string`,input): filename to write parallel domain decomposition data to

        """
        mesh.write_domain_decomposition(self.p,filename)

    def uniformly_refine(self,num_times=0):
        """
        Uniformly refines mesh specified number of times

        Args:
            num_times(:obj:`int`,input): number of times to uniformly refine mesh

        """
        mesh.uniformly_refine(self.p,num_times)

    def project(self,target):
        """ Project this mesh onto a target mesh and return.

        Args:
            target(:obj:`Mesh()`): libmesh mesh object to project onto.
        Returns:
            :obj:`Mesh()`: Projection onto target.
        """

        pointer = mesh.project_grid(self.p  ,target.p)
        output  = Mesh(pointer=pointer)
        return output

    def to_python(self):
        """ Return libmesh object data to python dictionaries.

        Returns:
            :obj:`dict`: Node X, Y, Z
            :obj:`dict`: Element connectivity
            :obj:`dict`: Element components
        """
        x =     mesh.VectorDouble ()
        y =     mesh.VectorDouble ()
        z =     mesh.VectorDouble ()
        elems = mesh.VectorInt ()

        self.get_data(x, y, z, elems)

        nodes = {'X' : [], 'Y' : [], 'Z' : []}
        ele   = {}
        comp  = {}
        for cn in range(len(x)):
            nodes['X'].append(x[cn])
            nodes['Y'].append(y[cn])
            nodes['Z'].append(z[cn])
        for ce in range (int(len(elems)/4)):
            ele[ce]  = []
            for i in range (3):
                ele[ce].append(elems[ce * 4 + i])
            comp[ce] = elems[ce * 4 + 3]

        return nodes, ele, comp

################################################################################

def read_mesh_and_get_data(filename,x,y,z,elems):
    """ Interfaces with libMesh to read mesh and return mesh data

    Args:
        filename(:obj:`string`,input): mesh filename to be read
        x(:obj:`list`,output): list of floats containing x-coordinate for each node
        y(:obj:`list`,output): list of floats containing y-coordinate for each node
        z(:obj:`list`,output): list of floats containing z-coordinate for each node
        elems(:obj:`list`,output): list of integers containing 3 node numbers for right-handed triangle element plus a 4th integer for component id.  List will be 4*(Number of Elements) long.
    """
    mesh.read_mesh_and_get_data(_lmp(),filename,x,y,z,elems)

################################################################################

class EquationSystems( LibObj ):
    """
    Class serves as python interface to libMesh::EquatonsSystems object
    It is intended to be a smart pointer that will deallocate memory newed 
    up in C++ when this python object goes out of scope.
    """

    def derived_init (self,in_mesh):
        """
        Constructor

        Args:
            in_mesh(:obj:`steam.libmesh.Mesh()`): Python mesh object which wraps |obj_Mesh| object, which is derived from libMesh::Mesh object.

        """

        # For garbage collection to work in the right order, save dependant
        self.parent = in_mesh
        self.p = equation_systems.get_new_pointer(in_mesh())
        
    def delete_pointer (self):
        """
        Destructor. Cleans up memory allocated in C++ code.

        """
        equation_systems.delete_pointer(self.p)

    def __call__(self):
        """
        Overloaded () to return pointer to C++ libMesh::EquationSystems object

        """
        return self.p

    def __str__(self):
        """
        Overloaded print() to wrap libMesh::EquationSystems::print_info()

        Returns:
            (:obj:`string`): Returns " " just for neatness

        """
        equation_systems.print_info(self.p)
        return " "

    def print_info(self):
        """
        Wrapper for libMesh::EquationSystems::print_info()

        """
        equation_systems.print_info(self.p)

    def add_transient_explicit_system(self,name):
        """
        Wrapper for libMesh::EquationSystems::add_system<TransientExplicitSystem>()

        Args:
            name(:obj:`string`): system name

        """
        equation_systems.add_transient_explicit_system(self.p,name)

    def init(self):
        """ Wrapper for libMesh::EquationSystems::init()

        """
        equation_systems.init(self.p)

    def reinit(self):
        """ Wrapper for libMesh::EquationSystems::reinit()

        """
        equation_systems.reinit(self.p)


    def n_systems(self):
        """ Wrapper for libMesh::EquationSystems::n_systems()

        Returns:
            :obj:`int`: Number of systems stored in an EquationSystems object.
        """
        equation_systems.n_systems(self.p)

    def get_mesh_ptr(self):
        """ Wrapper for libMesh::EquationSystems::get_mesh_ptr()

        Returns:
            :obj:`steam.libmesh.Mesh()`: Pointer to libMesh mesh object associated with the EquationSystems.
        """
        #mesh_ptr = equation_systems.get_mesh_ptr( self.p )
        #print( 'type(mesh_ptr):     ', type(mesh_ptr) )
        #return mesh_ptr
#        return equation_systems.BaseToMesh( 
#                                    equation_systems.get_mesh_ptr( self.p ) )
        return equation_systems.get_mesh_ptr( self.p )

    def write(self,type='PLT',filename="test.plt"):
        """
        Writes solution data associated with this EquationSystems object. Wraps libMesh::EquationSystems::write()

        Args:
            type(:obj:`string`): 'PLT' for tecplot, 'E' for ExodusII
            filename(:obj:`string`): filename to be written

        """
        if (type == "DAT"):
            itype = 0
        elif (type == "PLT"):
            itype = 1
        elif (type == "E"):
            itype = 2
        elif (type == "VTK"):
            itype = 3

        equation_systems.write(self.p,itype,filename)

    def from_python(self, py_soln, soln_name = 'nSoln'):
        """ Take all of the variables from the py_soln into libMesh.

        Args:
            py_soln(:obj:`steam.solution.Solution()`): Python solution object.
            soln_name(:obj:`string`): Name of the solution to be stored

        Returns:
            :obj:`string`: Name of the solution stored in libMesh

        """

        ### Do we need to save solnp and py_soln as attributes of self?
        ### I'm not sure that we do.  -- CJE
        self.solnp   = equation_systems.add_explicit_system(self.p,soln_name)
        self.py_soln = py_soln

        self.vars = self.py_soln.data.columns.tolist()


        ### Add all variables to the system
        self.var2num = {}
        self.var2num = { name : equation_systems.add_system_variable( 
                self.solnp, name, store_at_str2enum( self.py_soln.store_at ) )
                for name in self.vars}

        self.init()

        for name in self.vars:
            varray = equation_systems.VectorDouble ()
            v      = (self.py_soln.data[name]).tolist()
            for i in (v):
                varray.append(i)

            equation_systems.fill_system_variable(
                    self.solnp, name, varray, 
                    store_at_str2enum( self.py_soln.store_at ) )

        ### Return the name of the system that corresponds to this solution
        return system.get_name( system.ExpToSys(self.solnp) )
        #return system.get_name( self.solnp )

    def sys_to_python(self, sys_name = None, sys_num = None, 
                      local_or_entire='LOCAL_SOLN', output_rank=0 ):
        """ Return libmesh System object data to python dictionaries.

        Args:
            sys_name(:obj:`string`): Name of the system to be recalled from libMesh.
            sys_num(:obj:`string`): Number of the system to be recalled from libMesh.
            local_or_entire(:obj:`string`): Flag indicating whether to return the solution local to the rank or the entire solution.  Acceptable values: "LOCAL_SOLN" or "ENTIRE_SOLN"
            output_rank(:obj:`int`): The rank that will be returning the entire solution if local_or_entire=="ENTIRE_SOLN"
        Returns:
            :obj:`dict`: All variable values stored with variable names as keys.
        """

        ### Initialize the data storage vectors
        var_names = system.VectorStr()
        var_vals  = system.VectorDouble()

        ### Get a pointer to the system of interest
        if sys_name is not None:
            system.get_data_by_name( self.p, sys_name, var_names, var_vals, 
                                     store_at_str2enum(self.py_soln.store_at),
                                     soln_split_str2enum(local_or_entire),
                                     output_rank )
        elif sys_num is not None:
            system.get_data_by_num( self.p, sys_num, var_names, var_vals )
        else:
            raise ValueError( 
              'EquationSystems.sys_to_python must be given either sys_num or sys_name.')

        ### Put the data into dictionaries that will be easily converted
        ### to DataFrames
        out_data = {}
        for i, var in enumerate( var_names ):
            out_data[var] = list(var_vals[i::len(var_names)])

        return out_data

################################################################################
