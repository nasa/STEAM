""" Module for the Solution Class
"""
import logging
logger = logging.getLogger(__name__)
logger.debug("Top of file")

import pandas as pd
import numpy  as np
import steam.util

from code import interact
#interact( local = dict( globals(), **locals() ) )

### Import the solution storage enum if libmesh is available
#from steam import has_libmesh
#if steam.has_libmesh:
#from steam.libmesh.equation_systems import ElemNodeStorage as DataAt

class Solution():
    """Class definition for STEAM Solution

    The solution is stored in a DataFrame (data).

    """

    def __init__(self,mesh=None,file=None,filetype=None,copy=None):
        """Constructor for Solution class.

            copy=MESH will copy a solution (and reference its mesh)
        """

        self.data   = None            # This will store the actual data in a DataFrame
        self.info   = ""              # String descriptor of information
        self.mesh   = mesh            # The mesh that this solution is based on
        self.hash   = ""              # The mesh's hash

        self.store_at = None          # Data is at 'NODES' or 'ELEMENTS'

        if (file is not None):
            self.read(file,filetype)

        if (copy is not None):
            self.copy(copy)

    def __eq__( self, other ):
        """ Used for equality determination.

        Note that the meshes are not compared.  This prevents getting into
        an infinite loop when static solutions are compared between meshes.

        Args:
            other (:obj:`~steam.solution.Solution`): The solution with which this one is to be compared.
        """

        ### Check that enough stuff can at least be compared
        if ( self.__class__ != other.__class__ or
             type( self.data ) != type( other.data ) ):
           return False

        if isinstance( self.data, pd.DataFrame ):
            return (
                    self.data.equals( other.data ) and
                    self.store_at == other.store_at
                   )
        elif self.data is None:
            ### We already know that other.data must be NoneType
            return self.store_at == other.store_at
        else:
            raise TypeError('soln.data should be either pd.DataFrame or None')
            

    def __neq__( self, other ):
        """ Used in comparison of self != other

        Args:
            other (:obj:`~steam.solution.Solution`): The solution with which this one is to be compared.
        """
        return not self == other

    def __str__(self, print_mesh = True):
        """Print to screen."""

        string  = "Python Solution Object:\n"

        if (self.data is None):
            string += "   Data Size       : Not in memory\n"
        else:
            string += "   Data Size       : {}\n".format(len(self.data))
            vars = self.data.columns.tolist()
            var_string = " ".join(str(i) for i in vars)
            string += "   Variables       : {}\n".format(var_string)
            
            self.set_store_flag()

        if self.store_at is None:
            string += "   Neither Element- nor Node-based\n"
        elif self.store_at.upper().find('ELEM') == 0:
            string += "   Element-based\n"
        elif self.store_at.upper().find('NODE') == 0:
            string += "   Node-based\n"
        else:
            raise ValueError( 'solution.store_at must be either "ELEMENTS", '
                            + '"NODES", or None.' )

        if (self.mesh is None):
            string += "   Mesh            : None\n"
        elif print_mesh:
            if (self.mesh.is_point_cloud):
                string += "   Mesh (PointCloud): \n"
            else:
                string += "   Mesh             : \n"
            ### Format the mesh to be indented
            import re
            grid_string = self.mesh.__str__( print_static = False )
            grid_string = re.sub( '^','   '  ,grid_string)
            grid_string = re.sub('\n','\n   ',grid_string)
            #! Remove some lines we don't want
            grid_string = re.sub('^.*Mesh Object.*\n','',grid_string) 
            #grid_string = re.sub('^.*Hash.*\n','',grid_string) 
            string += grid_string
        else:
            string += "   Mesh is populated.\n"

        return string

    def copy(self,copy_from):
        """ Copy the data in 'copy_from' into self.
        
        This does not make a copy of the mesh."""
        from copy import deepcopy

        #! This is silly.  But I don't know how else to copy the values
        def df_deepcopy(df):
            #! Replace the data
            copy_df = pd.DataFrame(
                                   deepcopy(df.values),
                                   deepcopy(df.index),
                                   deepcopy(df.columns)
                                   )
            return copy_df

        #! Copy the Pandas objects
        if (copy_from.data   is not None):
            self.data   = df_deepcopy(copy_from.data  )

        #! Now copy the non-Pandas objects
        self.info       = deepcopy(copy_from.info      )
        self.hash       = deepcopy(copy_from.hash      )

        self.store_at   = deepcopy(copy_from.store_at   )

        #! The mesh is only referenced
        self.mesh       = copy_from.mesh
        return

    def init(self):
        """ Make solution data consistent after setting it."""

        if self.mesh is not None:
            self.set_store_flag()

        return

    def set_store_flag(self, mesh=None):
        """ Determine whether the solution is element- or node-based.
        
        This compares the data size to the mesh size.  Set self.store_at 
        flag accordingly.  self.store_at is intended to work in conjunction 
        with the swigged DataAt enum.

        Note that it's possible to have self.store_at set to neither 'ELEMENTS'
        nor 'NODES' because solutions can exist without an awareness of 
        where their data is.

        Args:
            mesh(:obj:`~steam.mesh.Mesh`): Mesh to compare against; default is self.mesh.
        """

        if (mesh is None):
            mesh = self.mesh

        if mesh is None:
            ### There's no way to know 
            self.store_at = None
            return

        ### Allow for the case when mesh is old and doesn't have the 
        ### 'is_point_cloud' flag
        #try:
        if mesh.is_point_cloud:
            self.store_at = 'NODES'
            return

        if (mesh.xyz_pt.shape[0] == mesh.conn.shape[0]):
            if mesh.conn.shape[1] == 1:
                # This should be a point cloud, connectivity is one-to-one
                self.store_at = 'NODES'
                return
            raise Exception("Mesh type is ambiguous - need better check!")

        if (self.data.shape[0] == mesh.xyz_pt.shape[0]):
            self.store_at = 'NODES'
        elif (self.data.shape[0] == mesh.conn.shape[0]):
            self.store_at = 'ELEMENTS'
        else:
            self.store_at = None
            logger.error( 'Inconsistent shape between mesh and solution.' +
                          '  Unable to determine and set self.store_at' )

        return 


    def read(self,file,filetype='CDAT'):
        """Method to read a new solution object.
 
        Args:
            file (:obj:`str`): Path to file.
            filetype (:obj:`str`): Type of file. See table below:

        Here are  the file types::

                DF    - DataFrame, space delinated file with column headers
                TRIQ  - TRIQ ASCII
                uTRIQ - TRIQ fortran unformated
                PKL   - Pickle file
                *    - Fail
        """
 
        ftu = filetype.upper()
        if   (ftu == "DF" or ftu == "CDAT"):
            self.read_cdat(file)
        elif (ftu == "TRIQ"):
            steam.io.read_triq_ascii(file,soln=self)
        elif (ftu == "UTRIQ"):
            try:
                steam.io.read_triq_uform(file,soln=self)
            except FloatingPointError:
                # I assume double precision, if single, then try again
                steam.io.read_triq_uform(file,soln=self,dp=False)
        elif (ftu == "PKL"):
            steam.io.read_pickle(file,soln=self)
        elif (ftu == "TECPLOT_ASCII"):
            steam.io.read_tecplot_ascii(file,soln=self)
        else:
            raise IOError( f'Solution filetype not recognized: {filetype}' )

        # Init solution variables
        self.init()

    def update_hash(self):
       """ Update the solution's hash using xxhash."""
       self.hash = steam.util.hasher(self.data.values.copy(order='C'))


    def vars(self):
        """ Return a list of the solution data variables.

        Returns:
            (:obj:`list`): List of variables
        """
        return self.data.columns.tolist()

    def rename_vars(self,vars,inplace=True):
        """ Rename variables in the solution data.

        Pass a dictionary with the entries being "BEFORE":"AFTER".
        Variables are renamed in place.

        Args:
            vars (:obj:`dict`): Dictionary of variables to rename
            inplace (:obj:`bool`): Default true, rename in place or create new solution data
        """
        self.data.rename(columns=vars,inplace=inplace)

        return

    def remove_vars(self,vars):
        """ Remove variables from solution data.

        Args:
            vars (:obj:`list`): List of variables to remove
        """

        for var in vars:
            try:
                self.data.drop(var, axis=1, inplace=True)
            except:
                raise Exception("Variable '{}' not found in solution data!".format(var))

        self.init()
        return

    def add_var(self, var, data, index_in=None):
        """ Add a new variable based on a list of data.

        Args:
            var  (:obj:`str`): Name of variable.
            data (:obj:`list`): List of data, needs to match existing size (point/element).
            index_in (:obj:`pandas.core.indexes.range.RangeIndex`): Optional index to apply to the `data` input.  This input is necessary when appending a variable to a solution with a discontinuous or when a variable will not be defined at all nodes/elements in the solution.
            
        """

        new_df = pd.DataFrame( list(data), columns=[var], index=index_in )

        if (self.data is None):
            self.data = new_df
            return

        if (var in self.data.columns):
            raise ValueError( 
                f"Variable {var} already exists.  Remove it before adding." )

        self.data = self.data.join(new_df)


    def check_is_point(self,mesh=None):
        """ Return True if this is point-based data.
        
            This compares the data size to the mesh size in order to check 
            what is in the .point and .element attributes of the object.

            Args:
                mesh(:obj:`~steam.mesh.Mesh`): Mesh to compare against; default is self.mesh.
        """

        check = False

        if (mesh is None):
            mesh = self.mesh

        if mesh.is_point_cloud:
            return True

        if (mesh.xyz_pt.shape[0] == mesh.conn.shape[0]):
            raise Exception("Mesh type is ambiguous - need better check!")

        if (self.data.shape[0] == mesh.xyz_pt.shape[0]):
            check = True

        return check

    def check_is_element(self,mesh=None):
        """ Return True if this is element-based data. 
        
            This compares the data size to the mesh size in order to check 
            what is in the .point and .element attributes of the object.
        
            Args:
                mesh(:obj:`~steam.mesh.Mesh`): Mesh to compare against; default is self.mesh.
        """

        check = False

        if (mesh is None):
            mesh = self.mesh
        
        if mesh.is_point_cloud:
            return False

        if (mesh.xyz_pt.shape[0] == mesh.conn.shape[0]):
            raise Exception("Mesh type is ambiguous - need better check!")

        if (self.data.shape[0] == mesh.conn.shape[0]):
            check = True

        return check
        

    def blank_copy(self):
        """ Return a copy of this solution with zero-ed data.

        Returns:
            (:obj:`Solution`): Copy of solution with zero-ed data
        """

        copy = Solution(self.mesh)

        shape = self.data.shape
        cols  = self.data.columns

        copy.data = pd.DataFrame(np.zeros(shape),columns=cols)
        copy.info   = "Copy of "+self.info

        # Init solution variables
        self.init()

        return copy

################################################################################
### IO
    def read_cdat(self,file, usecols = None):
        """ Read CDAT style dataframe.
        
        Args:
            file (:obj:`str`): Path to file.
            usecols (:obj:`list`,optional): Column names to read from file.  Default behavior is to read all columns.
        """
        self.data = pd.read_table(file,delim_whitespace=True, usecols=usecols)
        self.store_at = 'NODES'
        return


#!  -------- Format Conversions ---------

    def node_to_element(self,return_new=False):
        """ Convert node data to element data.

            This merely averages the data at the vertices to calcluate values at elements.

            Args:
                return_new (:obj:`bool`): Return a new object instead of replacing self? Default False.

            Returns:
                (:obj:`~steam.solution.Solution`): Option return of element-based solultion.
        """

        #! First check to make sure that sizes make sense:
        if (self.data.shape[0] != self.mesh.xyz_pt.shape[0]):
            raise Exception("Soln is not node based!")

        #! Look at each element of the grid and average the solution data
        #! at each vertex.
        elements  = self.mesh.conn.values;
        node_mat  = self.data.values;
        elem_mat  = np.zeros([self.mesh.conn.shape[0],self.data.shape[1]])

        for (e,elem) in enumerate(elements):
            for node in elem:
                elem_mat[e] += node_mat[node]
            elem_mat[e] /= float(elem.size)

        # I think this is faster:
        #elem_mat = np.sum(node_mat[elements],axis=1)/elements.shape[1]

        #! Replace the data
        vars = self.data.columns.tolist()
        elem_df  = pd.DataFrame(
                                elem_mat,
                                columns=vars)

        if (return_new):
            soln         = Solution();
            soln.mesh    = self.mesh
            soln.data    = elem_df
#            soln.point   = False
#            soln.element = True
            soln.store_at = 'ELEMENTS'
            soln.init()
            return soln

        self.data    = elem_df
        self.store_at = 'ELEMENTS'
        # Init solution variables
        self.init()

        return

 
    def element_to_node(self,return_new=False):
        """ Convert element data to node data.

            This merely averages the data at neighboring elements to the vertices.

            Args:
                return_new (:obj:`bool`): Return a new object instead of replacing self? Default False.

            Returns:
                (:obj:`~steam.solution.Solution`): Option return of node-based solultion.
        """

        #! First check to make sure that sizes make sense:
        if (self.data.shape[0] != self.mesh.conn.shape[0]):
            raise Exception("Soln is not element based!")

        #! Look at each element of the grid and average the solution data
        #! at each vertex.
        
        elements  = self.mesh.conn.values;
        num =  np.zeros(self.mesh.xyz_pt.shape[0])
        node_mat  = np.zeros([self.mesh.xyz_pt.shape[0],self.data.shape[1]])
        elem_mat  = self.data.values;

        for (e,elem) in enumerate(elements):
            for node in elem:
                node_mat[node] += elem_mat[e]
                num[node]      += 1;

        for node in range(self.mesh.xyz_pt.shape[0]):
            node_mat[node]     /= float(num[node])

        #! Has to be a faster way
        #for n in range(self.mesh.xyz_pt.shape[0]):
        #    node_mat[n] = sum(elem_mat[np.any(np.where(elements==n,True,False),axis=1)]) / \
        #                  np.sum(np.where(elements==n,1,0))
        

        #! Replace the data
        vars = self.data.columns.tolist()
        node_df   = pd.DataFrame(
                                 node_mat,
                                 columns=vars)
        if (return_new):
            soln = Solution();
            soln.data = node_df
            soln.mesh = self.mesh
            soln.init()
            soln.store_at = 'NODES'
            return soln
    
        self.data = node_df
        self.store_at = 'NODES'
        self.init()

        return

    def return_subset(self,el_list=[]):
        """Return components from the solution as a new solution.

        This will return the data for elements that are part of a
        list.  Indexing for the elements is not changed in order to 
        maintain compatability between the original mesh's elements.

        Args:
            el_list : List of elements to return.

        Returns:
            (:obj:`solution`): Requested subset of the solution.
        """
    
        #! This needs to be an element-based solution
        if ( self.store_at.upper().find('ELEM') != 0 and
                             not self.mesh.is_point_cloud ):
            raise Exception("The solution needs to be element-based")
            return

        soln = steam.solution.Solution()

        #! I'm using .index.isin since it does two things:
        #! - sorts the elments and ignores out-of-bounds
        def df_sub(df,list):
            return df[df.index.isin(list)]

        soln.data = df_sub(self.data,el_list)

        ### Include a subset of the static solution?
        #include = False if self.mesh.static_soln is None else True
        #soln.mesh = self.mesh.return_subset(el_list, include_static = include )
        soln.mesh = self.mesh.return_subset(el_list, include_static = False )
        soln.init()

        return soln


    def return_comp(self,comp_list=[]):
        """Return components from the solution as a new solution.

        This will return the data for elements that are part of a certain
        component or components from the mesh.  Indexing for the
        elements is not changed in order to maintain compatability
        between the original mesh's elements.

        Args:
            comp_list : List of components to return.

        Returns:
            (:obj:`solution`): Requested subset of the solution.
        """
    
        #! This needs to be an element-based solution
        if ( self.store_at.upper().find('ELEM') != 0 and 
                             not self.mesh.is_point_cloud ):
            raise Exception("The solution needs to be element-based")
            return

        soln = steam.solution.Solution()

        el_list   = []  # Elems to get
    
        for comp in comp_list:
            elems = self.mesh.get_comp(comp)
            el_list.extend(elems)
    
        el_list = list(set(el_list))

        return self.return_subset(el_list)

## Begin HDF5 Routines

    def write_hdf5(self,hdf5,root="/",options=None):
        """ Write the solution in memory to disk.
        
        This will overwrite any data that is in the HDF5 file.

        Args:
            hdf5 (`~pandas.io.pytables.HDFStore`): HDFStore object.
            root (:obj:`str`,optional): Location in HDF5 to write database.
                Default is to write to root location "/".
            options (:obj:`dict`, optional): options to vary behavior.

        Current options are:
            ioMesh  = True/False : Should the meshes be written? [T]
        """

        try:
            ioMesh = options['ioMesh']
        except:
            ioMesh = True

        hdf5.put(
                 root,
                 self.data,
                 format='t',
                 )

        ### Make objects that will store attributes
        iroot = root
        hdf5.get_storer(iroot).attrs.info     = self.info
        hdf5.get_storer(iroot).attrs.store_at = self.store_at
        hdf5.get_storer(iroot).attrs.hash     = self.hash

        if (ioMesh):
            self.mesh.write_hdf5(hdf5,root+"/mesh")

        return

    def read_hdf5(self,hdf5,root="/",options=None):
        """ Read the solution from the HDF5 file to memory.

        Args:
            hdf5 (`~pandas.io.pytables.HDFStore`): HDFStore object.
            root (:obj:`str`,optional): Location in HDF5 to write database.
               Default is to write to root location "/".
            options (:obj:`dict`, optional): options to vary behavior.

        Current options are:
            ioMesh  = True/False : Should the meshes be read? [T]
            readData  = True/False : Should the data   be read? [T]
       """

        try:
            ioMesh = options['ioMesh']
        except:
            ioMesh = True

        try:
            readData = options['readData']
        except:
            readData = True


        ### Check to make sure that this database has a solution
        exist = hdf5.get_node(root)
        if exist is None:
            #print(" -- No solution found in HDF5 File:",root)
            logger.warn(" -- No solution found in HDF5 File:",root)
            return None

        if (readData):
            self.data = hdf5.get(root)

        if (ioMesh):
            self.mesh = steam.mesh.Mesh()
            self.mesh.read_hdf5(hdf5,root+"/mesh")

        iroot = root
        self.info    = hdf5.get_storer(iroot).attrs.info
        self.store_at= hdf5.get_storer(iroot).attrs.store_at
        #! Solutions in databases have hashes
        try:
            self.hash       = hdf5.get_storer(iroot).attrs.hash
        except AttributeError:
            pass

        return self

    @classmethod
    def from_libmesh(cls,l_es, name = None, num = None, mesh = None, 
                     local_or_entire="LOCAL_SOLN", output_rank=0 ):
        """ Constructor to convert from libmesh object.
        
        Args:
            l_es (:obj:`~steam.libmesh.equation_systems()`): libmesh equation_systems object.
            name(:obj:`string`): Name of the system to be recalled from libMesh.
            num(:obj:`string`): Number of the system to be recalled from libMesh.
            mesh (:obj:`~steam.mesh.Mesh`): Python mesh object.
            local_or_entire(:obj:`string`): Flag indicating whether to return the solution local to the rank or the entire solution.  Acceptable values: "LOCAL_SOLN" or "ENTIRE_SOLN"
            output_rank(:obj:`int`): The rank that will be returning the entire solution if local_or_entire=="ENTIRE_SOLN"
        Returns:
            :obj:`~steam.solution.Solution`: python solution object.
        """

        data_dict = l_es.sys_to_python( sys_name = name, sys_num = num, 
                                        local_or_entire=local_or_entire,
                                        output_rank=output_rank )

        soln = cls()
        soln.data = pd.DataFrame.from_dict( data_dict )

        ### Add the mesh to the solution if possible
        if mesh is None:
            ### Get the mesh from the provided EquationSystems
            #l_mesh = steam.libmesh.Mesh( pointer = l_es.get_mesh_ptr() )
            l_mesh = steam.libmesh.Mesh( pointer = l_es.parent.p )
            soln.mesh = steam.mesh.Mesh.from_libmesh( l_mesh )
        else:
            ### Take a python Mesh object as an input
            n_mesh_vals = -999
            if soln.store_at.upper().find( 'ELEM' ) == 0:
                n_mesh_vals = len( mesh.conn.index )
            elif soln.store_at.upper().find( 'NODE' ) == 0:
                n_mesh_vals = len( mesh.xyz_pt.index )
            else:
                raise ValueError( 
                      'soln.store_at must match either ELEM or NODE' )

            ### Confirm that data is the right size for the mesh
            if len( soln.data ) == n_mesh_vals:
                soln.mesh = mesh
            else:
                raise ValueError( 
                      'mesh size is not compatible with solution.\n mesh ' +
                      'has {:d} {:s} and solution has {:d} {:s}\n\n'.format(
                      n_mesh_vals, soln.store_at.lower(), soln.data.shape, 
                      soln.store_at.lower() ))

        soln.init()

        return soln

################################################################################

def uniform_soln( mesh, node_elem_flag = 'ELEMENTS', values = [1.0] ):
    """Create a uniform solution for a given mesh.  This function is 
    valuable for unit testing.

    Args:
        mesh (:obj:`~steam.mesh.Mesh`): Python mesh object.
        node_elem_flag (:obj:`str`): Flag indicating whether the output 
            solution will store data at elements ('ELEM') or nodes ('NODE')
        values (:obj:`list`): List of uniform values to be applied everywhere.
               NOTE: one or more entries can be "rand" -- indicating that that
               variable should be populated with random numbers between 0 and
               1.  This option can be useful for debugging or unit testing.
    Returns:
        :obj:`~steam.solution.Solution`: Python solution object.
    """

    out_soln = steam.solution.Solution( mesh = mesh )

    npts = len( mesh.xyz_pt.index )
    nele = len( mesh.conn.index )
    nq   = len( values )
    qs   = {}

    if node_elem_flag.upper().find( 'ELEM' ) == 0:
        ### Storing data at elements
        nloop = nele
        out_soln.store_at = 'ELEMENTS'
        mesh.get_xyz_el()
        index_source = mesh.xyz_el
    elif node_elem_flag.upper().find( 'NODE' ) == 0:
        ### Storing data at nodes
        out_soln.store_at = 'NODES'
        nloop = npts
        index_source = mesh.xyz_pt
    else:
        raise ValueError( f'node_elem_flag must be either ELEM or NODE, '
                        + f'not {node_elem_flag}' )

    ### Set solution variables
    for j in range(nq):
        num = j + 1
        #qs["q{}".format(num)] = []
        try:
            qs["q{}".format(num)] = list( np.ones(nloop) * values[j] )
        except TypeError:
            ### Replace uniform value with random numbers
            if isinstance(values[j], str) and values[j].lower() == 'rand':
                qs["q{}".format(num)] = list( np.random.rand(nloop) )
            else:
                raise ValueError( "Entries in 'values' must be numbers " +
                                  'or "rand".' )

#    interact( local = dict( globals(), **locals() ) )
    out_soln.data = pd.DataFrame.from_dict(qs)
    out_soln.data.set_index( index_source.index, inplace=True )

    # Init solution quantities
    out_soln.init()

    return out_soln
            
def interp_half_body_soln_to_full_body_grid( soln, mesh, vars=[], plane=None, tolerance=1.0e-6 ):
    """Interpolate half body solution data onto a full body mesh.

    This is currently implemented to handle element based data onto mesh elements, 
    so it will not handle node based data or a point cloud.


    Args:
        soln (:obj:`~steam.solution.Solution`): STEAM solution object
        mesh (:obj:`~steam.mesh.Mesh`):         STEAM mesh object
        vars (:obj:`list`):                     List of variables names that need to be mirrored (ex: velocity or shear component)
        tolerance (:obj:`float`):               Float value used to determine the differences between soln and mesh to automatically find the mirror plane.

    Returns:
        :obj:`~steam.solution.Solution`: STEAM solution object.
    """

    mesh.get_xyz_el()

    # Check that each of the vars is in the soln
    soln_vars = soln.data.columns.tolist()
    for var in vars:
        assert (var in soln_vars), "Variable ({}) was not found in the solution variables".format (var)

    # Determine the plane to mirror across
    if plane is None:
        xyz_vars = soln.mesh.xyz_pt.columns.tolist ()
        soln_bounds = [0 for i in range (6)]
        soln_bounds [0] = soln.mesh.xyz_pt.min()[xyz_vars[0]]
        soln_bounds [1] = soln.mesh.xyz_pt.max()[xyz_vars[0]]
        soln_bounds [2] = soln.mesh.xyz_pt.min()[xyz_vars[1]]
        soln_bounds [3] = soln.mesh.xyz_pt.max()[xyz_vars[1]]
        soln_bounds [4] = soln.mesh.xyz_pt.min()[xyz_vars[2]]
        soln_bounds [5] = soln.mesh.xyz_pt.max()[xyz_vars[2]]
    
        mesh_xyz_vars = mesh.xyz_pt.columns.tolist ()
        mesh_bounds = [0 for i in range (6)]
        mesh_bounds [0] = mesh.xyz_pt.min()[mesh_xyz_vars[0]]
        mesh_bounds [1] = mesh.xyz_pt.max()[mesh_xyz_vars[0]]
        mesh_bounds [2] = mesh.xyz_pt.min()[mesh_xyz_vars[1]]
        mesh_bounds [3] = mesh.xyz_pt.max()[mesh_xyz_vars[1]]
        mesh_bounds [4] = mesh.xyz_pt.min()[mesh_xyz_vars[2]]
        mesh_bounds [5] = mesh.xyz_pt.max()[mesh_xyz_vars[2]]
    
        comparisons = [True for i in range (6)]
        for i in range (6):
            if abs (soln_bounds [i] - mesh_bounds [i]) > tolerance: comparisons [i] = False
            else: comparisons [i] = True
        if comparisons.count (False) != 1:
            return
        else:
            mirror_plane_index = comparisons.index (False)
    
        if mirror_plane_index == 0:
            plane = xyz_vars [0] + str ('<') + str (soln_bounds[0])
        elif mirror_plane_index == 1:
            plane = xyz_vars [0] + str ('>') + str (soln_bounds[1])
        elif mirror_plane_index == 2:
            plane = xyz_vars [1] + str ('<') + str (soln_bounds[2])
        elif mirror_plane_index == 3:
            plane = xyz_vars [1] + str ('>') + str (soln_bounds[3])
        elif mirror_plane_index == 4:
            plane = xyz_vars [2] + str ('<') + str (soln_bounds[4])
        elif mirror_plane_index == 5:
            plane = xyz_vars [2] + str ('>') + str (soln_bounds[5])

#    print (plane)
    # Determine elements on mesh that are on opposite side of plane from soln    
    mesh_half_mirror_elems = mesh.xyz_el.query (plane).index
    
    try:
        plane = plane [:plane.index ('<')]
    except ValueError:
        plane = plane [:plane.index ('>')]

    # Mirror the mesh element xyz values
    mesh_half = steam.mesh.Mesh (copy=mesh)
    mesh_half.xyz_el[plane].values[mesh_half_mirror_elems] *= -1.0

    # Interpolate the solution to the new mesh
    trans = steam.interpolate.Transform (source=soln.mesh, target=mesh_half)
    trans.inverse_distance ()
    new_soln = trans.apply(soln)

    # Mirror the vector component variables
    for var in vars:
        new_soln.data[var].values[mesh_half_mirror_elems] *= -1.0
    new_soln.mesh = mesh

    return new_soln

###############################################################################

def half_to_full( full_mesh, half_soln, flip_vars=[] ):
    """ Take a half-body Solution object and return a full-body version.

    NOTE: The difference between this function and 
          `~steam.solution.interp_half_body_soln_to_full_body_grid` is that 
          this function does no interpolation.  Instead it mirrors the solution
          and exactly copies it onto to a symmetrical mesh.

    Many CFD solutions and databases are constructed using only half the 
    vehicle in order to save disk space and computational time by leveraging
    symmetry.  When a full-body solution is needed from a half-body one, this
    function (along with `~steam.mesh.half_to_full`) can be used to deliver it.

    Args:
        full_mesh (:obj:`~steam.mesh.Mesh): The full body mesh that will be associated with the solution.  This mesh should have been generated by `~steam.mesh.half_to_full`.
        half_soln (:obj:`~steam.solution.Solution`): Mesh to mirror
        flip_vars (:obj:`list`): Variables whose values should be multiplied by -1.0 on the mirrored side of the vehicle.  Useful for components of shear and similar parameters.

    Returns:
        :obj:`~steam.solution.Solution`: The full-body solution
    """

#    raise NotImplementedError( 'Not done yet.' )

    if half_soln.store_at.upper().startswith( 'ELEM' ):
#        interact( local = dict( globals(), **locals() ) )
        map_series = full_mesh._elem_map
#        raise NotImplementedError( 'working on it.' )
    elif half_soln.store_at.upper().startswith( 'NODE' ):
        map_series = full_mesh._node_map
#        raise NotImplementedError( 'working on it.' )
    else:
        raise ValueError( 'half_soln must have data stored at either elements '
                        + 'or nodes.  \n\nCurrently half_soln.store_at = '
                        + f'{half_soln.store_at}' )

    full_soln = steam.solution.Solution( full_mesh )

    ### Constructed the mirrored side of the solution
    mirror_data = half_soln.data.loc[ map_series.index ]
    #   The index needs to be replaced with the new element IDs
    mirror_data.set_index( map_series.values, inplace=True )
    for var in flip_vars:
        mirror_data[var] *= -1.0

    ### Original data and mirrored data need to be reassembled
    full_soln.data = pd.concat( (half_soln.data, mirror_data) )
    full_soln.set_store_flag()

    return full_soln

###############################################################################
