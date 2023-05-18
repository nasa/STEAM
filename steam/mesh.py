""" Module for the Mesh Class
"""
import logging
logger = logging.getLogger(__name__)
logger.debug("Top of file")

import pandas as pd
import numpy  as np
import steam.util
import steam.libmesh

from code import interact
#interact(local= dict( globals(), **locals()) )

# This turns on some extra debugging in this module if you need it
module_debug = False


class Mesh:
    """Class definition for STEAM Mesh

    The mesh is stored in a DataFrame with minimal values of
    X, Y, and Z.

    A mesh can be converted into a point cloud (based on verticies).
    Element connectivity is forced to be a simple one-to-one to the nodes.

    """

    def __init__(self,file=None,filetype=None,copy=None,is_point_cloud=False,
                 **kwargs):
        """Constructor for Mesh class.

        Args:
            file (:obj:`str`, optional): Path to file from which to read the mesh.
            filetype (:obj:`str`, optional): Flag indicating the type of file to read. See table below for allowable values.
            copy (:obj:`~steam.mesh.Mesh`, optional): Mesh to copy.
            is_point_cloud(:obj:`bool`): Flag indicating whether or not this mesh is a point cloud (without connectivity).  Defaults to "False".

        Here are  the file types::

                DF   - DataFrame, space delinated file with column headers
                CDAT - DataFrame, space delinated file with column headers
                DAT  - Tecplot ASCII File
                TRI  - ASCII .tri  File
                TRIQ - ASCII .triq File
                PKL  - Pickled Mesh File 
                *    - Pass to libmesh
        """

        self.xyz_pt     = None          # This will store the node centroid DataFrame
        self.xyz_el     = None          # This will store the element centroid DataFrame
        self.conn       = None          # This will store the element connectivity for FE grids (0-based)

        self.comp       = None          # This will store the element components   for FE grids
        self.comp_table = dict()        # Table of component associations.
        self.info       = ""            # String descriptor of information
        self.hash       = None

        self.static_soln = None         # House information about the points in the mesh if necessary


        ### Set is_point_cloud before reading
        self.is_point_cloud = is_point_cloud

        if (file is not None):
            self.read(file,filetype, **kwargs)

        if (copy is not None):
            self.copy(copy)

        if is_point_cloud:
            self.to_point_cloud()

    def __str__(self, print_static = True ):
        """Print to screen."""
        if (self.is_point_cloud):
            string  = "Python Mesh Object (PointCloud):\n"
        else:
            string  = "Python Mesh Object:\n"
        string += "   Points    : {}\n".format(len(self.xyz_pt))
        string += "   Elements  : {}\n".format(len(self.conn))
        comps = set(self.comp[0])
        string += "   # of Comps: {}\n".format(len(comps) + len(self.comp_table))
        comp_string = " ".join(str(i) for i in comps)
        string += "   Int Components: {}\n".format(comp_string)
        comp_string = " ".join(str(i) for i in self.comp_table)
        string += "   Str Components: {}\n".format(comp_string)

        ### Print the static solution
        if self.static_soln is  None:
            string += "   No Static Solution.\n"
        elif print_static:
            string += "   Static Solution:\n"

            ### Format the static solution to be indented
            import re
            soln_string = self.static_soln.__str__( print_mesh = False )
            soln_string = re.sub( '^','      '  ,soln_string)
            soln_string = re.sub('\n','\n      ',soln_string)
            #! Remove some lines we don't want
            #soln_string = re.sub('^.*Mesh Object.*\n','',soln_string) 
            string += soln_string

        else:
            string += "   Static Solution Exists.\n"

        return string

    def __eq__( self, other, compare_static = True, compare_xyz_el = False ):
        """ Used for equality determination.

        There are several comparisons possible with this function, but only
        when called explicitly, rather than with the "==" operator.
        """

        ### Two meshes with identical xyz_pt and conn should match, so
        ### it's not necessary for xyz_el to match.
        if compare_xyz_el:
            ### It's possible to have a mesh that doesn't contain xyz_el
            if self.xyz_el is None:
                if other.xyz_el is not None:
                    return False
            else:
                if other.xyz_el is None:
                    return False
                if not self.xyz_el.equals( other.xyz_el ):
                    return False

        ### Main comparison
        if ( self.__class__ == other.__class__ and
             self.xyz_pt.equals( other.xyz_pt ) and
             self.conn.equals( other.conn ) and
             self.comp_table == other.comp_table and
             self.is_point_cloud == other.is_point_cloud ):

            if compare_static:
                return self.static_soln == other.static_soln
            else:
                return True
        else:
            return False

    def __neq__( self, other, **kwargs ):
        """ Used in comparison of self != other

        Args:
            other (:obj:`~steam.mesh.Mesh()`): The mesh with which this one is to be compared.
        """
        return not steam.mesh.Mesh.__eq__( self, other, **kwargs )

    def copy(self,copy_from):
        """ Copy the data in 'copy_from' to self.
        
        Args:
            copy_from (:obj:`~steam.mesh.Mesh): The mesh to copy
        """

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
        if (copy_from.xyz_pt is not None):
            self.xyz_pt = df_deepcopy(copy_from.xyz_pt)
        if (copy_from.xyz_el is not None):
            self.xyz_el = df_deepcopy(copy_from.xyz_el)
        if (copy_from.conn   is not None):
            self.conn   = df_deepcopy(copy_from.conn  )

        #! Now copy the non-Pandas objects
        self.comp       = deepcopy(copy_from.comp)
        self.comp_table = deepcopy(copy_from.comp_table)
        self.info       = deepcopy(copy_from.info      )
        self.hash       = deepcopy(copy_from.hash      )

        ### Copy the static solution
        if copy_from.static_soln is not None:
            if self.static_soln is None:
                self.static_soln = steam.solution.Solution()
            self.static_soln.copy( copy_from.static_soln )

        return

    def update_hash(self):
        """ Update the mesh's hash using xxhash.
        
        Returns:
            :obj:`str`: A hash of the xyz_pt and elements."""

        #! This is the guy right here:
        #! - I'd check xyz_el and xyz_pt, but I'm not 100% sure that we will
        #!   always have xyz_el populated
        xyz_hash  = steam.util.hasher(self.xyz_pt.values.copy(order='C'))
        conn_hash = steam.util.hasher(  self.conn.values.copy(order='C'))
        self.hash = xyz_hash+conn_hash

        ### A mesh with a static_soln differs from one without, even if
        ### all other things are equal
        if self.static_soln:
            self.static_soln.update_hash()
            self.hash += self.static_soln.hash

    @classmethod
    def from_libmesh(cls,lmesh):
        """ Constructor to convert from libmesh object.
        
        Args:
            lmesh (:obj:`~steam.libmesh.Mesh`): libmesh mesh object.
        Returns:
            :obj:`~steam.mesh.Mesh()`: python mesh object.
        """

        nodes, elems, comps = lmesh.to_python()

        mesh = cls()
        mesh.xyz_pt = pd.DataFrame.from_dict (nodes)
        mesh.conn   = pd.DataFrame.from_dict (elems, 'index')
        mesh.comp   = pd.DataFrame.from_dict (comps, 'index')

        return mesh

    def read(self,file,filetype=None, **kwargs):
        """Method to read a new mesh object.
 
        Args:
            file (:obj:`str`): Path to file.
            filetype (:obj:`str`): Type of file. See table below:

        Here are  the file types::

                DF   - DataFrame, space delinated file with column headers
                CDAT - DataFrame, space delinated file with column headers
                DAT  - Tecplot ASCII File
                TRI  - ASCII .tri  File
                TRIQ - ASCII .triq File
                PKL  - Pickled Mesh File 
                *    - Pass to libmesh
        """
        ftu = filetype.upper()
        if (ftu == "CDAT" or ftu == "DF"):
            self.read_cdat(file, **kwargs)
            if self.is_point_cloud:
                self.static_soln = steam.solution.Solution( mesh = self, 
                                        file = file, filetype = filetype )
        elif (ftu == "DAT"):
            self.read_tecplot_ascii(file)
        elif (ftu == "TRI" or ftu == "TRIQ"):
            if self.is_point_cloud:
                self.static_soln = steam.io.read_triq_ascii(file,
                                                    ioSoln=True,mesh=self )[1]
            else:
                steam.io.read_triq_ascii(file,ioSoln=False,mesh=self)
        elif (ftu == "UTRI" or ftu == "UTRIQ"):
            if self.is_point_cloud:
                self.static_soln = steam.io.read_triq_uform(file,
                                                    ioSoln=False,mesh=self)[1]
            else:
                steam.io.read_triq_uform(file,ioSoln=False,mesh=self)
        elif (ftu == "PKL"):
            steam.io.read_pickle(file,mesh=self)
        else:
            self.read_libmesh(file)

        if (self.comp is None):
            print("* Error: mesh filetype not recognized")
            raise IOError

    def read_cdat(self,file, coord_vars = ["X","Y","Z"], **kwargs ):
        """ Read CDAT style dataframe.
        
        Args:
            file (:obj:`str`): Path to file.
            coord_vars (:obj:`list`): List containing names of X, Y, and Z coordinates in file.
        """

        final_vars = ["X","Y","Z"]

        self.xyz_pt = pd.read_table(file,delim_whitespace=True,
                                    usecols = coord_vars,
                                    )

        ### Reformatted variable checking in more Pythonic style
        for var in coord_vars:
            error_message = '{:s} missing in grid file: {:s}'.format(var, file)
            try:
                assert var in self.xyz_pt.columns, error_message
            except AssertError:
                logger.error( error_message )
                raise
        
        ### Swap the variable names if necessary
        if coord_vars != final_vars:
            swap_dict = dict( zip( coord_vars, final_vars ) )
            self.xyz_pt.rename( index=int, columns=swap_dict, inplace=True )

        # Check for bad values
        bcheck = pd.isnull(self.xyz_pt).values.any()
        if bcheck:
            logger.error("Null values found in grid dat file!")
            raise Exception("Null values found in grid dat file!")


        ### Any mesh read from a CDAT or DF will be a point cloud, since
        ### these file formats can't store connectivity
        self.to_point_cloud( keep_nodes = True, replace_comp = True )

        return
    
    def read_tecplot_ascii(self,file):
        """ Read tecplot ASCII file.
        
        Args:
            file (:obj:`str`): Path to file.
        """
        #this is temporary until we can get up and running with libmesh
        nodes = {}
        elem  = {}
        comp  = {}
        tecfile = open (file, 'r')
        varline = tecfile.readline ()
        varline = varline.replace ('"', ' ')
        vars = varline.split ('=')[1].split ()
        for cv in vars:
            nodes [cv] = []
        zoneline = tecfile.readline ()
        npoints = int (zoneline.split (',') [0].split () [-1])
        nele    = int (zoneline.split (',') [1].split () [-1])
        for cn in range (npoints):
            nodeline = tecfile.readline ()
            nodes [vars [0]].append (float (nodeline.split ()[0]))
            nodes [vars [1]].append (float (nodeline.split ()[1]))
            nodes [vars [2]].append (float (nodeline.split ()[2]))
        self.xyz_pt = pd.DataFrame.from_dict (nodes)

        for ce in range (nele):
            eleline = tecfile.readline ()
            comp[ce] = 1
            elem[ce] = []
            for cnode in eleline.split ():
                elem [ce].append (int (cnode))
        self.conn = pd.DataFrame.from_dict(elem, 'index')
        self.comp= pd.DataFrame.from_dict(comp, 'index')
        return
            
    def read_libmesh(self,file):           
        """ Use libmesh to read file.
        
        Args:
            file (:obj:`str`): Path to file.
        """
        try:
            tmp = steam.libmesh.Mesh(file)
        except:
            return

        nodes, elems, comps = tmp.to_python()

        self.xyz_pt = pd.DataFrame.from_dict (nodes)
        self.conn    = pd.DataFrame.from_dict (elems, 'index')
        self.comp   = pd.DataFrame.from_dict (comps, 'index')
        return



    def write_hdf5(self,hdf5,root="/"):
        """ Write the mesh in memory to disk.
        
        This will overwrite any data that is in the HDF5 file.

        Args:
            hdf5 (:obj:`~pandas.io.pytables.HDFStore`): HDFStore object.
            root (:obj:`str`,optional): Location in HDF5 to write database.
                Default is to write to root location "/".
        """

        ### The mesh is made queryable so that we can select points
        ### by their mesh coordinates.

        hdf5.put(
                 root+"/xyz_pt",
                 self.xyz_pt,
                 format='t',
                 data_columns = True
                 )
        if (self.xyz_el is not None):
            hdf5.put(
                     root+"/xyz_el",
                     self.xyz_el,
                     format='t',
                     data_columns = True
                     )
        hdf5.put(
                 root+"/ele",
                 self.conn,
                 format='t',
                 )
        hdf5.put(
                 root+"/comp",
                 self.comp,
                 format='t',
                 )

        ### Have to stick the attributes somewhere
        iroot = root+'/xyz_pt/table'
        hdf5.get_storer(iroot).attrs.info           = self.info
        #hdf5.get_storer(iroot).attrs.comps          = self.comp_table
        hdf5.get_storer(iroot).attrs.is_point_cloud = self.is_point_cloud
        hdf5.get_storer(iroot).attrs.hash           = self.hash

        ### comp_table is stored as a pickle object allowing it to be
        ### arbitrarily large
        steam.util.write_hdf5_pnode( hdf5, f'{root}/comp_table',
                                     self.comp_table )

        ### Write the static solution
        static_root = f'{root}/static_soln'
        if self.static_soln is not None:
            self.static_soln.write_hdf5(
                                   hdf5,
                                   root=static_root,
                                   options={"ioMesh":False}
                                   )
        else:
            ### Ensure that there is no static_soln
            try:
                hdf5.remove( static_root )
            except KeyError:
                pass
            #raise NotImplementedError( 'Not done here yet.' )

        return

    def read_hdf5(self,hdf5,root="/"):
        """ Read the mesh from the HDF5 file to memory.

        Args:
            hdf5 (`~pandas.io.pytables.HDFStore`): HDFStore object.
            root (:obj:`str`,optional): Location in HDF5 to write database.
               Default is to write to root location "/".
       """


        ### Check to make sure that this database has a mesh
        exist = hdf5.get_node(root)
        if exist is None:
            print(" -- No mesh found in HDF5 File")
            return None

        self.xyz_pt = hdf5.get(root+"/xyz_pt")
        try:
            self.xyz_el = hdf5.get(root+"/xyz_el")
        except:
            self.xyz_el = None

        self.conn    = hdf5.get(root+"/ele")
        self.comp   = hdf5.get(root+"/comp")

        ### Where I stuck the attributes
        iroot = root+'/xyz_pt/table' 
        try:
            self.comp_table = steam.util.read_hdf5_pnode( hdf5, 
                                                      f'{root}/comp_table' )
        except TypeError:
            ### Ensure backwards compatibility with old way 
            ### of storing comp_table
            self.comp_table = hdf5.get_storer(iroot).attrs.comps
        self.info       = hdf5.get_storer(iroot).attrs.info

        #! Meshes in databases have hashes
        try:
            self.hash       = hdf5.get_storer(iroot).attrs.hash
        except AttributeError:
            pass

        #! If I have a is_point_cloud attribute, read it.  Otherwise default false
        try:
            self.is_point_cloud = hdf5.get_storer(iroot).attrs.is_point_cloud
        except AttributeError:
            self.is_point_cloud = False


        ### Check whether there is a static_soln and read it if there is one
        static_root = root + '/static_soln'
        if hdf5.get_node( static_root ) is not None:

            ### Create a solution if it doesn't yet exist
            if self.static_soln is None:
                self.static_soln = steam.solution.Solution()

            self.static_soln.read_hdf5(
                        hdf5,
                        root = static_root,
                        options={"ioMesh":False,"readData":True}
                        )
            ### Manually set static_soln.mesh; otherwise we'd get in an
            ### infinite loop of reading meshes and static solutions
            self.static_soln.mesh = self

        return self

    #! Component manipulation

    def set_comp(self,name,comp_list):
        """ Create a new component entry.
        
        The comp_list can contain other components (as strings) or element 
        id numbers (as integers).  If you're refering to one of the 
        component ids tracked in libMesh (an integer), then you
        need to pass it as a string.  For example, "'1'" is component 
        number 1, "1" is element 1.

        self.comp_table includes two items comp, and elem.  'comp' list 
        other components that need expanded.  'elem' list element ids.  
        Intent is to make their use faster.

        Args:
            name (:obj:`str`): The compment name.
            comp_list (:obj:`list`): List of component numbers or names to assign to component.
        """
        
        subcomp = list()
        elemids = list()
        for item in comp_list:
            if type(item) == type("string"):
                subcomp.extend([item])
           #elif type(item) == type(1):
            elif (int(item) == item):
                elemids.extend([item])
            else:
                print("Type:",type(item))
                raise TypeError("Component type not found for {} !".format(
                                                                        item))

        self.comp_table[name] = dict()
        self.comp_table[name]['comp'] = subcomp
        self.comp_table[name]['elem'] = elemids

    def set_comp_by_xyz(self,name, xmin=None, xmax=None, ymin=None, ymax=None, zmin=None, zmax=None):
        """ Create a new component entry based on xyz min/max values.
        
        Args:
            name (:obj:`str`): The name of the new compment.
            xmin (:obj:`float`):  The x coordinate minimum value for the bounding box
            xmax (:obj:`float`):  The x coordinate maximum value for the bounding box
            ymin (:obj:`float`):  The y coordinate minimum value for the bounding box
            ymax (:obj:`float`):  The y coordinate maximum value for the bounding box
            zmin (:obj:`float`):  The z coordinate minimum value for the bounding box
            zmax (:obj:`float`):  The z coordinate maximum value for the bounding box
        """
        elemids = list()
        # If the element coordinates have not been set then calculate them
        if self.xyz_el is None: self.get_xyz_el()

        for i in range (self.xyz_el.shape[0]):
            inside_bbox = True

            try:
                if self.xyz_el['X'][i] < xmin: inside_bbox = False
            except TypeError:
                pass

            try:
                if self.xyz_el['X'][i] > xmax: inside_bbox = False
            except TypeError:
                pass

            try:
                if self.xyz_el['Y'][i] < ymin: inside_bbox = False
            except TypeError:
                pass

            try:
                if self.xyz_el['Y'][i] > ymax: inside_bbox = False
            except TypeError:
                pass

            try:
                if self.xyz_el['Z'][i] < zmin: inside_bbox = False
            except TypeError:
                pass

            try:
                if self.xyz_el['Z'][i] > zmax: inside_bbox = False
            except TypeError:
                pass

            if inside_bbox: elemids.append(i)

        self.comp_table[name] = dict()
        self.comp_table[name]['comp'] = []
        self.comp_table[name]['elem'] = elemids

    def set_comp_by_x_and_theta( *args, **kwargs ):
        """ NOTE: This method is very poorly documented and 
        seems to fill a niche role that will be needed by almost no one.  
        As a result, Cam Ebner edited it to
        remove all functionality on 3/31/21.  The only thing it should
        ever do is raise a NotImplementedError when called.

        Original docstring:
        Create a new component entry based on x min/max and theta min/max.
        
        Args:
            name (:obj:`str`): The name of the new compment.
            xmin (:obj:`float`):  The x coordinate minimum value for the bounding box
            xmax (:obj:`float`):  The x coordinate maximum value for the bounding box
            tmin (:obj:`float`):  The circumferencial direction minimum
            tmax (:obj:`float`):  The circumferencial direction maximum
        """
        raise NotImplementedError( 'This function is constructed assuming a '
                + 'left-handed, undocumented definition of Theta.  Because '
                + 'it is lacking both documentation and wise coding '
                + 'practices, it has been rendered inert.' )

        import math

        pi = 3.14159265358979323846
        elemids = list()
        # If the element coordinates have not been set then calculate them
        if self.xyz_el is None: self.get_xyz_el()

        for i in range (self.xyz_el.shape[0]):
            inside_bbox = True

            x = self.xyz_el['X'][i]
            y = self.xyz_el['Y'][i]
            z = self.xyz_el['Z'][i]

            try:
                if x < xmin: inside_bbox = False
            except TypeError:
                pass

            try:
                if x > xmax: inside_bbox = False
            except TypeError:
                pass

            if y == 0.0 and z > 0.0: theta = 0
            elif y == 0.0 and z < 0.0: theta = 180.0
            elif z == 0.0 and y > 0.0: theta = 90.0
            elif z == 0.0 and y < 0.0: theta = 270.0
            elif y > 0.0 and z > 0.0: theta = math.atan (y/z) * 180.0 / pi
            elif y > 0.0 and z < 0.0: theta = math.atan (y/z) * 180.0 / pi + 180.0
            elif y < 0.0 and z < 0.0: theta = math.atan (y/z) * 180.0 / pi + 180.0
            elif y < 0.0 and z > 0.0: theta = math.atan (y/z) * 180.0 / pi + 360.0


            if tmin < tmax:
                try:
                    if theta < tmin or theta > tmax: inside_bbox = False
                except TypeError:
                    pass

            else:
                try:
                    if theta < tmin and theta > tmax: inside_bbox = False
                except TypeError:
                    pass

            if inside_bbox: elemids.append(i)

        self.comp_table[name] = dict()
        self.comp_table[name]['comp'] = []
        self.comp_table[name]['elem'] = elemids

    def get_comp(self,names):
        """ Return a list of elements associated with component(s).

        List of returned elements is made unique.  Sub-components are 
        expanded and added to the list.  If a :obj:`~steam.mesh.Mesh` is 
        desired, rather than a list of elements, use 
        :obj:`~steam.mesh.Mesh.return_comp`.

        Args:
            names (:obj:`list`): List of compment names.  'all' returns everything.

        Returns:
            (:obj:`list`): List of element ids to assigned to component.
        """

        names_list = names

        if type(names) == type(str()):
            names_list = [names]
        if type(names) != type(list()):
            names_list = [names]

        out_list = list()
        for name in names_list:

            if (name == 'all'):
                return list(self.comp.index)

            # Integers are OK since they will map to a componet in self.comp
            try:
                i = int(name)
            except:
                pass
            else:
                subset   = self.comp[self.comp[0] == i]
                out_list.extend(list(subset.index))
                continue

            #! Make sure that the name exists
            try:
                self.comp_table[name]
            except:
                raise Exception("Component {} not found in component table and not int!".format(name))

            #! Going to run through the components associated with name and
            #! append their elements to the master list that we are assembling.
            for comp in self.comp_table[name]['comp']:
                ids = self.get_comp(comp)
                out_list.extend(ids)

            #! Append all element ids, too
            out_list.extend(self.comp_table[name]['elem'])

        #! Now that we have them all, make unique:
        out_list = list(set(out_list))

        return out_list

    #! Grid Manipulation

    def get_xyz_el(self):
        """ Populate element centroids from node centroids.

            This merely averages the data at the vertices to calcluate element centroids.

        """

        #! Look at each element of the grid and average the centroid
        #! at each vertex.
        
        vars = ["X","Y","Z"]
        elements  = self.conn.values;
        node_mat  = self.xyz_pt[vars].values;
        elem_mat  = np.zeros([self.conn.shape[0],self.xyz_pt.shape[1]])

        for (e,elem) in enumerate(elements):
            for node in elem:
                elem_mat[e] += node_mat[node]
            elem_mat[e] /= float(elem.size)

        #! Replace the data
        elem_df = pd.DataFrame( elem_mat, columns=vars, index=self.conn.index )
#        elem_df = pd.DataFrame( elem_mat, columns=vars )

        self.xyz_el = elem_df

        return

    def to_point_cloud(self, keep_nodes = True, replace_comp = False):
        """ This will convert the mesh to a point cloud.

        Connectivity is thrown out and replaced with a one-to-one with either 
        the nodes or the element centroids, dependinging on 'keep_nodes'.

        Components are all set to 1 if 'replace_comp' is True OR if self.conn
        is not defined.

        Args:
            keep_nodes (:obj:`bool`): Optional flag indicating whether to save the mesh's node locations or element centroids.  Defaults to True (saving node locations and throwing out element centroids).
            replace_comp (:obj:`bool`): Optional flag indicating whether or not to overwrite the component tags and comp_table.  Defaults to False.
        """
        
        logger.debug( 'Converting to point cloud' )

        if keep_nodes:
            ### Replace components if necessary
            try:
                if self.xyz_pt.shape[0] != self.xyz_el.shape[0]:
                    replace_comp = True
            except AttributeError:
                replace_comp = True

            self.xyz_el = self.xyz_pt
        else:
            self.xyz_pt = self.xyz_el
 
        #! Simple one-to-one connectivity
        npt         = self.xyz_el.values.shape[0]
        conn        = np.arange(npt)
        self.conn   = pd.DataFrame(conn)

        if replace_comp or (self.comp is None):
            logger.info( 'Overwriting/creating mesh components' )
            comp            = np.ones(npt,dtype='int')
            self.comp       = pd.DataFrame(comp)
            self.comp_table = dict()

        self.is_point_cloud = True

    def transform(self, ops, comps=[]):
        """ Mirror/Rotate/Translate grid. 
        
        Rotate about x, y, or z.  "r{DEG}{AXIS}": "r90x","x-90y"
        Translate in x, y, or z.  "t{DIST}{DIR}": "t-103.483x"
        Scale in x, y, or z.      "s{FACT}{AXIS}": "s0.0254x","s0.0254y"
        Mirror across plane of constant x, y, or z. "m{LOC}{AXIS}": "m0y", "m10z"

        Args:
            ops (:obj:`list`): List of operations to perform (in order!).
            comps (:obj:`list`): List of components to transform.

        """

        #! Set-up the Regex to match the operations
        import re
        parse = re.compile("([rtmsRTMS])([-0-9.]+)([xyzXYZ])")
       
        # Check if component list was passed
        # Transform entire mesh if no component list found
        # If component list was passed pull list of vertices from components 
        node_mat  = self.xyz_pt.values;
        if len(comps) == 0:
            i_verts = list(self.xyz_pt.index.values)
        else:
            i_verts = list(set(self.conn.values[self.get_comp(comps)].flatten()))

        #! Go through each operation and execute it in order
        for op in ops:
            #! Parse the operation
            match = parse.match(op)
            if (match is None):
                raise Exception("Operation '{}' was not recognized!".format(op))
            (type,val,axis) = match.groups()

            #print(type,val,axis)
            val = float(val)

            if (type.lower() == "t"):
                print("Translating in {} by {} units".format(axis,val))
                if (axis.lower() == "x"):
                    node_mat[i_verts,0] += val
                if (axis.lower() == "y"):
                    node_mat[i_verts,1] += val
                if (axis.lower() == "z"):
                    node_mat[i_verts,2] += val

            if (type == "r"):
                print("Rotating about {} by {} degrees".format(axis,val))
                rad = val/180.0 * np.pi
                if (axis.lower() == "x"):
                    x = node_mat[i_verts,0]
                    y = node_mat[i_verts,1]*np.cos(rad) - node_mat[i_verts,2]*np.sin(rad)
                    z = node_mat[i_verts,2]*np.cos(rad) + node_mat[i_verts,1]*np.sin(rad)
                if (axis.lower() == "y"):
                    x = node_mat[i_verts,0]*np.cos(rad) + node_mat[i_verts,2]*np.sin(rad)
                    y = node_mat[i_verts,1]
                    z = node_mat[i_verts,2]*np.cos(rad) - node_mat[i_verts,0]*np.sin(rad)
                if (axis.lower() == "z"):
                    x = node_mat[i_verts,0]*np.cos(rad) - node_mat[i_verts,1]*np.sin(rad)
                    y = node_mat[i_verts,1]*np.cos(rad) + node_mat[i_verts,0]*np.sin(rad)
                    z = node_mat[i_verts,2]
                node_mat[i_verts,0] = x
                node_mat[i_verts,1] = y
                node_mat[i_verts,2] = z

            if (type == "m"):
                print("Mirroring across {} = {} ".format(axis,val))
                if (axis.lower() == "x"):
                    node_mat[i_verts,0] -= val
                    node_mat[i_verts,0] *= -1.0
                    node_mat[i_verts,0] += val
                if (axis.lower() == "y"):
                    node_mat[i_verts,1] -= val
                    node_mat[i_verts,1] *= -1.0
                    node_mat[i_verts,1] += val
                if (axis.lower() == "z"):
                    node_mat[i_verts,2] -= val
                    node_mat[i_verts,2] *= -1.0
                    node_mat[i_verts,2] += val

            if (type.lower() == "s"):
                print("Scaling in {} = {} ".format(axis,val))
                if (axis.lower() == "x"):
                    node_mat[i_verts,0] *= val
                if (axis.lower() == "y"):
                    node_mat[i_verts,1] *= val
                if (axis.lower() == "z"):
                    node_mat[i_verts,2] *= val

        #! Update xyz_wl if I have it
        if (self.xyz_el is not None):
            self.get_xyz_el()

        return

    @steam.util.timer(logger)
    def slice(self,planes,component=None,delete=None,
              returnList=False,progress=False):
        """ Slice a mesh with cutting planes.

        planes are described as a list with three three-part tuples.
        Each tuple is a point on the plane.  For example:
        x = 45 plane:  [(45,0,0),(45,1,0),(45,0,1)]
        The arguement is a list of planes, so it should be a list
        of lists or tuples.

        Components can be replaced for each plane.  This needs to be
        a list equal in length to the planes list.  Each entry is a 
        tuple of the new component number in the (pos,neg) direction.
        -1 means do not change the components.  To change the component
        to 1 for everything in the positive normal direction but to not
        change those in the negative normal direction: [(1,-1)]

        Geometry can be deleted by specifying a delete list.  This is
        the same size as the planes list.  Each entry is a tuple to
        either keep (-1) or delete (1) the faces on the (pos,neg)
        direction.

        The returnList (optional) is the original ID of the elements
        that were cut.  This is needed when doing line loads to relate
        the new elements to their originals.  It is size...?

        Args:
            planes: list of planes, details above.
            component: Optional list of component changes.
            delete: Optional list of delete commands.
            returnList: Return the list of original elements? [F]
            progress: Show progress bar? [F]

        Returns:
            list of original elements for each element (Optional)
        """

        if self.is_point_cloud:
            raise AttributeError("PointCloud objects cannot be sliced!")
        #assert self.is_point_cloud == False, (
        #       'This mesh is a point cloud -- Should not get here' )

        debug = False # Write out a bunch of sub-timings?

        import copy
        if returnList:
            #! - copy since we overwrite some of new_conn
            originalElem = copy.deepcopy(self.conn.index.values)

        if component is None:
            component = []
        if delete is None:
            delete = []

        #! Fill component array if not enough provided
        delta = len(planes) - len(component)
        for i in range(delta):
            component.append((-1,-1))

        #! Fill delete array if not enough provided
        delta = len(planes) - len(delete)
        for i in range(delta):
            delete.append((-1,-1))
        
        #! If I actually want to delete something, then what I'm going to do is
        #! mark it as a component number past what exists or is being assigned.
        delcomp_min = self.comp.max()[0] + 1001
        if (len(component)>0):
            delcomp_min = max(delcomp_min,max(max(component)))
        delcomp_max = delcomp_min
        #print(delcomp_min,delcomp_max)
        for (i,d1) in enumerate(delete):
            values = [component[i][0],component[i][1]]
            for (j,d2) in enumerate(d1):
                #! We delete if d2 is greater than zero
                if (d2 > 0):
                    values[j] = delcomp_max
                    delcomp_max += 1
            component[i] = (values[0],values[1])

        #! Format for cleaner output
        str_fmt = "({:.5e}, {:.5e}, {:.5e})"
    
        #! Do the work
        for (ip,plane) in enumerate(planes):
            if (progress):
                steam.util.progress_bar(ip+1,len(planes),'   Cut Planes:')

            #! Check on the size of things
            if (len(plane) != 3):
                raise Exception("Planes must have three points!")

            for p in plane:
                if (len(p) != 3):
                    raise Exception("Points must have three components!")

            str_plane  = "[" +str_fmt.format(*plane[0])
            str_plane += ", "+str_fmt.format(*plane[1])
            str_plane += ", "+str_fmt.format(*plane[2])+"]"

            logger.debug(" Slicing with plane:")
            logger.debug("  ["+str_fmt.format(*plane[0])+",")
            logger.debug("   "+str_fmt.format(*plane[1])+",")
            logger.debug("   "+str_fmt.format(*plane[2])+"]")
            import time
            start = time.time()

            #! Get the plane verticies and the normal (to see what side things are on
            pv1 = np.array(plane[0])
            pv2 = np.array(plane[1])
            pv3 = np.array(plane[2])
            norm = np.cross(pv2-pv1,pv3-pv2)
            mag  = np.dot(norm,norm)
            #print(" - Plane normal is :",norm)
            #print(" - Norm magnitude  :",mag)
            norm = norm / mag
            logger.debug("  Unit normal is  :")
            logger.debug("   "+str_fmt.format(*norm))
            
            if debug:
                logger.debug(" -- Init    : {}".format(time.time()-start))
            start = time.time()
            #! Loop through the triangles and see which are cut by this plane
            node_xyz  = self.xyz_pt.values;
            tri_conn  = self.conn.values;
            if debug:
                logger.debug(" -- Pandas  : {}".format(time.time()-start))
            
            start = time.time()
            #! Find each node's distance to the plane
            #! - Numpy can do this all in the array.  No loops required
            node_dist = np.dot(node_xyz - pv1,norm)
            if debug:
                logger.debug(" -- NodeDist: {}".format(time.time()-start))
 
            #! Build tri array of node distances from plane
            #! - Numpy is smart enough to use a 2-D array (conn)
            #!   to describe the indexing in node_dist and return
            #!   the appropriate values into a 2-D array (tri_dist)
            start = time.time()
            tri_dist = node_dist[tri_conn]
            if debug:
                logger.debug(" -- TriDist : {}".format(time.time()-start))

            start = time.time()
            #! dpos/dneg contain True/False based on if a triangle has
            #!  a node that is either less than zero or greater than zero
            dneg  = tri_dist.min(axis=1) < 0
            dpos  = tri_dist.max(axis=1) > 0
            #! Make this a 2D array where each row tells us if a tri has
            #!  a node on the negative side and positive side
            #dbool = np.column_stack((dneg,dpos))
            #! If a node is split, then all elements in the row are True.  Map to new list
            #split = np.all(dbool,axis=1)

            split = np.all([dneg,dpos],axis=0)
            #split = list(map(all,dbool))
            #! This is a handy tool.  compress returns the indicies for values that are True.
            #!  These are the triangles that are split by the plane.
            #from itertools import compress
            #cut_tri = list(compress(range(len(split)),split))

            cut_tri = np.nonzero(split)[0]
            if debug:
                logger.debug(" -- Get Cuts: {}".format(time.time()-start))

        #   print(" -- - t1    :",tmp_time[0])
        #   print(" -- - t2    :",tmp_time[1])
            
            if (len(cut_tri) == 0):
                if debug:
                    logger.debug(" - No intersected triangles")
            #import code

            ### There is a way to combine the "do I cut this edge?" and "where to cut this edge"
            ### in one step.  Hopefully all in vector math, too.  Just find the intersection of
            ### the extrapolated lines and the plane.  If the distance that I just interpolated
            ### is outside of the the bounds of the actual edge, then I'm not cut and I don't need
            ### the point I just made.  If I need it, then I put an entry in the database if one
            ### wasn't there already. That last bit is hard to vectorize.  Maybe run a collapse step 
            ### later with a kd-tree where I remove duplicates.


            #print("Going to cut: ",len(cut_tri))
            
            #! Allocate temporary matricies (?)  Need to:
            #! - Cut each edge and build a mapping of those cuts (2*num_cut edges will be cut)
            #! - Make new triangles (3*num_cut new, but one replaces)
            
            start = time.time()
            #! This will house what edges get cut
            
            #! Cut the edges I need to cut
            cut_conn = tri_conn[cut_tri]
            cut_dist = tri_dist[cut_tri]
            cut_nxyz = node_xyz[cut_conn]

            # I think that this can be done by looking at the length of the edges
            # and comparing that to the distances from the plane...

            # #return(cut_conn,cut_dist,cut_nxyz,norm)
            # ##(conn,dist,nxyz,norm) = mesh.slice(planes,comps)
            with np.errstate(divide='ignore', invalid='ignore'):
              # Cut them all
              xyz0 = cut_nxyz[:,0]
              xyz1 = cut_nxyz[:,1]
              xyz2 = cut_nxyz[:,2]
              dis0 = cut_dist[:,0]
              dis1 = cut_dist[:,1]
              dis2 = cut_dist[:,2]

              #! Edge lengths in plane normal direction
              len01 = dis1 - dis0
              len12 = dis2 - dis1
              len20 = dis0 - dis2

              #! This all works, but it has issues if a line is parallel to the plane.
              #! Then lenXX is 0.0, so you end up dividing by zero.

              #! Location for intersects:
              #print(np.hstack((len01,len12,len20)).reshape(-1,3))
              int01 = xyz0 - (dis0/len01).reshape(-1,1) * (xyz1 - xyz0)
              int12 = xyz1 - (dis1/len12).reshape(-1,1) * (xyz2 - xyz1)
              int20 = xyz2 - (dis2/len20).reshape(-1,1) * (xyz0 - xyz2)

              #! Is this inside of the edge?
              bound0 = np.dot(int01-xyz0,norm) / len01
              bound1 = np.dot(int12-xyz1,norm) / len12
              bound2 = np.dot(int20-xyz2,norm) / len20

              #! Check to make sure it is inside of the edge and make a T/F array
              ecut0 = np.column_stack(( bound0 > 0 , bound0 < 1.0 )).all(axis=1)
              ecut1 = np.column_stack(( bound1 > 0 , bound1 < 1.0 )).all(axis=1)
              ecut2 = np.column_stack(( bound2 > 0 , bound2 < 1.0 )).all(axis=1)

              ecut = np.column_stack((ecut0,ecut1,ecut2))
              xcut = np.hstack((int01,int12,int20)).reshape(-1,3,3)

            if debug:
                logger.debug(" -- PreEdges: {}".format(time.time()-start))

            # So no I have ecut (if an edge is cut) and xcut (the new node along that cut)

            #! Add all of the new points where needed

            start = time.time()

            #! Some of these are duplicates, so get that figured out
            #! Get all of the new stuff (where we are cutting)
            all_new  = xcut[ecut]
            idx      = []
            for xyz in all_new:
                j=np.where(np.all((all_new - xyz) == 0.0,axis=1))[0][0]
                idx.append(j)
            idx = np.array(idx,dtype=int)
            
            idmap = np.array(idx) - np.array(range(len(idx)),dtype=int)
            idmap[idmap==0] = 1  # These are new nodes I       want
            idmap[idmap< 0] = 0  # These are new nodes I don't want
            
            #! Update idx based on the masking in idmap
            new_idx = np.cumsum(idmap) - 1
            #new_idx[idx[idmap==0]]
            new_idx = new_idx[idx]

            #! Indexes to new points
            old_node = node_xyz.shape[0]
            icut = np.zeros((xcut.shape[0],xcut.shape[1]),dtype=int) - 1
            #icut[ecut] = np.cumsum(idmap) + old_node - 1
            #icut[ecut] = idx + old_node - 1
            icut[ecut] = new_idx + old_node
            
            #! Add the new nodes
            node_xyz      = np.append(node_xyz,all_new[idx[idmap==1]],axis=0)

            if debug:
                logger.debug(" -- Edges   : {}".format(time.time()-start))

            start = time.time()
            #! Ok, fix connectivity and make new triangles if edges are cut
            #! - copy since we overwrite some of new_conn
            import copy
            new_conn = copy.deepcopy(self.conn.values)
            new_comp = copy.deepcopy(self.comp.values)
            
            # Make a conservative guess - every cut tri is cut into three (two new)
            new_tri = np.zeros((len(cut_tri)*2,3),dtype='uint32')
            new_com = np.zeros((len(cut_tri)*2,1),dtype='uint32')
            inew = 0

            #for (t,conn) in enumerate(tri_conn):
            for (it,t) in enumerate(cut_tri):
                node = tri_conn[t]
                ebool = ecut[it]   # Bools about if an edge is cut
                enode = icut[it]   # Nodes for the cut ones
                #for (e,cut) in enumerate(ebool):
                #    # Find the first instance of this point in the nodes_array
                #    if (cut):
                #        n = np.all(xcut[it,e] == node_xyz,axis=1).argmax()
                #        enode[e] = n
            
                if (np.all(enode == [-1,-1,-1])):
                    #! No cuts
                    continue
                elif (enode[1] == -1 and enode[2] == -1):
                #
                #    1 ----a---- 2
                #    |    /     /
                #    |   /    /
                #    |  /   /
                #    |  /  /
                #    | / /
                #    | /
                #    3
                    a = enode[0]
            
                    new_conn[t][1] = a

                    new_tri[inew][0] = a
                    new_tri[inew][1] = node[1]
                    new_tri[inew][2] = node[2]
                    new_com[inew] = new_comp[t]
                    inew += 1

                        
                elif (enode[0] == -1 and enode[2] == -1):
                #
                #    1 --------- 2
                #    |\         /
                #    |  \     /
                #    |    \ /
                #    |     a
                #    |   /
                #    | /
                #    3
                    a = enode[1]
            
                    new_conn[t][2] = a

                    new_tri[inew][0] = node[0]
                    new_tri[inew][1] = a
                    new_tri[inew][2] = node[2]
                    new_com[inew] = new_comp[t]
                    inew += 1
                    
                elif (enode[0] == -1 and enode[1] == -1):
                #
                #    1 --------- 2
                #    |  ______/ /
                #    | /      /
                #    a      /
                #    |     /
                #    |   /
                #    | /
                #    3
                    a = enode[2]
            
                    new_conn[t][2] = a
                    new_tri[inew][0] = a
                    new_tri[inew][1] = node[1]
                    new_tri[inew][2] = node[2]
                    new_com[inew] = new_comp[t]
                    inew += 1
                    
                elif (enode[2] == -1):
                #
                #    1 ----a---- 2
                #    |\    |    /
                #    |  \  |  /
                #    |    \|/
                #    |     b
                #    |   /
                #    | /
                #    3
                    a = enode[0]
                    b = enode[1]
            
                    #print("Cut 1",node,a,b)
                    new_conn[t][1] = a
                    new_conn[t][2] = b

                    new_tri[inew][0] = a
                    new_tri[inew][1] = node[1]
                    new_tri[inew][2] = b
                    new_com[inew] = new_comp[t]
                    inew += 1

                    new_tri[inew][0] = b
                    new_tri[inew][1] = node[2]
                    new_tri[inew][2] = node[0]
                    new_com[inew] = new_comp[t]
                    inew += 1

                elif (enode[1] == -1):
                #
                #    1 ----a----- 2
                #    |   //      /
                #    | / |     /
                #    b  /    /
                #    |  |  /
                #    | / /
                #    | /
                #    3
                    a = enode[0]
                    b = enode[2]
            
                    #print("Cut 2",node,a,b)
                    new_conn[t][1] = a
                    new_conn[t][2] = b

                    new_tri[inew][0] = a
                    new_tri[inew][1] = node[1]
                    new_tri[inew][2] = node[2]
                    new_com[inew] = new_comp[t]
                    inew += 1

                    new_tri[inew][0] = b
                    new_tri[inew][1] = a
                    new_tri[inew][2] = node[2]
                    new_com[inew] = new_comp[t]
                    inew += 1

                elif (enode[0] == -1):
                #
                #    1 --------- 2
                #    |\         /
                #    |  \     /
                #    |    \ /
                #    b-----a
                #    |   /
                #    | /
                #    3
                    a = enode[1]
                    b = enode[2]
            
                    #print("Cut 3",node,a,b)
                    new_conn[t][1] = a
                    new_conn[t][2] = b

                    new_tri[inew][0] = node[0]
                    new_tri[inew][1] = node[1]
                    new_tri[inew][2] = a
                    new_com[inew] = new_comp[t]
                    inew += 1

                    new_tri[inew][0] = a
                    new_tri[inew][1] = node[2]
                    new_tri[inew][2] = b
                    new_com[inew] = new_comp[t]
                    inew += 1

                else:
                    node = tri_conn[t]
                    nxyz = []
                    ndis = []
                    for i in [0,1,2]:
                        n = node[i]
                        nxyz.append(node_xyz[n])
                        ndis.append(np.dot(nxyz[i] - pv1,norm))
                    print(ndis)
                    raise Exception("I found an impossible cut:",enode)

            if returnList:
                anew = np.sum(ecut,axis=1)  # The number of cut edges in each
                                            # tri is the same as the number of
                                            # new trianglular elements
                new_ori = originalElem[np.repeat(cut_tri,anew)]
                #for i in range(inew):
                #    new_ori[i] = originalElem[t]
                originalElem = np.append(originalElem,new_ori,axis=0)

            if debug:
                logger.debug(" -- NewTris : {}".format(time.time()-start))
            start = time.time()
            new_conn = np.append(new_conn,new_tri[0:inew],axis=0)
            new_comp = np.append(new_comp,new_com[0:inew],axis=0)


            #Put it back into place
            conn = pd.DataFrame(new_conn)
            xyz  = pd.DataFrame(node_xyz,columns=self.xyz_pt.columns)
            comp = pd.DataFrame(new_comp)
            if debug:
                logger.debug(" -- ToPandas: {}".format(time.time()-start))
            
            #print("Pre:",self.conn.shape)
            #print("Pre:",self.comp.shape)
            #print("Pre:",self.xyz_pt.shape)
            self.conn   = conn
            self.xyz_pt = xyz
            self.comp   = comp

            if (self.xyz_el is not None):
                self.get_xyz_el()
            #print("Post:",self.conn.shape)
            #print("Post:",self.comp.shape)
            #print("Post:",self.xyz_pt.shape)

            start = time.time()
            #! Check to see if we need to update components
            if max(component[ip]) > 0:
                (cpos,cneg) = component[ip]

                #! Loop through every element.  At this point it will be
                #! on one side or the other.
                node_xyz  = self.xyz_pt.values
                tri_conn  = self.conn.values
                new_comp  = self.comp.values

                node_dist = np.dot(node_xyz - pv1,norm)
                tri_dist = node_dist[tri_conn]

                #! dpos/dneg contain True/False based on if a triangle has
                #!  a node that is either less than zero or greater than zero (with some tolerance)
                if (cpos > 0):
                    #dpos  = tri_dist.max(axis=1) >  1e-5
                    #tris = list(compress(range(len(dpos)),dpos))
                    tris = np.where(tri_dist.max(axis=1)>1e-5)[0]
                    new_comp[tris] = cpos

                if (cneg > 0):
                    #dneg  = tri_dist.min(axis=1) < -1e-5
                    #tris = list(compress(range(len(dneg)),dneg))
                    tris = np.where(tri_dist.max(axis=1)<1e-5)[0]
                    new_comp[tris] = cneg

                #for (t,conn) in enumerate(tri_conn):
                #    side = []
                #    for (n,node) in enumerate(conn):
                #        n = node_xyz[node]
                #        dist = n - pv1
                #        dist = np.dot(dist,norm)
                #        side.append(dist)
                #    if   (cpos > 0 and max(side) >  1e-5):
                #        new_comp[t] = cpos
                #    elif (cneg > 0 and min(side) < -1e-5):
                #        new_comp[t] = cneg
            if debug:
                logger.debug(" -- Comps   : {}".format(time.time()-start))

            #! Check if we need to delete components
            if max(delete[ip]) > 0:
                logger.debug(" -- Deleting sides...")
                comps = []
                for (i,side) in enumerate(delete[ip]):
                    if (side >= 0):
                        comps.append(component[ip][i])
                #print("Remove comp: ",comps)
                self.rm_comp(comps)

#        if (self.xyz_el is not None):
#            self.get_xyz_el()

        if returnList:
            return originalElem
        else:
            return

        return

    def rm_comp(self,comp_list=[]):
        """Remove components from the mesh.

        This will remove the elements that are part of a certain
        component or components from the mesh.  Indexing for the
        elements is not changed in order to maintain compatability
        between the original grid elements.  Node indexing is
        reset to be contiguous.

        Args:
            comp_list : List of components to remove.
        """
    
        elements  = self.conn.values
        el_drop   = []  # Elems to drop
        pt_drop   = []  # Nodes to drop
    
        for comp in comp_list:
            elems = self.get_comp(comp)
            el_drop.extend(elems)
    
        el_drop = list(set(el_drop))
        self.conn = self.conn.drop(el_drop)
        self.comp = self.comp.drop(el_drop)
        if (self.xyz_el is not None):
            self.xyz_el = self.xyz_el.drop(el_drop)
    
        #! Make an array that will tell us if the nodes are used
        #! and also give us their new mapping for them.
        used = np.zeros([self.xyz_pt.shape[0],1])
        elements  = self.conn.values
        for (e,elem) in enumerate(elements):
            for node in elem:
                used[node] = 1
    
        node_map  = np.full( [self.xyz_pt.shape[0],1],-1.0)
        count = 0
        for (n,node) in enumerate(used):
            #! Is this node used?
            if (node == 1):
                node_map[n] = count
                count += 1
            else:
                pt_drop.extend([n])
    
        self.xyz_pt = self.xyz_pt.drop(pt_drop)
        self.xyz_pt.reset_index(drop=True,inplace=True)
    
        #! Renumber the vertex indicies in the conn array
        for (e,elem) in enumerate(elements):
            for (n,node) in enumerate(elem):
                elements[e][n] = node_map[node]

        return

    def return_subset(self,el_list, include_static=True):
        """Return elements from the mesh as a new mesh.

        This will return the elements list as input.
        Indexing for the elements is not changed in order 
        to maintain compatability between the original grid 
        elements.  Node indexing is reset to be contiguous.

        Args:
            el_list (:obj:`list`): List of elements to keep.
            include_static (:obj:`bool`): Flag indicating whether or not a consistent subset of the static_solution should be included as the static solution of the returned mesh

        Returns:
            (:obj:`~steam.mesh.Mesh`): Requested subset of the mesh.
        """

        mesh = steam.mesh.Mesh()
        pt_list   = []  # Nodes to drop

        #! I'm using .index.isin since it does two things:
        #! - sorts the elments and ignores out-of-bounds
        def df_sub(df,list):
#            return df[df.index.isin(list)]
            return df.loc[df.index.isin(list)]

        mesh.conn = df_sub(self.conn,el_list)
        mesh.comp = df_sub(self.comp,el_list)
        import copy
        mesh.comp_table = copy.deepcopy(self.comp_table)
        if (self.xyz_el is not None):
            mesh.xyz_el = df_sub(self.xyz_el,el_list)
    
        #! Make an array that will tell us if the nodes are used
        #! and also give us their new mapping for them.
        used = np.zeros([self.xyz_pt.shape[0],1], dtype = bool )
        elements  = mesh.conn.values
        for (e,elem) in enumerate(elements):
            for node in elem:
                used[node] = True
    
        node_map  = np.full( [self.xyz_pt.shape[0],1],-1.0)
        count = 0
        for (n,node) in enumerate(used):
            #! Is this node used?
            if node:
                pt_list.extend([n])
                node_map[n] = count
                count += 1
    
        mesh.xyz_pt = df_sub(self.xyz_pt,pt_list)
        mesh.xyz_pt.reset_index(drop=True,inplace=True)
    
        #! Renumber the vertex indicies in the conn array
        for (e,elem) in enumerate(elements):
            for (n,node) in enumerate(elem):
                elements[e][n] = node_map[node]

        mesh.is_point_cloud = self.is_point_cloud

        if include_static and self.static_soln is not None:
            mesh.static_soln = self.static_soln.return_subset( el_list )

        return mesh

    def return_comp(self,comp_list, include_static=True ):
        """Return components from the mesh as a new mesh.

        This will return the elements that are part of a certain
        component or components from the mesh.  Indexing for the
        elements is not changed in order to maintain compatability
        between the original grid elements.  Node indexing is
        reset to be contiguous.

        Args:
            comp_list (:obj:`list`): Components to keep.
            include_static (:obj:`bool`): Flag indicating whether or not a consistent subset of the static_solution should be included as the static solution of the returned mesh

        Returns:
            (:obj:`~steam.mesh.Mesh`): Requested subset of the mesh.
        """
    
        mesh = steam.mesh.Mesh()

        el_list   = []  # Elems to get
    
        for comp in comp_list:
            elems = self.get_comp(comp)
            el_list.extend(elems)
    
        el_list = list(set(el_list))

        return self.return_subset(el_list, include_static)

    def get_areas(self,comps=None,wet=False):
        """ Return a mesh's projected area in X,Y,Z directions.
    
        When wet=False, these are signed areas.
        The area is signed consisant with the mesh elements' normal directions.
        The output is a numpy array with [X_area, Y_area, Z_area].

        If wet=True, this is wetted area.
        The output is a numpy aray with [X_wet, Y_wet, Z_wet, Total]
        X,Y,Z_wet is the sum of the magnitudes of the X,Y,Z projected areas.
        Total is the total area (RSS of element X,Y,Z projections).
        
        """

        #! Loop over every triangle and get the normal and area
        if comps is not None:
            elements = self.conn.values[self.get_comp(comps)]
        else:
            elements  = self.conn.values
        nodes     = self.xyz_pt.values

        # If wet, then I'm returning four values (X, Y, Z, Total)
        # otherwise, just three (X, Y, Z):w
        if wet:
            areas     = np.array([0.0,0.0,0.0,0.0])
        else:
            areas     = np.array([0.0,0.0,0.0])
    
        for (e,elem) in enumerate(elements):
            n1 = nodes[elem[0]]
            n2 = nodes[elem[1]]
            n3 = nodes[elem[2]]
    
            v1 = n2 - n1
            v2 = n3 - n2
            v3 = np.cross(v1,v2)
            area = v3 / 2.0

            if wet:
                areas[0:3] += abs(area)
                areas[3]   += np.dot(area,area)**0.5
            else:
                areas += area

        return areas

#! Simple mesh creation utilities

def mk_cube (length, dims, file=0, filename="cube"):
    """ Make a cube mesh.

    The cube is centered on 0,0,0.

    Components are made for each of the six sides.  The mapping is:
        xmin = 1   xmax = 2
        ymin = 3   ymax = 4
        zmin = 5   zmin = 6
    
    Args:
        length (:obj:`int`): cube edge length.
        dims   (:obj:`int`): Number of points along each edge
    Returns:
        :obj:`~steam.mesh.Mesh()`: Mesh object
    """

    nodes,elems,comps = build_cube(length, dims)
    write_simple_grid(nodes,elems,comps,file,filename)
    mesh = simple_to_mesh(nodes,elems,comps)

    return mesh

def mk_sphere (radius, npts, file=0, filename="sphere"):
    """ Make a sphere mesh.

    The sphere is centered on 0,0,0.
    
    Args:
        radius (:obj:`int`): sphere radius
        npts   (:obj:`int`): Number of points along a hemisphere

    Returns:
        :obj:`~steam.mesh.Mesh()`: Mesh object
    """

    nodes,elems,comps = build_sphere(radius, npts)
    write_simple_grid(nodes,elems,comps,file,filename)
    mesh = simple_to_mesh(nodes,elems,comps)

    return mesh

def simple_to_mesh(nodes,elems,comps):
    """ Convert simple grid to Mesh() object.

    Args:
        nodes (:obj:`dict`):  Node X, Y, and Z array
        elems (:obj:`list`):  2D list of element connectivity
        comps (:obj:`list`):  1D list of `int`s to be used as element components

    Returns:
        :obj:`~steam.mesh.Mesh()`: Mesh object
    """

    ### Only `int` inputs are allowed as input components
    if not all( [isinstance(comp, int) for comp in comps] ):
        raise TypeError( "All input components in `comps` must be `int`" )

    if not len(elems) == len(comps):
        raise ValueError( 
                    '"elems" and "comps" inputs must be of equal length.' )

    mesh = Mesh()
    mesh.xyz_pt = pd.DataFrame.from_dict(nodes)
    mesh.conn   = pd.DataFrame( elems )
    mesh.comp   = pd.DataFrame( comps )

    return mesh

def build_cube(length,dims):
    """ Build the cube grid by hand.

    The cube is centered on 0,0,0.
    
    Args:
        length (:obj:`int`): cube edge length.
        dims   (:obj:`int`): Number of points along each edge
        
    Returns:
        :obj:`dict`: Node X, Y, and Z array
        :obj:`list`: 2D list of element connectivity
        :obj:`list`: 1D list of element components
    """

    nodes = {'X' : [], 'Y' : [], 'Z' : []}
    elems = []
    comps = []
    
    min = 0.0 - length / 2.0
    max = 0.0 + length / 2.0
    
    for i in range (dims):
        for j in range (dims):
            for k in range (dims):
    
                if (i == 0) or (i == dims-1) or (j == 0) or (j == dims-1) or (k == 0) or (k == dims-1):
                    nodes ['X'].append (min + (i * length /(dims-1)))
                    nodes ['Y'].append (min + (j * length /(dims-1)))
                    nodes ['Z'].append (min + (k * length /(dims-1)))
    
    def rev_norm(inl):
        outl = []
        outl.append(inl[2])
        outl.append(inl[1])
        outl.append(inl[0])
        return outl

    # Bottom Face
    c = 0
    for i in range (dims-1):
        for i in range (dims-1):
            # Norms are good
            elems.append ([c + i, c + i+1, c + i+dims+1]       )
            elems.append ([c + i+dims+1, c + i+dims, c + i]    )
            comps.append(1)
            comps.append(1)
        c += dims
    
    # Front Face
    node_list = list(range (dims-1))
    c = (dims) * (dims)
    for i in range (dims-2):
        node_list.extend (list(range (c, c+dims-1)))
        c += (dims-1) * 4
    
    for i in node_list:
        # Norms are good once reversed
        if i < dims-1:
            elems.append (rev_norm([i, i+1, i+1+(dims*dims)]))
            elems.append (rev_norm([i+1+(dims*dims), i+(dims*dims), i]))
        if i > dims:
            elems.append (rev_norm([i, i+1, i+1+((dims-1)*4)]))
            elems.append (rev_norm([i+1+((dims-1)*4), i+((dims-1)*4), i]))
        comps.append(3)
        comps.append(3)
    
    # Top Face
    c = (dims) * (dims)
    c += ((dims-1) * 4) * (dims - 2)
    for i in range (dims-1):
        for i in range (dims-1):
            # Norms are good (when reversed)
            elems.append (rev_norm([c + i, c + i + 1, c + i + dims + 1]     ))
            elems.append (rev_norm([c + i + dims + 1, c + i + dims, c + i]  ))
            comps.append(2)
            comps.append(2)
        c += dims
        
    
    # Back Face
    c = (dims-1) * (dims)
    node_list = list(range (c,c+dims-1))
    for i in range (dims-2):
        c += (dims-1) * 4
        node_list.extend (list(range (c, c+dims-1)))
    
    c = 0
    for i in node_list:
        # These normals are good
        if c < ((dims - 1) * (dims - 2)):
            elems.append ([i, i+1, i+1+((dims-1)*4)])
            elems.append ([i+1+((dims-1)*4), i+((dims-1)*4), i])
        if c >= ((dims - 1) * (dims - 2)):
            elems.append ([i, i+1, i+1+(dims*dims)])
            elems.append ([i+1+(dims*dims), i+(dims*dims), i])
        comps.append(4)
        comps.append(4)
        c += 1
    
    # Left Face
    node_list = list(range (0,(dims*(dims-1)),dims))
    c = 0
    for i in range (dims-2):
        if i == 0:
            c += dims*dims
        if i > 0:
            c += (dims-1) * 4
        node_list.extend ([c])
        node_list.extend (list(range (c+dims, c+dims+(dims-2)*2, 2)))
    
    c = 0
    for i in range (dims-1):
        for j in range (dims-1):
            current = node_list [c]
            # Below are good
            if i == 0:
                if j == 0:
                    elems.append ([current, current+dims, (dims*dims)+dims])
                    elems.append ([dims+(dims*dims), (dims*dims), current])
                if j > 0:
                    elems.append ([current, current+dims, (dims*dims)+dims+(j*2)])
                    elems.append ([(dims*dims)+dims+(j*2), (dims*dims)+dims+((j-1)*2), current])
            if i > 0 and i < dims-2:
                if j == 0:
                    elems.append ([current, current+dims, current+dims+((dims-1)*4)])
                    elems.append ([current+dims+((dims-1)*4), current+((dims-1)*4), current])
                if j > 0:
                    elems.append ([current, current+2, current+2+((dims-1)*4)])
                    elems.append ([current+2+((dims-1)*4), current+((dims-1)*4), current])
            if i == dims-2:
                ref = (dims * dims) + (((dims - 1) * 4) * (dims - 2))
                if j == 0:
                    elems.append ([current, current+dims, current+((dims-1)*4)+dims])
                    elems.append ([current+((dims-1)*4)+dims, current+((dims-1)*4), current])
                if j > 0:
                    elems.append ([current, current+2, ref+dims+(j*dims)])
                    elems.append ([ref+dims+(j*dims), ref+dims+((j-1)*dims), current])
    
            comps.append(5)
            comps.append(5)
            c += 1
    
    
    # Right Face
    node_list = list(range (dims-1,(dims*(dims-1)),dims))
    c = 0
    for i in range (dims-2):
        if i == 0:
            c += dims*dims-1+dims
        if i > 0:
            c += (dims-1) * 4
        node_list.extend (list(range (c, c+((dims-1)*2), 2)))
    
    c = 0
    for i in range (dims-1):
        for j in range (dims-1):
            current = node_list [c]
            if i == 0:
            # Norms are good (when reversed)
                if j < (dims-2):
                    elems.append (rev_norm([current, current+dims, (dims*dims)+dims+1+(j*2)]))
                    elems.append (rev_norm([(dims*dims)+dims+1+(j*2), (dims*dims)+dims+1+((j-1)*2), current]))
                if j == (dims-2):
                    elems.append (rev_norm([current, current+dims, (dims*dims)-1+((dims-1)*4)]))
                    elems.append (rev_norm([(dims*dims)-1+((dims-1)*4), (dims*dims)+dims+1+((j-1)*2), current]))
            if i > 0 and i < dims-2:
                if j < (dims-2):
                    elems.append (rev_norm([current, current+2, current+2+((dims-1)*4)]))
                    elems.append (rev_norm([current+2+((dims-1)*4), current+((dims-1)*4), current]))
                if j == (dims-2):
                    elems.append (rev_norm([current, current+dims, current+dims+((dims-1)*4)]))
                    elems.append (rev_norm([current+dims+((dims-1)*4), current+((dims-1)*4), current]))
            if i == dims-2:
                ref = (dims * dims) + (((dims - 1) * 4) * (dims - 2)) -1
                if j < dims-2:
                    elems.append (rev_norm([current, current+2, ref+((j+2)*dims)]))
                    elems.append (rev_norm([ref+((j+2)*dims), ref+((j+1)*dims), current]))
                if j == dims-2:
                    elems.append (rev_norm([current, current+dims, ref+(dims*dims)]))
                    elems.append (rev_norm([ref+(dims*dims), ref+(dims*dims)-dims, current]))
            comps.append(6)
            comps.append(6)
    
            c += 1

    # Currently the ele object is 0 based
    return nodes, elems, comps

    
def build_sphere(radius, npts):
    """ Build the sphere grid by hand.

    The sphere is centered on 0,0,0.
    
    Args:
        radius (:obj:`int`): sphere radius
        npts   (:obj:`int`): Number of points along a hemisphere
        
    Returns:
        :obj:`dict`: Node X, Y, and Z array
        :obj:`list`: 2D list of element connectivity
        :obj:`list`: 1D list of element components
    """

    
    nodes = {'X' : [], 'Y' : [], 'Z' : []}
    elems = []
    comps = []
    
    #! Do the first arc from top to botton (npts)
    #! This is in the X/Y plane and each point is
    #! pi/(npts-1) radians apart
    for i in range(npts):
        arc = float(i)*np.pi/(float(npts)-1)
        nodes['X'].append(np.sin(arc) * radius)
        nodes['Y'].append(np.cos(arc) * radius)
        nodes['Z'].append(    0.0     * radius)
        
    #! OK, now we're going to march around the sphere using
    #! an identical spacing
    
    np_z = 2*npts - 1
    
    ### The i=0 case is addressed above, so start at i=1
    for i in range(1, np_z):
        arc = float(i)*2.0*np.pi/(float(np_z)-1)
        #print(i,np.degrees(arc))
        #! OK, we don't need the first and last points since they
        #! are axes.
        for j in range(npts-2):
            #! Get the elements sorted out
            my_i   = npts + (npts-2)*(i-1)  + j
            prev_i = my_i - (npts-2)

            ### The last arc duplicates the i=0 arc
            if i==np_z - 1:
                my_i = j + 1
            else:
                nodes['X'].append(nodes['X'][j+1] * np.cos(arc))
                nodes['Y'].append(nodes['Y'][j+1]              )
                nodes['Z'].append(nodes['X'][j+1] * np.sin(arc))

            if (i == 1):
                prev_i -= 1
    
            if (j == 0):
                #! Top Axis
                elems.append([my_i,prev_i,0])
            elif (j == npts-3):
                #! Bottom Axis
                elems.append([prev_i,my_i,npts-1])
                elems.append([my_i-1,my_i,prev_i])
                elems.append([prev_i,prev_i-1,my_i-1])
            else:
                elems.append([my_i-1,my_i,prev_i])
                elems.append([prev_i,prev_i-1,my_i-1])

    ### All components are set to '1' and the number of comps should equal
    ### the number of elements
    comps = [1 for elem in elems]

    return nodes, elems, comps
    
def write_simple_grid(nodes,elems,comps,type,base):
    """ Manually write out simple grid.

    Args:
        nodes (:obj:`dict`):   Node X, Y, and Z array
        elems (:obj:`list`):   2D list of element connectivity
        comps (:obj:`list`):   list of element components
        type  (:obj:`int`):    Write as tecplot ASCII (1) or ASCII .tri (2)
        base  (:obj:`str`): Basename for file ('base'.{tec,tri})
    """

    if type == 1:
        tec = open ("{}.dat".format(base), "w")
        tec.write ('VARIABLES = "X" "Y" "Z"\n')
        tec.write ("ZONE N = " + str (len(nodes['X'])) + \
                      ", E = " + str (len (elems))     + \
                   ", DATAPACKING = POINT, ZONETYPE = FETRIANGLE\n")
        for i in range (len (nodes ['X'])):
            tec.write("{} {} {}\n".format(nodes['X'][i],nodes ['Y'][i],nodes ['Z'][i]))
        for i in range (len (elems)):
            tec.write ("{} {} {}\n".format(elems[i][0]+1,elems[i][1]+1,elems[i][2]+1))
        tec.close ()

    if type == 2:
        tri = open ("{}.tri".format(base), "w")
        tri.write("{} {}\n".format(len(nodes['X']),len(elems)))
        for i in range (len (nodes ['X'])):
            tri.write("{} {} {}\n".format(nodes['X'][i],nodes ['Y'][i],nodes ['Z'][i]))
        for i in range (len (elems)):
            tri.write ("{} {} {}\n".format(elems[i][0]+1,elems[i][1]+1,elems[i][2]+1))
        for i in range (len (comps)):
            tri.write ("{}\n".format(comps[i]))
        tri.close ()

    return

###############################################################################

def half_to_full( half_mesh, mirror_plane, mirror_coord=0.0, tol=1e-8,
                  mirror_comp_name=None, return_mirror_plane=False ):
    """ Take a half-body Solution object and return a full-body version.

    Many CFD solutions and databases are constructed using only half the 
    vehicle in order to save disk space and computational time by leveraging
    symmetry.  When a full-body solution is needed for any reason, this 
    function can be used to deliver it.

    Args:
        half_mesh (:obj:`~steam.mesh.Mesh`): Mesh to mirror
        mirror_plane (:obj:`str`): Axis defining the mirror plane.  For example, to mirror a mesh about the x=10 plane, this input would simply be 'X'
        mirror_coord (:obj:`float`): Coordinate defining the mirror plane.  For example, to mirror a mesh about the x=10 plane, this input would simply be 10
        tol (:obj:`float`): Tolerance for how close a point must be to the mirror plane in order to be treated as on the plane and not be mirrored across it.
        mirror_comp_name (:obj:`str`): Optional name for the component containing all of the mirrored elements.
        return_mirror_plane (:obj:`bool`): Instead of returning the full-body mesh, return a PointCloud containing the points that fall within `tol` of the mirror plane.

    Returns:
        :obj:`~steam.mesh.Mesh`: The full-body mesh
    """

    mirror_plane = _check_plane_input( mirror_plane )

    ### Operating on the input mesh will alter it in place, so use a copy
    half = steam.mesh.Mesh( copy=half_mesh )
    
    ### Shift all pts so that the cut plane is at 0
    if mirror_coord:
        half.transform( [f't{-1.0*mirror_coord}{mirror_plane}'] )

    ### Get indices of pts on the cut plane
    on_plane_ind, off_plane_ind = _on_plane_indices( half, mirror_plane, 
                                                     mirror_coord, tol )
    if return_mirror_plane:
        plane = steam.mesh.Mesh( copy=half_mesh, is_point_cloud=True )
        plane.xyz_pt = half_mesh.xyz_pt.loc[on_plane_ind]
        return plane
    
    ### Create the points that will be mirrored across the plane
    mirror_pts = pd.DataFrame( half.xyz_pt.loc[ off_plane_ind ] )
    mirror_pts[ mirror_plane ] *= -1.0
    #   New indices are needed for these points
    mirror_pts[ 'orig_index' ] = mirror_pts.index
    max_orig_ind = max( half.xyz_pt.index )
    mirror_pts.index = range( max_orig_ind+1, 
                              max_orig_ind + mirror_pts.shape[0] + 1)
    
    ### Create the connectivity for the mirrored elements and reassign indices
    mirror_conn = half.conn.copy()
    max_orig_ind = max( half.conn.index )
    mirror_conn.index = range( max_orig_ind+1, 
                               max_orig_ind + mirror_conn.shape[0] + 1)
    #   Every value in mirror_conn that appears in off_plane_ind must be
    #   replaced by the new index of that mirrored point.
    alt_index = pd.Series( mirror_pts.index, index = mirror_pts.orig_index )
    for col in [0,1,2]:
        #   Column of entries that need to be swapped to a mirrored point
        swap_index = mirror_conn[ col ].index[ 
                                    mirror_conn[ col ].isin( off_plane_ind ) ]
        swap_col = mirror_conn.loc[swap_index][ col ]
        #   Perform the swap
        mirror_conn[ col ][swap_index] = alt_index[ swap_col ]
    
    ### All new elements should be added to a new component
    new_comp_num = int( max( half.comp.values ) ) + 1
    mirror_comp = np.ones( (mirror_conn.shape[0], 1), 
                           dtype=int ) * new_comp_num
    
    ### Points and connectivity need to be recombined and then the "orig_index"
    #   column needs to be removed
    full_pts = pd.concat( (half.xyz_pt, mirror_pts), sort=False )
    full_pts.drop( columns = ('orig_index'), inplace=True )
    full_conn = pd.concat( (half.conn, mirror_conn) )
    full_comp = pd.DataFrame( np.vstack( (half.comp.values, mirror_comp) ),
                              index=full_conn.index )
    
    ### Assemble all the pieces into a mesh
    full_mesh = steam.mesh.Mesh()
    full_mesh.xyz_pt = full_pts
    full_mesh.conn = full_conn
    full_mesh.comp = full_comp

    if mirror_comp_name:
        full_mesh.set_comp( mirror_comp_name, new_comp_num )

    ### Undo the coordinate shift that was done on the half mesh
    if mirror_coord:
        full_mesh.transform( [f't{mirror_coord}{mirror_plane}'] )

    ### Store information inside the mesh for use in solution.half_to_full
    #   We need a "decoder ring" for each new element to tell us the element
    #   ID from which it is mirrored.
#    full_mesh.elem_map = dict( zip(half.conn.index, mirror_conn.index) )
    full_mesh._elem_map = pd.Series( mirror_conn.index, index=half.conn.index)
    full_mesh._node_map = pd.Series( mirror_pts.index, index=off_plane_ind )
#    interact( local = dict( globals(), **locals() ) )
#    full_mesh.alt_elem_map = pd.Series( half.conn.index, 
#                                        index=mirror_conn.index)

    return full_mesh

###############################################################################

def _check_plane_input( plane ):
    """ Confirm that a given value for a plane coordinate is valid.

    This function simply checks that it's a `:obj:str` that matches "X", "Y",
    or "Z" and returns the upper case version of that letter.
    """ 

    ### `plane` has to be a string
    try:
        out_plane = plane.upper()
    except AttributeError:
        raise TypeError( f'"plane" input must be str.\nValueGiven:'
                       + f'              {plane}\n' )

    ### Don't allow bad values of plane
    if out_plane not in ('X', 'Y', 'Z'):
        raise ValueError( "'plane' input must be either 'X', 'Y', "
                          + f"or 'Z'.\nValue given:     {plane}" )

    return out_plane

###############################################################################

def _on_plane_indices( mesh, plane, plane_coord, tol=1e-8 ):
    """ Find which points on a Mesh lie on a plane.

    This function, which is primarily intended for use with 
    `~steam.mesh.half_to_full` and `~steam.mesh.half_to_full`, takes a Mesh,
    a plane axis and coordinate, and a tolerance and returns the indices that
    are on and off of that plane.


    Args:
        mesh (:obj:`~steam.mesh.Mesh`): The mesh of interest
        plane (:obj:`str`): Axis defining the plane.
        plane_coord (:obj:`float`): Coordinate defining the plane.  For example, to find points on the x=10 plane, this input would simply be 10
        tol (:obj:`float`): Tolerance for how close a point must be to the plane in order to be treated as "on" it.

    Example:  To find points on the x=10 plane, one would call:
    "_on_plane_indices( mesh, 'X', 10 )"
    """

    plane = _check_plane_input( plane )

    ### Get indices of pts on the cut plane
    on_plane_pts = abs( mesh.xyz_pt[ plane ] - plane_coord ) <= tol 
    on_plane_ind = mesh.xyz_pt.index[ on_plane_pts ]
    off_plane_ind = mesh.xyz_pt.index[ ~on_plane_pts ]

    return on_plane_ind, off_plane_ind
    
###############################################################################

#class PointCloud( Mesh ):
#    """ Class definition for STEAM PointCloud
#
#    A PointCloud is a derived class from Mesh, so that it maintains much of
#    the same functionality as a mesh.  In addition to the standard mesh 
#    attributes, a PointCloud has a solution attached to it that stores fixed
#    data that won't be interpolated.  These can be body point names/numbers,
#    for example.
#
#    The mesh is stored in a DataFrame with minimal values of
#    X, Y, and Z.
#    """
#
#    def __init__(self,file=None,filetype=None,copy=None):
#        """Constructor for Mesh class.
#
#            copy=MESH will copy a mesh
#        """
#
#        ### Call the constructor for Mesh class.
#        Mesh.__init__(self, file = file, filetype = filetype, copy = copy )
#
#        ### Convert the mesh to a PointCloud
#        self.to_point_cloud()
#
#        ### Define a static solution to store the data like point names
#        import steam.solution
#        self.static_soln = steam.solution.Solution(mesh = self)
#
#    ###############################################
#
#    def to_point_cloud(self):
#        """ This will convert a standard mesh to a PointCloud-style mesh.
#
#        Element connectivity is thrown out and replaced with
#        a one-to-one with the points.
#
#        Components are all set to 1 (TBR)
#
#        Certain operations are no longer possible.
#        """
#        
#        self.xyz_el = self.xyz_pt
# 
#        #! Simple one-to-one connectivity
#        npt         = self.xyz_el.values.shape[0]
#        conn        = np.arange(npt)
#        self.conn   = pd.DataFrame(conn)
#
#        comp            = np.ones(npt,dtype='int')
#        self.comp       = pd.DataFrame(comp)
#        self.comp_table = dict()
#
#        self.is_point_cloud = True
#
#        self.init()
#
#    ###############################################
#
#    def read(self,file,filetype=None):
#        """Method to read a new mesh object.
# 
#        Args:
#            file (:obj:`str`): Path to file.
#            filetype (:obj:`str`): Type of file. See table below:
#
#        Here are  the file types::
#
#                DF   - DataFrame, space delinated file with column headers
#                CDAT - DataFrame, space delinated file with column headers
#                DAT  - Tecplot ASCII File
#                TRI  - ASCII .tri  File
#                TRIQ - ASCII .triq File
#                PKL  - Pickled Mesh File 
#                *    - Pass to libmesh
#        """
#        
#        ftu = filetype.upper()
#        if   (ftu == "CDAT" or ftu == "DF"):
#            self.read_cdat(file)
#            self.static_soln.read_cdat(file)
#        elif (ftu == "DAT"):
#            self.read_tecplot_ascii(file)
#        elif (ftu == "TRI" or ftu == "TRIQ"):
#            steam.io.read_triq_ascii(file,ioSoln=True,mesh=self, 
#                                     soln = self.static_soln)
#        elif (ftu == "UTRI" or ftu == "UTRIQ"):
#            steam.io.read_triq_uform(file,ioSoln=True,mesh=self)
#        elif (ftu == "PKL"):
#            steam.io.read_pickle(file,mesh=self)
#        else:
#            self.read_libmesh(file)
#
#        if (self.comp is None):
#            print("* Error: mesh filetype not recognized")
#            raise IOError
#
#        # Populate mesh information
#        self.init()
#        
#    ###############################################
#
#    def slice(self, *args, **kwargs):
#        """ Method to replace Mesh.slice because PointCloud objects 
#        cannot be sliced.
#        
#        This method simply raises an AttributeError alerting the user that
#        PointCloud objects cannot be sliced.
#        """
#
#        raise AttributeError("PointCloud objects cannot be sliced!")
#
