""" Module for the Database class
"""
import logging
logger = logging.getLogger(__name__)
logger.debug("Top of file")

import pandas as pd
import numpy  as np
import steam
import os
from copy import deepcopy

### Import "code" module for debugging
#from code import interact
#interact( local = dict( globals(), **locals() ) )

# This turns on some extra debugging in this module if you need it
module_debug = False

class Database(steam.table.Table):
    """Class definition for STEAM Database

    The Database class is an inherited class from Table.

    A Database is a collection of data that is indexed by a number
    of independent variables.  The data can be dependent variables
    or mesh/solution objects.  Each database can have a Delaunay 
    interpolator object.

    Each database has a self.data DataFrame object that contains
    the independent and dependent variables as well as pointers
    to any surface maps that exist.  Similarly, each database can
    contain any number of meshes.  They are stored in a dictionary 
    based on their unique meshID, self.mesh.  
    The solution data is stored on disk or, optionally, in
    memory.
    """

    def __init__(self,data=None,description=''):
        """Constructor for Database class.

        Args:
          data (:obj:`str`): CDAT File to read or DataFrame to inheret
          description(:obj:`str`): string description of the database
        """

        self.description = description    # Metadata or description

        ### Run the Table class __init__ method
        try:
            steam.table.Table.__init__( self, data=data)
        except TypeError:
            ### data may be a list of DBPoint objects
            ### Assume this is the case and allow errors to be raised if not
            steam.table.Table.__init__( self, data=None)

        self.meshes     = {}             # Mesh Object(s)
        self.solns      = {}             # Dictionary of map DataFrames
        self.intspace   = {}             # Interpolation Spaces

        ### If this is going to be a large database, then I might want
        ### to immediately store it on disk instead of putting it into memory.
        self.onDisk         = False
        self.__check_opened = False # Did onDisk_check_open open this store?

        #! If I passed in data, store the base directory since all paths should
        #! be relative to that.  I only need this when loading data, I don't save it.
        import os.path
        self.fileBasePath = "."
        if (data is None) or isinstance(data, pd.DataFrame):
            pass
        elif isinstance(data, list):
            for db_pt in data:
                self.add_point( db_pt )
        elif os.path.isfile( data ):    
            path = data.replace("//","/")
            self.fileBasePath = "{}".format(os.path.dirname(path))
            if (self.fileBasePath == ''):
                self.fileBasePath = "."

    def __str__(self):
        """Return the print string of the Database object."""

        string  = "Database Object:"
        if self.description != "":
            string += "\n "+self.description

        string += "\n Independent Variables: "+str(self.indep)

        string += "\n Dependent Variables  : "+str(self.dep)

        string += "\n DataFrame            : "
        if (self.data is not None):
            string +=str(self.data.shape)
        else:
            string +="None"

        string += "\n Meshes               : {} meshes".format(len(self.meshes))
        for (key,val) in self.meshes.items():
            string += "\n    {} : {:7d} points".format(
                key,
                self.meshes[key].xyz_pt.shape[0]
                )

        string += "\n Solution Maps        : {} solutions".format(len(self.solns))
        # Only print the first and last XXX solutions as an example
        nprint = 1
        i = 0
        for (key,val) in self.solns.items():
            i += 1;

            if (i == nprint+1):
                string += "\n    ...";
            if i > nprint and i < (len(self.solns)):
               continue

            if (self.solns[key].data is not None):
                string += "\n    {} : {:d} points, {:3d} vars".format(
                    key,
                    self.solns[key].data.shape[0],
                    self.solns[key].data.shape[1]
                    )
            else:
                string += "\n    {} : Data not loaded".format(
                    key)

        string += "\n HDF5 File            : "
        if (self.hdf5_file is not None):
            string +=self.hdf5_file
        else:
            string +="None"

        string += "\n Work on Disk         : "+str(self.onDisk)

        string += "\n Interpolation Spaces : {} spaces".format(len(self.intspace))
        for (key,val) in self.intspace.items():
            string += "\n    {} : {:7d} points".format(
                key,
                len(self.intspace[key].out_index)
                )

        return string
    
    def __repr__(self):
        return "STEAM Database Object"

    ### Iterators
    def __iter__(self):
        """ Iterator initialization.

        Keeps track of the objects by using an iter_count variable.
        It also limits the iterations to the points that are in the
        database at the time of the call
        """

        self._index_iterator = self.data.index.__iter__()
        return self

    def __next__(self):
        """ Return the next object as a DBPoint."""

        ### Should automatically raise StopIteration at the right time
        next_index = self._index_iterator.__next__()

        point = self.get_point( next_index )

        return point

    def __del__(self):
        """ Some extra delete calls."""

        if self.onDisk:
            self.onDisk_check_close()

    def __eq__( self, other ):
        """ Used for equality determination.
        """
    
        if isinstance( other, steam.database.Database ):
            raise NotImplementedError( 'Database equality not comparable ' +
                  'with "==" operator at this time.  See Gitlab issue #27.' )

        return False

    ### onDisk coordination
    @steam.util.timer(logger)
    def onDisk_check_open(self):
        """ Simple check that opens the store if it was closed.

        If the store is already open, it does nothing.  This routine
        and the close routine are called every time you try and work with
        onDisk data.  So they can see a lot of use!
        """

        if self.hdf5 is None:
            raise IOError("I need to have an HDF5 file initialized "+
                          "when reading straight to disk!")

        #! If the HDF5 file is closed, then I need to reopen it
        if (self.hdf5._handle is None):
            self.hdf5.open(mode=self.hdf5_mode)
            if module_debug:
                logging.debug("Opening HDF5, mode {}".format(self.hdf5_mode))
            self.__check_opened = True

    def onDisk_check_close(self):
        """ Simple check that closes the store if it was opened by .onDisk_check_open.

        If the container is open, then keep things open.  However, if the database
        module opened the store, then it should also close it.
        """
        if self.__check_opened:
            if self.hdf5._handle is not None:
                self.hdf5.close()

    def onDisk_to_mem(self):
        """ Convert from an onDisk database to one in memory.

        This will read in all of the solution data from disk and
        complete the database in memory.  It can then be manipulated
        in memory or written to a new container.

        ToDo: Create an onDisk_copy_to() method or similar to facilitate
        copying without having to read into memory.
        """

        for sid in self.solns:
            self.solns[sid].data = self.get_soln_df(sid)

        self.onDisk_check_close()
        self.hdf5      = None
        self.hdf5_mode = None
        self.hdf5_root = None
        self.onDisk    = False
        

    ### ## Pandas Wrappers / Coordination

    def add_var(self, var, data, index_in=None):
        """ Add a new variable based on a list of data.

        Args:
            var  (:obj:`str`): Name of variable.
            data (:obj:`list`): List of data, needs to match existing size (point/element).
            index_in (:obj:`pandas.core.indexes.range.RangeIndex`): Optional index to apply to the `data` input.  This input is necessary when appending a variable to a database with a discontinuous index or when a variable will not be defined at all points in the database.
        """

        new_df = pd.DataFrame( list(data), columns=[var], index=index_in )

        if (self.data is None):
            self.data = new_df
            return

        if (var in self.data.columns):
            raise Exception("Variable {} already exists.  Remove it before adding.".format(var))

        self.data = self.data.join(new_df)

    def get_soln_df(self,solnid,cols=None,index=None):
        """ Returns a DataFrame of a specified solnID.
        
        This provides a combined syntax that can be used for both
        in-memory and on-disk datasets.

        **Note:**  This returns the dataframe, only!  Not the entire solution.

        **Note:**  This is a copy of the actual objects, the only way that
        it would be possible to return the actual object would be
        to have the database in-memory and not use either of the
        keyword arguments (cols,index).  We desire consistency.
        
        Args:
         solnid (:obj:`str`): solutionID to get.  e.g. 'a000528'
         cols  (:obj:`list`, optional): List of columns from the soln.
         index (:obj:`list`, optional): List of indexes to get from the soln.

        Returns:
         (:obj:`~pandas.DataFrame`): Specified solution map with specified or all variables.
        """
        
        #! If the HDF5 file is closed, then I need to reopen it
        if self.onDisk:
            self.onDisk_check_open()

        ### Return it all if vars is set to none
        if cols is None:
            if index is None:
                if module_debug:
                    logger.debug("get_soln_df {} None None".format(solnid))
                # Return everything
                if self.onDisk:
                    output_df = self.hdf5.select(
                                    self.hdf5_root+'/solns/soln_'+solnid
                                           )
                else:
                    output_df =  self.solns[solnid].data
            else:
                if module_debug:
                    logger.debug("get_soln_df {} None {}".format(solnid,len(index)))
                # Return a subset at index locations
                if self.onDisk:
                    output_df =  self.hdf5.select(
                                    self.hdf5_root+'/solns/soln_'+solnid,
                                    where="index == index"
                                           )
                else:
                    s = self.solns[solnid].data
                    output_df =  s[s.index.isin(index)]

        ### If vars was set, then return a subset *as a dataframe*
        else:
            if index is None:
                if module_debug:
                    logger.debug("get_soln_df {} {} None".format(solnid,cols))
                # Return everything, but only certain variables
                if self.onDisk:
                    output_df =  self.hdf5.select(
                                    self.hdf5_root+'/solns/soln_'+solnid,
                                    columns=cols,
                                           )
                else:
                    output_df =  pd.DataFrame(self.solns[solnid].data[cols])
            else:
                if module_debug:
                    logger.debug("get_soln_df {} {} {}".format(solnid,cols,len(index)))
                # Return a subset at index locations, but only certain variables
                if self.onDisk:
                    output_df =  self.hdf5.select(
                                    self.hdf5_root+'/solns/soln_'+solnid,
                                    where="index == index",
                                    columns=cols
                                           )
                else:
                    s = self.solns[solnid].data
                    output_df =  pd.DataFrame(s[s.index.isin(index)][cols])

        # Return a copy
        return pd.DataFrame(output_df.values.copy(), 
                            output_df.index.copy(), 
                            output_df.columns.copy())

    def get_soln(self,solnid,get_data=True):
        """ Returns a Solution object into memory.

        This is useful since we need to disambiguate between solution data that
        is stored on disk and data that is stored into memory.  The result from
        this function will be a complete solution stored in memory.

        get_data really only has meaning for onDisk databases.  If the solutions
        are in memory, then you will always get data.

        Args:
         solnid (:obj:`str`): solutionID to get.  e.g. 'a000528'
         get_data (:obj:`bool`): Do you want populated data? [T]

        Returns:
         (:obj:`~steam.solution.Solution`): Specified solution.
        """

        if (solnid not in self.solns):
            raise ValueError(
                    "Solution ID '{}' not found in database".format(solnid) )
        #! This is a two-step process because we may store the DF on disk 
        #! or in memory and get_soln_df tracks that.
        soln = steam.solution.Solution(copy=self.solns[solnid])
        if get_data:
            soln.data = self.get_soln_df(solnid)  # This is a copy
        
        return soln

    def get_mesh(self,meshid):
        """ Returns a Mesh object from database to memory.

        Args:
         meshid (:obj:`str`): solutionID to get.  e.g. 'a000528'

        Returns:
         (:obj:`~steam.solution.Mesh`): Specified mesh.
        """

        if (meshid not in self.meshes):
            raise Exception("Mesh ID '{}' not found in database".format(
                                                                    meshid))
        #! This is a two-step process because we may store the DF on disk 
        #! or in memory and get_mesh_df tracks that.
        mesh = steam.mesh.Mesh(copy=self.meshes[meshid])
        
        return mesh

    ### ## Data Input Methods

    def read_data(self,data=None,mesh=None,ioSoln=True):
        """ Read database input data file(s).

        data file is a columnized file with a single header line.
        The headers should include the independent and dependent
        variables and associated solutions.  If there are solution
        maps, then the 'Soln_Path' and 'Soln_Type' variables are
        required.  'Mesh_Name' must also be specified if the solution
        does not have mesh stored in the same file.

        mesh is a columnized file that must have the following headers 
        and an entry for each mesh to be read.  The Mesh_Name must
        match the Mesh_Name specified for the solutions and it must
        begin with a letter 'g001', for example. ::

                Mesh_Name Mesh_Type Mesh_Path
                   ...       ...       ...
                   ...       ...       ...

        Mesh types are outlined in :meth:`Mesh.read <steam.mesh.Mesh.read>`. 

        Either file can have comments prior to the header line.
        Any line that begins with a '#' is skipped during the load.

        Following the data read, a consistency check is made to 
        ensure that all grids are loaded and that necessary
        variables are defined.

        Args:
            data (:obj:`str`): CDAT file name
            mesh (:obj:`str`): CDAT file name
            ioSoln (:obj:`bool`): Load the solutions automatically? [T]

        """

        #! Read the data
        if (data is not None):
            self.read_cdat(data)

        #! Read the meshes
        if (mesh is not None):
            try:
                meshlist = pd.read_table(mesh,delim_whitespace=True,comment='#')
            except:
                raise IOError("Cannot open mesh list file!\n"+" File: "+mesh)
       
            ### Who knows what capitalization was used, suss it out first
            try:
                pathvar = steam.util.get_var_cap(meshlist.columns,"mesh_path")
                typevar = steam.util.get_var_cap(meshlist.columns,"mesh_type")
                namevar = steam.util.get_var_cap(meshlist.columns,"mesh_name")
            except:
                raise IOError("Cannot find all required mesh headers in mesh file!")

            # Now load all of the meshes
            for (i,item) in meshlist.iterrows():
                filepath = self.fileBasePath+"/"+item[pathvar]
                if os.path.isabs(item[pathvar]):
                    filepath = item[pathvar]
                meshid_out = self.read_mesh(
                                        filepath,
                                        item[typevar],
                                        meshid=item[namevar]
                                        )
                # If a meshid changes, then let the user know since that means
                # that their mapping might be incorrect.
                if (meshid_out != item[namevar]):
                    print("MeshID '{}' was changed to '{}' due to namespace conflicts!"
                          .format(item[idvar],meshid_out))


        #! Check for solution variables.  If they exist, then
        #! check to make sure that we have all of the grid inputs
        #! that we need.
        if ioSoln:
            self.load_soln()

    def read_cbaero_adb( self, adb_file, alt_mesh = None, 
                         trans_method = 'INV_DIST', just_mesh=False ):
        """ Read a CBAero \*.adb binary database file directly into a database.

        This method is based on the adb2tecplot function provided with CBAero
        version 3.9.0.  The source code for that function is here:
        /software/x86_64/cbaero/cbaero.3.9.0.distro.1.28.2015/tools/adb2tecplot/adb2tecplot.C

        The resulting database has the following columns: ::

                mach    qbar(bar)  alpha(deg) press_inf   fs_qbar    sym_flag   stag_tri 
                ...        ...        ...        ...        ...        ...        ...    

        Args:
            adb_file (:obj:`str`): The CBAero file to read
            alt_mesh (:obj:`~steam.mesh.Mesh`): Alternate mesh on which to build the database.  Useful for building database on body points.
            trans_method (:obj:`str`): String flag indicating how the transformation to alt_mesh should be done.  Possible values: 'INV_DIST' (Default) for inverse distance and 'QUAD' for inverse distance quad
            just_mesh (:obj:`bool`): Read only the mesh from the CBAERO file and return?
        """
        
        logger.info('Read database from CBAero file:   {:s}'.format(adb_file))

        ### Define dtypes to store all data associated with a tri or a
        ### node so that we can read tris and nodes all at once.
        ### dtypes do not need to be defined within the context manager
        tri_type = [('n1',np.int32), ('n2',np.int32), ('n3',np.int32), 
             ('surf_id',np.int32), ('area',np.float32), 
             ('mat_id',np.int32), ('emiss',np.float32)]
        node_type = [('x',np.float32), ('y',np.float32), ('z',np.float32)]
        soln_header_type = [('stag_tri', np.int32),
                            ('fs_press', np.float32), 
                            ('fs_qbar',  np.float32) ]
        ### For the soln_type, an attempt was made to make variable names 
        ### consistent with BLayer
        soln_type = [('Tw (K)',        np.float32), # Wall Temperature, K
                     ('pw (Pa)',       np.float32), # Wall Pressure, Pa
                     ('rhow (kg/m^3)', np.float32), # Wall Density, Pa
                     ('Te (K)',        np.float32), # Edge Temperature, K
                     ('rhoe (kg/m^3)', np.float32), # Edge Density, Pa
                     ('Me',            np.float32), # Edge Mach number
                     ('Ae (m/s)',      np.float32), # Edge Speed of sound, m/s
                     ('CH (kg/m^2.s)', np.float32), # Film Coefficient, kg/(m^2-s)
                     ('Hr (J/kg)',     np.float32), # Recovery Enthalpy, J/kg
                     ('qw (W/m^2)',    np.float32), # Convective Qdot, Watts/m^2
                     ('qsr (W/m^2)',   np.float32), # Shock Radiation Qdot, Watts/m^2
                     ('taux (Pa)',     np.float32), # Shear stress X, Pa
                     ('tauy (Pa)',     np.float32), # Shear stress Y, Pa
                     ('tauz (Pa)',     np.float32), # Shear stress Z, Pa
                     ('s (m)',         np.float32), # Running Length, m
                     ('delta (m)',     np.float32), # Boundary layer thickness, m (also transition code based on sign)
                     ('Re-theta',      np.float32)  # Momentum Thickness Rynolds number
                    ]
        ### open the file and read it
        with open( adb_file, 'rb' ) as f:
            ### Read the endian-ness check digit:
            dum_int, num_nodes, num_tris, num_mach, num_qbar, num_alpha = (
                    np.fromfile( f, dtype = np.int32, count = 6 ) )

            ### Check that the dummy integer is within 10 of Dave Kinney's
            ### pre-specified magic number.  It has to be within 10 because
            ### apparently there are different values depending on the version.
            if not abs(dum_int - (-123789456)) < 10:
                raise ValueError( 'The *.adb file is likely written in the ' +
                                  'wrong endianness for this machine.' )
            
            ### Determine adb_version based on dum_int
            ### This logic is copied directly from adb2tecplot
            adb_version = 1 if ( dum_int == -123789456 ) else 2

            logger.info("Number of Nodes:     {:d}".format(num_nodes) )
            logger.info("Number of Tris:      {:d}".format(num_tris) )
            logger.info("Number of Machs:     {:d}".format(num_mach) )
            logger.info("Number of QBars:     {:d}".format(num_qbar) )
            logger.info("Number of Alphas:    {:d}".format(num_alpha) )

            ### Skip 6 unnecessary floats
            f.seek( np.float32().itemsize * 6, 1 )
            
            ### Read Machs, QBars, and Alphas
            #   qbar is in bars
            #   alpha is in degrees
            mach_list = np.fromfile( f, dtype = np.float32, count = num_mach )
            qbar_list = np.fromfile( f, dtype = np.float32, count = num_qbar )
            alpha_list = np.fromfile(f, dtype = np.float32, count = num_alpha)

            ### Read the mesh
            #   Read all tris
            tris = np.fromfile( f, dtype = tri_type, count = num_tris )
            const_tri_data = { var: tris[var] for var in 
                               ['surf_id', 'area', 'mat_id', 'emiss'] }
            #   Read all nodes
            nodes = np.fromfile( f, dtype = node_type, count = num_nodes )

            ### Make a Mesh object with node and tri data
            grid = steam.mesh.Mesh()
            node_dict = {}
            node_dict['X'] = nodes['x']
            node_dict['Y'] = nodes['y']
            node_dict['Z'] = nodes['z']
            tri_dict = {i: [tri['n1']-1, tri['n2']-1, tri['n3']-1] for i, tri
                            in enumerate(tris) }    # -1 to get to 0-index
            surf_dict = {i: tri['surf_id'] for i, tri in enumerate(tris) }
            grid.xyz_pt = pd.DataFrame.from_dict( node_dict )
            grid.conn   = pd.DataFrame.from_dict(tri_dict, 'index')
            grid.comp   = pd.DataFrame.from_dict(surf_dict, 'index')
            grid.get_xyz_el()

            ### Return the mesh if that's all that's desired
            if just_mesh:
                return grid

            ### If desired, construct a transformation to the alternate mesh
            if alt_mesh is not None:
                logger.debug( 'Constructing database on alternate mesh ' + 
                              'instead of mesh contained in *.adb file' )
                transform = steam.interpolate.Transform( source = grid, 
                                                         target = alt_mesh )
                if trans_method.upper() == 'INV_DIST':
                    transform.inverse_distance()
                elif trans_method.upper() == 'QUAD':
                    transform.inverse_distance_quad()
                else:
                    raise ValueError( 'trans_method must be either ' +
                                      '"INV_DIST" or "QUAD"' )

            ### Skip a dummy integer and a dummy float
            f.seek( np.int32().itemsize + np.float32().itemsize, 1 )
            ### Symmetry flag
            sym_flag = np.fromfile( f, dtype = np.int32, count = 1 )[0]

            num_solns = len(alpha_list) * len(mach_list) * len(qbar_list)
            soln_num = 1
            ### Loop over all freestream independent variables 
            for alpha in alpha_list:
              for mach in mach_list:
                for qbar in qbar_list:
                    steam.util.progress_bar( soln_num, num_solns,
                                             prefix="Read CBAERO Solns :" )
                    soln_num += 1

                    ### Read solution data
                    stag_tri = np.fromfile( f, dtype = np.int32, count=1)[0]
                    fs_press = np.fromfile( f, dtype = np.float32, count=1)[0]
                    fs_qbar  = np.fromfile( f, dtype = np.float32, count=1)[0]
                    tri_data = np.fromfile( f, dtype = soln_type,
                                               count = num_tris )

                    ### Build a solution object
                    soln = steam.solution.Solution(mesh = grid)
                    soln.store_at = 'ELEMENT'
                    data_dict = { var[0]: tri_data[var[0]] 
                                  for var in soln_type }

                    ### Calculate boundary layer thickness and transition
                    ### again based on the logic in adb2tecplot
                    if adb_version == 1:
                        from math import floor
                        turbulence = floor( data_dict['delta (m)'] ) /1000.0
                        turbulence[ turbulence < 0.1 ] = 0.0
                        turbulence[ turbulence >= 0.1 ] = 1.0
                        data_dict['turbulence'] =  turbulence
                        data_dict['delta (m)'] = np.zeros_like( 
                                                    data_dict['delta (m)'])
                    else:
                        turbulence = np.zeros_like( data_dict['delta (m)'] )
                        turbulence[data_dict['delta (m)'] < 0.] = 1.0
                        data_dict['turbulence'] =  turbulence
                        data_dict['delta (m)'] = abs( data_dict['delta (m)'] )

                    ### Add constant parameters from the mesh
                    data_dict.update( const_tri_data )

                    soln.data = pd.DataFrame.from_dict( data_dict )
                    soln.mesh = grid
                    soln.init()

                    ### Move the solution to the alternate mesh if necessary
                    out_mesh = grid
                    out_soln = soln
                    if alt_mesh is not None:
                        out_mesh = alt_mesh
                        out_soln = transform.apply( soln )

                    ### Create a DBPoint and add to the current database
                    data = pd.Series( {'mach': mach, 'qbar(bar)': qbar, 
                                   'alpha(deg)': alpha, 'press_inf': fs_press,
                                   'fs_qbar': fs_qbar, 'sym_flag': sym_flag,
                                   'stag_tri': stag_tri } )
                    point = steam.database.DBPoint( data=data, mesh=out_mesh, 
                                                    soln = out_soln)
                    self.add_point( point )

    def load_soln(self,remove_cols=False):
        """ Look through the data DataFrame and find any solns that
        are specified.  Read them into memory or straight to the
        HDF5 file if 'onDisk' is non-empty in database.

        This depends on the following columns being populated in the
        data table::

            soln_path soln_type

        If there is more than one non-contained mesh (not included in the
        solution format) used, then each entry in the
        data table also needs ``Mesh_Name`` defined.

        Capitalization for soln_path, soln_type, and mesh_name do not
        matter and they are renamed to all lower_case upon load.  This
        is to make things more consistent for subsequent operations.

        Args:
            remove_cols (:obj:`bool`): Remove the soln_path, soln_type, and mesh_name columns? [F]

        """
        
        ### Add a meshid variable to DataFrame if it isn't there.
        ### If it isn't there, then we assume that there is only
        ### one mesh and that all maps use it.  Check this assumption.
        try:
            meshvar = steam.util.get_var_cap(self.data.columns,"mesh_name")
        except:
            meshvar = "meshid"

        self.data = self.data.assign(meshid = lambda x: None)

        ### Disambiguate capitalization for variables I need
        try:
            pathvar = steam.util.get_var_cap(self.data.columns,"soln_path")
            typevar = steam.util.get_var_cap(self.data.columns,"soln_type")
        except:
            raise IOError("Cannot find all required soln headers!")

        ### Rename the variable to make things more consistent in memory
        self.data.rename(
                         columns={
                                  pathvar:"soln_path",
                                  typevar:"soln_type",
                                  },
                         inplace=True
                         )

        pathvar = "soln_path"
        typevar = "soln_type"

        ### Add a solnid variable to DataFrame
        self.data = self.data.assign(solnid = lambda x: None)

        ### Iterate through all of the rows of data and load each
        ### solution into memory or onto disk.
        for (i,item) in self.data.iterrows():
            # print(i,item[typevar],item[pathvar])
            steam.util.progress_bar(i+1,self.data.shape[0],prefix="Loading Solns :")

            # Load the solution
            filetype = item[typevar]
            filepath = self.fileBasePath+"/"+item[pathvar]
            if os.path.isabs(item[pathvar]):
                filepath = item[pathvar]

            try:
                soln = steam.solution.Solution(file=filepath,filetype=filetype)
            except:
                raise IOError("Cannot open solution file!\n"+" File: "+file)

            # Check to see if we loaded a mesh, too:
            meshid = item[meshvar]
            if (soln.mesh is not None):
                meshid = self.add_mesh(soln.mesh)
            else:
                if (meshid is None):
                    if len(self.meshes) != 1:
                        raise IOError("Implicit mesh to solns match requires that" +
                               " there is one and only one mesh!\n   You probably" +
                               " need to specify a MeshID column in database input.")
                        return
                    else:
                        meshid = list(self.meshes)[0]
                mesh     = self.meshes[meshid]
                soln.mesh = mesh

            ### Append the solution to disk or memory
            solnid = self.add_soln(soln)
            self.data.at[i,"solnid"] = solnid
            self.data.at[i,"meshid"] = meshid

        return
 

    def read_mesh(self, file, filetype, meshid=None):
        """ Read a mesh of the appropriate file type.  This will append
        it to the dictionary of meshes for use with solution maps.

        Args:
            file (:obj:`str`):     Path to mesh to read
            filetype (:obj:`str`): See :meth:`Mesh.read <steam.mesh.Mesh.read>` for file types.
            meshid (:obj:`str`,optional): MeshID to use when storing
                the mesh.  Otherwise, unique default will be found.
        Returns:
            :obj:`str`: The MeshID that was used to store the mesh.

        """

        ### Read in the mesh.
        try:
            mesh = steam.mesh.Mesh(file,filetype)
        except:
            raise IOError("Cannot open mesh file!\n"+" File: "+file)
        
        mesh.info = file

        #! add_mesh handles MeshID conflicts, appending to the database's
        #! list of meshes and returns the unique MeshID
        meshid_out = self.add_mesh(mesh,meshid)

        # Conclude
        print("Read mesh: "+file)

        return meshid_out

    def compare_hash(self,obj1,obj2,compare_obj=False):
        """ Simple routine that will compare two obj hashes for equality.

        This is its own routine since we have the capability to check hashes and look for
        a hash conflict.


        Args:
            mesh1 (:obj:`~steam.mesh.Mesh` or :obj:`~steam.solution.Solution`): First Object
            mesh2 (:obj:`~steam.mesh.Mesh` or :obj:`~steam.solution.Solution`): Second Object
            compare_obj (:obj:`bool`): Perform hash conflict compare on matches? [F]

        Returns:
            :obj:`bool`: True / False for equlity
        """
    
        if type(obj1) == steam.mesh.Mesh:
            pass
        elif type(obj1) == steam.solution.Solution:
            pass
        else:
            raise TypeError ("Only a Mesh or Solution can be compared!")
    
        if type(obj2) == steam.mesh.Mesh:
            pass
        elif type(obj2) == steam.solution.Solution:
            pass
        else:
            raise TypeError ("Only a Mesh or Solution can be compared!")
    
        equal = False
        if obj1.hash == obj2.hash:
            ##! Check to make sure we don't have a hash conflict
            ##! For the record, I have never seen this happen AMS (6/20/18)
            if compare_obj:
                if obj1 != obj2:
                    raise AssertionError ("Hash conflict error when comparing objects!")
            equal = True
        return equal
 

    def find_obj_hash(self,obj,compare_obj=False):
        """ Simple routine that will check to see if a mesh/soln's hash is already in the database.

        We have the capability to check hashes and look for a hash conflict.

        The return will be None if it is new, or the object id if it is a duplicate.

        Args:
            obj (:obj:`~steam.mesh.Mesh` or :obj:`~steam.solution.Solution`): Object to look for
            compare_obj (:obj:`bool`): Perform hash conflict compare on matches? [F]

        Returns:
            :obj:`str` or :obj:`None`
        """
    
        if type(obj) == steam.mesh.Mesh:
            db_dict = self.meshes
            db_meth = self.get_mesh
        elif type(obj) == steam.solution.Solution:
            db_dict = self.solns
            db_meth = self.get_soln
        else:
            raise TypeError ("Only a Mesh or Solution can be compared!")

        # If I'm empty, then return
        if len(db_dict) <= 0:
            return None

        (keys,hashes) = zip(*[(key,obj2.hash) for (key,obj2) in db_dict.items()])

        if obj.hash in hashes:
            index = hashes.index(obj.hash)
            key   = keys[index]
            ##! Check to make sure we don't have a hash conflict
            ##! For the record, I have never seen this happen AMS (6/20/18)
            if compare_obj:
                obj2 = db_meth(key)
                if obj != obj2:
                    raise AssertionError ("Hash conflict error when comparing objects!")
            return key
        return None
 

    def add_mesh(self, mesh, meshid=None):
        """ Append a mesh to the list of meshes and give unique GridID. 

        Args:
            mesh (:obj:`steam.mesh`): The mesh to store in database.
            meshid (:obj:`str`,optional) MeshID to use when storing
                the mesh.  Otherwise, unique default will be found.
        Returns:
            :obj:`str`: The MeshID that was used to store the mesh.

        """

        ### Update hash
        mesh.update_hash()

        ### Check to see if this mesh already exists
        old_meshid = self.find_obj_hash(mesh)
        if old_meshid is not None:
            return old_meshid

        ### Get a unique id for the mesh
        mesh_num = len(self.meshes)
        if meshid is None:
            meshid = "{:s}{:06d}".format("g",mesh_num)

        ### If it exists, increment the mesh_num until it doesn't conflict
        while meshid in self.meshes:
            mesh_num += 1
            meshid    = "{:s}{:06d}".format("g",mesh_num)

        # Add a copy of the mesh
        mesh = steam.mesh.Mesh(copy=mesh)

        ### Load into memory and/or drop to disk.  Meshes are kept in memory
        ### and on disk, so either way we keep it in self.meshes
        if self.onDisk:
            self.onDisk_check_open()
            mesh.write_hdf5(self.hdf5,root="/meshes/mesh_"+meshid)

        self.meshes[meshid] = mesh

        return meshid
 
    def update_meshid(self, meshid, mesh):
        """ Update the mesh at mesh_id in the database with mesh.

        This will replace the grid for all cases in the database that
        use this meshid.  For this reason, it can touch several solutions
        or the entire database.  It is recommended that this method is
        reserved for cosmetic changes only (translation of a grid, rotation,
        or scaling).  Changing the dimensionality of the grid could result
        in solutions being incompatable with the new mesh.

        The id is reused and associated with the new mesh.  If the new mesh
        is already in the database as another meshid, then the existing mesh's
        meshid is used instead of the one passed into this method.  The meshid
        passed into this method is removed from the database in that instance 
        and all solutions as associated with the existing mesh's meshid.

        All solutions that use the meshid passed into this method are updated
        to point to the new mesh or, in the case of it being a duplicate of an
        existing mesh, the existing mesh and associated meshid.

        Args:
            meshid (:obj:`str`) : Existing meshid to replace.
            mesh (:obj:`~steam.mesh.Mesh`): The mesh to store in database.

        Returns:
            (:obj:`str`): meshid that cases with the existing meshid now have.  This could be different than what was passed in if the input duplicated an existing mesh.

        """

        if type(mesh) is not steam.mesh.Mesh:
            raise TypeError ("'mesh' input is not a STEAM Mesh: {}".format(type(mesh)))

        ### Update hash
        mesh.update_hash()

        ### Check to see if this meshid is legitimate
        if meshid not in self.meshes:
            raise LookupError ("Mesh id, {}, not found in database!".format(meshid))

        ### Check to see if this mesh is already the mesh in question
        if self.compare_hash(mesh,self.get_mesh(meshid)):
            return meshid

        old_meshid = self.find_obj_hash(mesh)
                
        # Store what solutions were pointing to this old mesh
        old_solns = list(self.data[self.data.meshid == meshid].solnid)

        if old_meshid is not None:
            logging.debug(" Updating mesh '{}' to point to duplicate, '{}'".format(meshid,old_meshid))
            # Replace everyone that is using the old mesh with the duplicate
            self.data.replace({'meshid':{meshid:old_meshid}},inplace=True)

            # Remove the old mesh
            self.remove_mesh(meshid)

            meshid = old_meshid
        else:
            logging.debug(" Updating mesh '{}' to point to new mesh".format(meshid))

            # Add a copy of the mesh
            mesh = steam.mesh.Mesh(copy=mesh)

            ### Load into memory and/or drop to disk.  Meshes are kept in memory
            ### and on disk, so either way we keep it in self.meshes
            if self.onDisk:
                self.onDisk_check_open()
                mesh.write_hdf5(self.hdf5,root="/meshes/mesh_"+meshid)

            self.meshes[meshid] = mesh

        # Point solutions to the new mesh
        for skey in old_solns:
            self.solns[skey].mesh = self.meshes[meshid]

        return meshid

 
    def update_meshidx(self, index, mesh):
        """ Update the mesh at the selected index in the database

        This will replace the grid for the selected case in the database.
        If you are not also updating the solution, then please reserve
        this method for cosmetic changes only (translation of a grid, rotation,
        or scaling).  Changing the dimensionality of the grid could result
        in the solution being incompatable with the new mesh.

        If the new mesh is already in the database, then the existing mesh's
        is used instead of adding a new one.  The mesh that was at this location
        is removed from the databse in that instance.

        Args:
            index (:obj:`int`) : Database index to update
            mesh (:obj:`~steam.mesh.Mesh`): The mesh to store in database.

        Returns:
            (:obj:`str`): meshid used at this index.

        """

        if type(mesh) is not steam.mesh.Mesh:
            raise TypeError ("'mesh' input is not a STEAM Mesh: {}".format(type(mesh)))

        ### Update hash
        mesh.update_hash()

        ### Check to see if this index is legitimate
        try:
            old_meshid = self.data.meshid[index]
        except:
            raise LookupError ("Could not find mesh at index {}, in database!".format(index))

        ### Check to see if this mesh is already the mesh in question
        if self.compare_hash(mesh,self.get_mesh(old_meshid)):
            return old_meshid

        meshid = self.add_mesh(mesh)
         
        # Update the database
        self.data.at[index,"meshid"] = meshid

        # Remove the old mesh.  If it's still in use somewhere then it will stick around
        self.remove_mesh(old_meshid)

        # Point solution to the new mesh
        self.solns[self.data.solnid[index]].mesh = self.meshes[meshid]

        return meshid

    def remove_mesh(self,meshid):
        """ Remove a mesh from the database if no cases use it.

        Checks to see if the meshid is used in the 'meshid' column of
        self.data.  If it isn't, then it removes it from the self.meshes[]
        dictionary.  If the meshid is still in use, then this does nothing.

        **Note:**  This does not clean up the hdf5 file so if you did this a lot,
        then you might want to write the database to a new container in the
        future.

        Args:
            meshid (:obj:`str`): ID of the mesh to be removed.

        Returns:
            (:obj:`bool`): Flag indicating whether or not the mesh was actually removed.
        """

        # Is the meshid used?
        used = (self.data.meshid == meshid).any()

        if not used:
            del self.meshes[meshid]
            return True
        else:
            logging.debug(" Meshid called in .remove_mesh, '{}', still in use!".format(meshid))

        return False

    def add_soln(self, soln):
        """ Simple routine to append a soln to the dictionary
        of solns and take care of some bookkeeping. The uniqure solutionID 
        is returned.

        Args:
            soln (:obj:`~steam.solution.Solution`): The solution to add.

        Returns:
            :obj:`str`: The solnID that was used to store the
                solution.
        """

        ## First, update the hash
        soln.update_hash()

        ### Check to see if this soln already exists
        old_solnid = self.find_obj_hash(soln)
        if old_solnid is not None:
            return old_solnid

        ### Get a unique id for the mesh
        soln_num = len(self.solns)
        solnid = "{:s}{:06d}".format("s",soln_num)

        ### If it exists, increment the soln_num until it doesn't conflict
        while solnid in self.solns:
            soln_num += 1
            solnid    = "{:s}{:06d}".format("s",soln_num)

        # Add a copy of the solution
        soln = steam.solution.Solution(copy=soln)

        ### Load into memory or drop to disk
        if self.onDisk:
            
            self.onDisk_check_open()

            #print(" - Adding to disk: ",self.hdf5_root+"/solns/hash_"+solnid)
            soln.write_hdf5(
                           self.hdf5,
                           root=self.hdf5_root+"/solns/soln_"+solnid,
                           options={"ioMesh":False}
                           )
            soln.data = None

        self.solns[solnid] = soln

        return solnid


    def update_solnid(self, solnid, soln):
        """ Simple routine to replace a solution in the database.

        This will replace the soln for all cases in the database that
        use this solnid.  For this reason, it can touch several solutions
        or the entire database.  The database will be updated with the mesh 
        associated with the new solution in order to maintain consistency.

        The solnid is reused and associated with the new solution.  If the new
        solution is already in the database as another solnid, then the 
        existing solution's solnid is used instead of the one passed in.  
        The solnid passed into this method is removed from the databse in 
        that instance and all cases that use the old solnid are update to 
        point to the existing one.

        Args:
            solnid (:obj:`str`) : Existing solnid to replace.
            soln (:obj:`~steam.solution.Solution`): The solution to store in the database.

        Returns:
            (:obj:`str`): solnid that cases with the existing solnid now have.  This could be different than what was passed in if the input duplicated an existing solution.

        """

        if type(soln) is not steam.solution.Solution:
            raise TypeError ("'soln' input is not a STEAM Solution: {}".format(type(soln)))

        ## Update hash
        soln.update_hash()

        ### Check to see if this solnid is legitimate
        if solnid not in self.solns:
            raise LookupError ("Solution id, {}, not found in database!".format(solnid))
         
        ### Check to see if the soln is a duplicate of the existing one
        if self.compare_hash(soln,self.get_soln(solnid)):
            return solnid

        # Update the mesh to reflect what is stored in the database
        soln_ind = self.data.index[ self.data.solnid == solnid ]
        for ind in soln_ind:
            meshid = self.update_meshidx( ind, soln.mesh )
        soln.mesh = self.meshes[ meshid ]
        #soln.mesh = self.solns[solnid].mesh
 
        old_solnid = self.find_obj_hash(soln)

        if old_solnid is not None:
            logging.debug(" Updating soln '{}' to point to duplicate, '{}'".format(solnid,old_solnid))
            # Replace everyone that is using the old soln with the duplicate
            self.data.replace({'solnid':{solnid:old_solnid}},inplace=True)

            # Remove the old soln
            self.remove_soln(solnid)

            solnid = old_solnid
        else:
            logging.debug(" Updating soln '{}' to point to new soln".format(solnid))

            # Add a copy of the solution
            soln = steam.solution.Solution(copy=soln)

            ### Load into memory or drop to disk
            if self.onDisk:
                
                self.onDisk_check_open()
                soln.write_hdf5(
                               self.hdf5,
                               root=self.hdf5_root+"/solns/soln_"+solnid,
                               options={"ioMesh":False}
                               )
                soln.data = None

            self.solns[solnid] = soln

            logging.debug(" Successfully updated soln '{}' to point to new soln".format(solnid))
        return solnid

 
    def update_solnidx(self, index, soln):
        """ Update the soln at the selected index in the database

        This will replace the solution for the selected case in the database.
        The database will be updated with the mesh associated with the new 
        solution in order to maintain consistency.

        If the new solution is already in the database, then the existing 
        solution is used instead of adding a new one.  The solution that was 
        at this location is removed from the database in that instance.

        Args:
            index (:obj:`int`) : Database index to update
            soln (:obj:`~steam.solution.Solution`): The solution to store in database.

        Returns:
            (:obj:`str`): solnid used at this index.

        """

        if type(soln) is not steam.solution.Solution:
            raise TypeError ("'soln' input is not a STEAM Solution: {}".format(type(soln)))

        ### Update hash
        soln.update_hash()

        ### Check to see if this index is legitimate
        try:
            old_solnid = self.data.solnid[index]
        except:
            raise LookupError ("Could not find soln at index {}, in database!".format(index))

        ### Check to see if this soln is already the soln in question
        if self.compare_hash(soln,self.get_soln(old_solnid)):
            return old_solnid

        # Update the mesh to point to the database's mesh
        meshid = self.update_meshidx( index, soln.mesh )
        soln.mesh = self.meshes[ meshid ]

        solnid = self.add_soln(soln)
         
        # Update the database
        self.data.at[index,"solnid"] = solnid

        # Remove the old soln.  If it's still in use somewhere then it will stick around
        self.remove_soln(old_solnid)

        return solnid


    def remove_soln(self,solnid):
        """ Remove a soln from the database if no cases use it.

        Checks to see if the solnid is used in the 'solnid' column of
        self.data.  If it isn't, then it removes it from the self.solns[]
        dictionary.  If the solnid is still in use, then this does nothing.

        **Note:**  This does not clean up the hdf5 file so if you did this a lot,
        then you might want to write the database to a new container in the
        future.

        Args:
            solnid(:obj:`str`): ID for the solution to be removed
        """

        # Is the solnid used?
        used = (self.data.solnid == solnid).any()

        if not used:
            del self.solns[solnid]
        else:
            logging.debug(" Solnid called in .remove_soln, '{}', still in use!".format(solnid))

        return


    def add_point(self, point):
        """ Add DBPoint to database.

        Check to see if the mesh already exists before adding it.
        solnID and meshID are made unique and updated.
        
        TODO: If independant variables are defined, then check to make
        sure that this is not a duplicate.

        Args:
            point (:obj:`steam.database.DBPoint`):
                DBPoint to add or a list of DBPoint objects to add.
        """

        ### Allow for the 'point' input to be a list of DBPoints
        if hasattr( point, '__iter__' ):
            for pt in point:
                self.add_point( pt )
            return

        ### Drop some legacy items that might still be around from another database
        drop = ['meshid', 'solnid',
                'soln_type', 'soln_path']

        for key in drop:
            if key in point.data.keys():
                point.data.drop([key],inplace=True)

        ### If we don't have any data, then init the frame with this
        if self.data is None:
            i = 0
            point.data.name = i
            if type(point.data) is pd.DataFrame:
                self.data = pd.DataFrame(point.data)
            else:
                ### Because of the way Pandas does things, need to 
                ### transpose a series.  Dumb dumb dumb.
                self.data = point.data.to_frame().T
            
            #! Check to see what datatypes things should be.
            #!  In verson 0.22 or newer we could use .infer_types(), but we have
            #!  to do it manually here.  If we don't do this, then everything is
            #!  type 'object' unless they were all one type to begin with.
            self.data = self.data.apply(pd.to_numeric, errors='ignore')

        else:
            ### Index for the next item (since it's 0 based)
            i = self.data.shape[0]

            ### Add the data
            self.data = pd.concat( (self.data, 
                                   pd.DataFrame(point.data.astype(float)).T),
                                   ignore_index=True )

        ### Add the mesh
        if (point.mesh is not None):
            meshid = self.add_mesh(point.mesh)
            self.data.at[i,"meshid"] = meshid
            if point.soln is not None:
                ### If the mesh already existed, rereference it instead of
                ### creating another instance
                point.soln.mesh = self.meshes[meshid]

        ### Append the solution
        if (point.soln is not None):
            solnid = self.add_soln(point.soln)
            self.data.at[i,"solnid"] = solnid

    def remove_pt( self, point ):
        """ Remove a point or list of points from the database.

        Args:
            point (:obj:`~steam.database.DBPoint` or :obj:`list`): The point (or a list of points) to be removed.
        """
        
        ### Accommodate point being either a DBPoint or a list of DBPoints
        points = np.asarray( point )
        if points.ndim == 0:
            points = points[None] # Makes points 1-dimensional

        assert all([pt.get_database() is self for pt in points]), ( 
               'Attempting to remove a DBPoint from a Database that ' +
               'does not contain it.' )

        ### Collect the mesh and solution IDs for points to be removed.
        ### These may or may not be removed, depending on whether
        ### other points us them
        indices = [pt.get_index() for pt in points]
        #drop_rows = self.data.loc[indices]
        self.data.drop( indices, axis='rows', inplace=True )
        meshids = [pt.get_meshid() for pt in points]
        solnids = [pt.get_solnid() for pt in points]

        ### Remove meshes from database if they aren't used anymore
        for meshid in set(meshids):
            self.remove_mesh( meshid )

        ### Remove solutions from database if they aren't used anymore
        for solnid in set(solnids):
            self.remove_soln( solnid )

    def get_point(self, index):
        """ Return a DBPoint at the specified data index.

        This is a copy of the original mesh and solution.

        Args:
            index (:obj:`int`): The (zero-based) index in data.

        Returns:
            :obj:`steam.database.DBPoint`:
                DBPoint of the specified index.
        """

        point = DBPoint()

        point.set_database(self)
        point.set_index(index)

        ### Get the data
        point.data = deepcopy(self.data.loc[index])

        ### Try and get the mesh
        try:
            point.mesh = self.get_mesh(point.data['meshid'])
            point.set_meshid(point.data['meshid'])
        except:
            point.mesh = None

        ### Try and get the soln
        try:
            point.soln = self.get_soln(point.data['solnid'])
            point.set_solnid(point.data['solnid'])
        except:
            point.soln = None
    
        return point

    ### ## HDF5 Methods

    def write_hdf5(self,hdf5,root="/",options=None):
        """ Simple routine to write all of this database to disk.
        
        Args:
            hdf5 (`~pandas.io.pytables.HDFStore`): HDFStore object.
            root (:obj:`str`,optional): Location in HDF5 to write database.
                Default is to write to root location "/".

        Current options are:

           * ioMesh  = True/False : Should the meshes be writen?    [T]
           * ioSoln  = True/False : Should the solutions be writen? [T]
           * ioSpace = True/False : Should the intspaces be writen? [T]
           * onDisk    = True/False : Orverwrites the database's onDisk attribute.

        """
        #! Store hdf5 for later look-up from disk
        self.hdf5      = hdf5
        self.hdf5_mode = hdf5._mode
        self.hdf5_root = root

        #! Process options
        ioMesh  = True
        ioSoln  = True
        ioSpace = True
        if options is not None:
            if ('ioMesh' in options):
                ioMesh = options['ioMesh']
            if ('ioSoln' in options):
                ioSoln = options['ioSoln']
            if ('ioSpace' in options):
                ioSpace = options['ioSpace']
            if ('onDisk' in options):
                self.onDisk = options['onDisk']

        ### Data is always stored in memory and on disk, so it needs
        ### to be updated since it might have changed.  Call the
        ### method from the table parent.
        if (self.data is not None):
            logging.info(" - Writing data...")
            super().write_hdf5(hdf5,root=root)

        ### Save the interpolation spaces
        if self.intspace and ioSpace:
            logging.info(" - Writing IntSpace(s)...")
            for (i,spaceid) in enumerate(self.intspace.keys()):
                steam.util.progress_bar(i+1,len(self.intspace.keys()),prefix="Writing Spaces:")
                space = self.intspace[spaceid]
                space_key = steam.util.string_to_hdf5_key(spaceid)
                space.write_hdf5(
                               hdf5,
                               root=root+"/spaces/"+space_key
                               )
            ### Write out the list of them, too
            spaceids = list(self.intspace)
            steam.util.write_hdf5_pnode(
                               hdf5,
                               root+"/space_ids",
                               spaceids
                               )


        ### Mesh data is always stored in memory and on disk, so it needs
        ### to be updated since it might have changed.
        meshids = []
        if self.meshes and ioMesh:
            logging.info(" - Writing mesh(s)...")
            for (i,meshid) in enumerate(self.meshes.keys()):
                steam.util.progress_bar(i+1,len(self.meshes.keys()),prefix="Writing Meshes:")
                mesh = self.meshes[meshid]
                mesh.write_hdf5(
                               hdf5,
                               root=root+"/meshes/mesh_"+meshid
                               )
            ### Write out the list of them, too
            meshids = list(self.meshes)
            steam.util.write_hdf5_pnode(
                               hdf5,
                               root+"/mesh_ids",
                               meshids
                               )

        ### If this it held on-disk, then it should have been manipulated
        ### onDisk.  Otherwise, go ahead and write the maps out.
        solnids = []
        if self.solns and ioSoln:
            if (not self.onDisk):
                logging.info(" - Writing solution map(s)...")
                for (i,solnid) in enumerate(self.solns.keys()):
                    steam.util.progress_bar(i+1,len(self.solns.keys()),prefix="Writing Solns :")
                    soln = self.solns[solnid]
                    soln.write_hdf5(
                                   hdf5,
                                   root=root+"/solns/soln_"+solnid,
                                   options={"ioMesh":False}
                                   )
            ### Write out the list of them, too
            solnids = list(self.solns)
            steam.util.write_hdf5_pnode(
                               hdf5,
                               root+"/soln_ids",
                               solnids
                               )

        ### Need to turn a lot of meta data into something
        ### that we can store in an HDF5.
        ### Attributes are an option.  I think that there is
        ### only 64kb allowed on each group node, so it might
        ### be wise to make dummy "meta","var","tri" groups
        logging.info(" - Writing meta data information...")

        ### Make objects that will store attributes
        info = pd.Series(np.zeros(1))
        iroot = root+'/info' 
        hdf5.put(iroot,info)
        hdf5.get_storer(iroot).attrs.description  = self.description
        hdf5.get_storer(iroot).attrs.indep        = self.indep
        hdf5.get_storer(iroot).attrs.dep          = self.dep

        return

    def read_hdf5(self,hdf5,root="/",options=None):
        """ Simple routine to read all necessary items from disk.
        
        Args:
            hdf5 (`~pandas.io.pytables.HDFStore`): HDFStore object.
            root (:obj:`str`,optional): Location in HDF5 to write database.
                Default is to write to root location "/".
            options (:obj:`dict`, optional): options to vary behavior.

        Current options are:

            * ioMesh  = True/False : Should the meshes be loaded?    [T]
            * ioSoln  = True/False : Should the solutions be loaded? [T]
            * ioSpace = True/False : Should the intspaces be loaded? [T]
            * onDisk  = True/False : Orverwrites the database's onDisk attribute.

        """
        #! Store hdf5 for later look-up from disk
        self.hdf5      = hdf5
        self.hdf5_mode = hdf5._mode
        self.hdf5_root = root

        #! Process options
        ioMesh  = True
        ioSoln  = True
        ioSpace = True
        if options is not None:
            if ('ioMesh' in options):
                ioMesh = options['ioMesh']
            if ('ioSoln' in options):
                ioSoln = options['ioSoln']
            if ('ioSpace' in options):
                ioSpace = options['ioSpace']
            if ('onDisk' in options):
                self.onDisk = options['onDisk']

        #! Location for info
        iroot = root+'/info' 

        ### Data table is always stored in memory and on disk, so it needs
        ### to be updated since it might have changed.
        exist = hdf5.get_node(root+"/data")
        if not exist is None:
            logging.info(" - Reading Data")
            super().read_hdf5(hdf5,root=root)


        ### read the interpolation space data
        try:
            spaceids = steam.util.read_hdf5_pnode(
                               hdf5,
                               root+"/space_ids"
                               )
        except:
            spaceids = []

        if (ioSpace and len(spaceids) > 0):
            for (i,spaceid) in enumerate(spaceids):
                steam.util.progress_bar(i+1,len(spaceids),prefix=" - Reading Spaces:")
                self.intspace[spaceid] = steam.interpolate.IntSpace()
                space_key = steam.util.string_to_hdf5_key(spaceid)
                self.intspace[spaceid].read_hdf5(
                    hdf5,
                    root=root+"/spaces/"+space_key,
                    )


        ### Mesh data is always stored in memory and on disk, so it needs
        ### to be updated since it might have changed.
        try:
            meshids = steam.util.read_hdf5_pnode(
                               hdf5,
                               root+"/mesh_ids"
                               )
        except:
            meshids = []

        if (ioMesh and len(meshids) > 0):
            for (i,meshid) in enumerate(meshids):
                steam.util.progress_bar(i+1,len(meshids),prefix=" - Reading Meshes:")
                self.meshes[meshid] = steam.mesh.Mesh()
                self.meshes[meshid].read_hdf5(
                    hdf5,
                    root=root+"/meshes/mesh_"+meshid
                    )

        ### If I'm not running onDisk, then I need to load all of this
        ### into memory.  Otherwise, opening the file handle is ehough.
        try:
            solnids = steam.util.read_hdf5_pnode(
                               hdf5,
                               root+"/soln_ids"
                               )
        except:
            solnids = []

        #! If there aren't any solutions, then skip.  If this is
        #! an onDisk load, then don't read the data here.
        ioSolnData = True
        if self.onDisk:
            ioSolnData = False

        if (ioSoln and len(solnids) > 0):
            for (i,solnid) in enumerate(solnids):
                steam.util.progress_bar(i+1,len(solnids),prefix=" - Reading Solns :")
                self.solns[solnid] = steam.solution.Solution()
                self.solns[solnid].read_hdf5(
                    hdf5,
                    root=root+"/solns/soln_"+solnid,
                    options={"ioMesh":False,"readData":ioSolnData}
                    )

            #! Re-associated the meshes with each solution
            for (i,row) in self.data.iterrows():
                solnid = row['solnid']
                meshid = row['meshid']
                mesh = self.meshes[meshid]
                soln = self.solns[ solnid]
                soln.mesh = mesh

        ### Need to turn a lot of meta data into something
        ### that we can store in an HDF5.
        ### Attributes are an option.  I think that there is
        ### only 64kb allowed on each group node, so it might
        ### be wise to make dummy "meta","var","tri" groups
        logging.info(" - Reading Metadata")

        ### The info object stores attributes
        try:
            self.description     = hdf5.get_storer(iroot).attrs.description   
        except:
            # This is what it used to be.  Should be removed
            logging.info(" This is an out-of-date database w/ meta flag")
            self.descritpion     = hdf5.get_storer(iroot).attrs.meta
        self.indep           = hdf5.get_storer(iroot).attrs.indep  
        self.dep             = hdf5.get_storer(iroot).attrs.dep    

        return

    ### ## Intepolation Interfaces

    def new_intspace(self,key='entire',subset=None,
                     indep=None,
                     subvars=None,subvals=None,
                     scale=None,opts=None):
        """ Create new interpolation space inside of database.
        
        When using a subset, a list of Pandas.Index values can be used.  If more
        convenient, a Pandas DataFrame query string can also be used.

        Examples:

        Use a subset list :    ```subset=list(db.data[db.data.Mach > 0.5].index)```
        Use a subset query:    ```subset='Mach > 0.5'```
        Use a subset query:    ```subset='Mach > 0.5 & AMCT == 10 & Alpha > 2.0'```

        Perform interpolation in 3-space: (p1^3, p2, log10(p3)): ```opts={'p1':"**3",'p3':"log10"}```

        Args:
            key (:obj:`str`,optional): Name of space. Defaults to "entire".
            subset (:obj:`str`, :obj:`~pandas.Index`, or :obj:`list`), optional) : string, pandas.Index or list of index locations to use. Default to all.
            indep (:obj:`list`, optional): independant parameters to use. Default to database's.
            scale (:obj:`str` or :obj:`dict`, optional): Scale parameter to pass to interpspace.  Defaults to None.
            opts (:obj:`dict`, optional): Operations to perform on the indep parameters.  See example above.
        """

        if subset is None:
            ### Cam ran into a problem where the intspace was built with 
            ### indices that were strings.  That was really hard to debug.
            sub = list(self.data.index.map(int))
        elif type(subset) == str:
            ### This is a pandas query call
            sub = list(self.data.query(subset).index.map(int))
        else:
            sub = subset

        if indep is None:
            param = self.indep
        else:
            param = indep

        #! Get a matrix of the point data that we're going to use

        self.intspace[key] = steam.interpolate.IntSpace(
                indep=param,
                name=key,
                subvars=subvars,
                subvals=subvals,
                )
        # Get the indicies and np array for interpolation space
        sub_idx = sub
        # This looks at only the index locations specified above
        # and only uses the columns from the indep parameters.
        points  = self.data.loc[sub_idx][param].values

        self.intspace[key].build_space(
                points    = points,
                out_index = sub_idx,
                scale=scale,
                opts=opts
                )
        

    def interpolate(self,  cases,  space='entire',
                           cols=None, soln_cols=None,
                           interp_solns=True, subset=None,
                           full_stencil=False, use_dep=True,
                           quiet=False, return_empty=True,
                           DANGER_MODE=False):
        """Method to use the Delaunay triangulation to interpolate
        the dependent variables.  Uses linear interpolation via
        calculated barycentric weights.
        
        Args: 
            cases (:obj:`DataFrame`): Indep parameters with values as interpolation target. Can be a dataframe, dictionary, or STEAM table/database.
            space (:obj:`str`): Name of interpolation space. Optional ('entire').
            cols (:obj:`list`,optional): list of columns to interpolate.  
                Default is dep vars if use_dep is True, otherwise [].
            soln_cols(:obj:`list`,optional) : list of solution map variables to interpolate.
                Default is all solution variables.
            subset (:obj:`list`): Subset of elements to interpolate to.
            full_stencil (:obj:`bool`): Keep all indicies in stencil (instead of only >1e-10 ones) [F]
            use_dep (:obj:`bool`): Include database dependent variables in interpolated variables [T]
            quiet (:obj:`bool`): Do not print out-of-database warnings. [F]
            return_empty (:obj:`bool`): Return empty DBPoints when attempting to interpolate to points outside of the interpolation space. [T]
            DANGER_MODE (:obj:`bool`): "Danger Mode" is faster to run, though less forgiving.  soln_cols must be explicitly defined in Danger Mode.

        Returns:
            :obj:`steam.database.DBPoint`:
                DBPoint populated with appropriate information.
        """

        if DANGER_MODE:
            case_list = cases
        else:
            ### If the user doesn't provide specified variables, then
            ### assume all of them. Provided, we're interpolating them.

            ### Get the names of all arguments and their current values
            import inspect
            frame = inspect.currentframe()
            args, _, _, vals = inspect.getargvalues(frame)

            # Turn whatever I was given into a DataFrame
            case_list = None
            if isinstance(cases, dict):
                case_list = pd.DataFrame([cases])
            elif isinstance(cases, pd.DataFrame):
                case_list = cases
            elif (hasattr(cases,'data')):
                case_list = cases.data

            ### I can't set 'cols=[]' as a keyword argument since that would save the
            ### values with subsequent calls.  Google "Python mutable default arguments".
            if cols is None:
                cols = list()

            ### Include database dependent variables in interpolation
            if use_dep:
                cols.extend( self.dep )

            ### Remove duplicate varibles
            cols = list(set(cols))

            ### If there aren't any maps, then we're not interpolating them
            if self.solns is None:
                interp_solns = False
            elif len(self.solns) == 0:
                interp_solns = False

            ### NOTE: Eventually, we may want to find the intersecting set of
            ###       variables at each interpolated point, which would
            ###       allow interpolation between solutions without common
            ###       sets of variables.  Obviously this would be slow, though.
            if soln_cols is None and interp_solns == True:
                soln_cols = self.get_soln(self.data.solnid.values[0]).vars()

            # ### Print the DANGER_MODE input that would yield the same output
            # command = ( 'steam.database.Database.interpolate( <Database>, ' +
            #         '<DataFrame containing target cases>, DANGER_MODE=True' )
            # ditch_args = ['self', 'cases', 'DANGER_MODE']
            # out_args = [arg for arg in args if arg not in ditch_args]
            # local = dict(globals(), **locals()) # Necessary for vals to work
            # for arg in out_args:
            #     if isinstance(vals[arg], str):
            #         vals[arg] = '"{:s}"'.format( vals[arg] )
            #     command += ', {:s}={:s}'.format( arg, str(vals[arg]) )
            # command += ' )'
            # command = command.replace( '\n      ', '' )
            # logger.info( 'To achieve same output more quickly using "Danger'
            #              + ' Mode", use this command:\n{:s}'.format(command) )

        #############################
        ### Now that all parameters are established, perform interpolation
        #############################
        #! Get the verticies and weights
        intspace = self.intspace[space]
        (vert_mat,weight_mat) = intspace.get_weights(case_list)

        DBPoints = []

        for i, (index, case) in enumerate( case_list.iterrows() ):

            #! Start interpolating
            verts   =   vert_mat[i]
            weights = weight_mat[i]

            ### We're going to pack it all into a DBPoint for output
            outP = steam.database.DBPoint()

            stencil_string = ""
            string = "   {:1s} {:^5s} {:^13s}".format("i","Ind","Weight")
            for var in intspace.indep:
                string += " {:^10s}".format(var)
            
            stencil_string += string + "\n"

            string = "   {:1s} {:^5s} {:^13s}".format('-','---','TARGET')
            for var in case[intspace.indep]:
                string += " {:10f}".format(var)
            stencil_string += string + "\n"


            for (j,p) in enumerate(verts):
                # Remove verticies that aren't really used to clean up output
                if (not full_stencil) and (abs(weights[j]) < 1e-10):
                        continue
                string = "   {:1d} {:5d} {:13.5e}".format(j,p,weights[j])
                for var in self.data[intspace.indep].loc[p].values:
                    string += " {:10f}".format(var)
                stencil_string += string + "\n"

            outP.set_stencil(stencil_string)

            #! Are we inside the database?  Look for a negative weight.
            #! We will allow negative weights as long as those values are
            #! much smaller than the maximum weight.
            ### NOTE: The 4e-8 threshold was determined by Cam's work 
            ###       comparing STEAM with CBAERO.  We used the smallest
            ###       threshold that would not cause any errors.
            if (min(weights) < 0) and abs(min(weights) / max(weights)) > 4e-8:

                if not quiet:
                    logger.warning("Case outside of database!  Returning no data!")

                conds = "Conds:"
                case_dict = case[intspace.indep].to_dict()
                for key in sorted(case_dict):
                    value = case_dict[key]
                    conds += " '{}'={:.3f}".format(key,value)

                if not quiet:
                    logger.warning(conds)

                string = "   {:1s} {:^5s} {:^13s}".format('!','OUT',
                                                          'OF DATABASE!')
                # Append to the stencil string
                outP.set_stencil(outP.get_stencil() + string + "\n")

                if return_empty:
                    DBPoints.append(outP)
                continue

            #!  Put the input case data (indep plus whatever else)
            #!  and the interpolated dependent variables together.
            #!
            #!  Interpolate the dependent variables
            outP.data = pd.Series(data=case)

            if (len(cols) > 0):
                interp_dep = 0
                for (j,p) in enumerate(verts):
                    if (not full_stencil) and (weights[j] < 1e-10):
                            continue
                    interp_dep += self.data[cols].loc[p]*weights[j]
                outP.data = interp_dep.combine_first( outP.data )

            ### Interpolate soln variables
            soln = None
            mesh = None
            if interp_solns:
            
                #! since shapes and the mesh should be consistent across everything
                # Get the mesh and also check a solution to get row size
                mesh_id = self.data['meshid'][verts[0]]
                soln_id = self.data['solnid'][verts[0]]
                if subset is None:
                    mesh = self.meshes[mesh_id]
                else:
                    mesh = self.meshes[mesh_id].return_subset(subset)

                soln_shape = self.get_soln_df(self.data['solnid'][verts[0]],
                                              cols=soln_cols,
                                              index=subset)
                soln_index = soln_shape.index
                soln_shape = soln_shape.shape

                soln       = steam.solution.Solution(mesh=mesh)

                soln.data  = pd.DataFrame(np.zeros(soln_shape),
                                          columns=soln_cols,
                                          index=soln_index)

                for (k,p) in enumerate(verts):
                    if (not full_stencil) and (weights[k] < 1e-10):
                        continue
                    old_soln  = self.get_soln_df(self.data['solnid'][p],
                                                 cols=soln_cols,
                                                 index=subset)
                    #soln.data += old_soln*weights[k]
                    ###### DEBUG
                    try:
                        soln.data += old_soln*weights[k]
                    except TypeError:
                        logger.warning( 'Unable to multiply weights.  There may'
                             + ' be a non-number variable multiplied here.' )
                        raise
                    ###### END DEBUG
            
                #print(" - new_soln:\n")
                #print(soln.iloc[:5])

                soln.init()

            outP.soln = soln
            outP.mesh = mesh

            DBPoints.append(outP)

        return (DBPoints)

    def indep_range( self, show = True ):
        """ Find the range of the database's independent parameters and 
        print it to the screen if desired (default).

        Args:
            show (:obj:`bool`): Print the variable ranges to screen? [T]
        Returns:
            :obj:`dict`: Dictionary containing independent variable names as keys and their minimum and maximum values in a tuple as values
        """

        indep_dict = {}

        ### Print the output header if desired
        if show:
            print('Independent parameter ranges:\n' )
            print('{:<30s} {:>12s} {:>12s}'.format('VarName', 
                                                   'MinVal', 'MaxVal'))
            print( ''.join( ['-']*72 ) )

        for var in self.indep:
            minval = min(self.data[var])
            maxval = max(self.data[var])
            indep_dict[ var ] = (minval, maxval)
                        
            if show:
                print( '{:<30s}: {:12g} {:12g}'.format(var, minval, maxval) )

        return indep_dict

    def to_common_mesh( self, target_mesh, nproc=None, inplace=True, k=None ):
        """ Given a mesh, convert all solutions in the database to that mesh.

        This method can run in serial or in parallel on a user-specified 
        number of processors.  
        
        The transformation(s) built by this method use the 
        :obj:`steam.interpolate.Transform.inverse_distance` method.  Users can
        specify the number of nearest points to use with the "k" keyword 
        argument.

        WARNING: This method runs very slowly, even on many processors.

        Args:
            target_mesh (:obj:`~steam.mesh.Mesh`): Target mesh for conversion
            nproc (:obj:`int`): Number of processors on which to run.  Default is to run on all available.
            inplace (:obj:`bool`): Alter the database in place.  Defaults to True.
            k (:obj:`int`): Number of nearest points to use for inverse_distance transformation.  Defaults to the default value in :obj:`steam.interpolate.Transform.inverse_distance`
        """
    
        ### Call the serial version if parallel is not needed
        if nproc==1:
            return self.to_common_mesh_serial( target_mesh, inplace=inplace, 
                                               k=k )

        logger.info( 'Database.to_common_mesh: ' + 
                     'Converting solutions to common mesh.')

        try:
            steam.libmesh._lmp
        except AttributeError:
            ### LibMesh is not initialized, which is what we want
            pass
        else:
            raise RuntimeError( 'Cannot use ' 
                       + 'steam.database.Database.to_common_mesh when libmesh'
                       + ' has been initialized.  This is due to conflicts '
                       + "between LibMesh and the multiprocessing module." )

        if not inplace:
            raise NotImplementedError( 
                              'Currently only modifies database in place.' )
        else:
            logger.debug('Database.to_common_mesh modifying database in place')

        from multiprocessing import Pool, Manager
        man = Manager()

        ### Get the mesh_id for the target mesh.  add_mesh checks if the new
        ### mesh is already stored in the database, so we can simply compare
        ### the ids 
        target_id = self.add_mesh( target_mesh )

        ### Build dict of solnids and meshids.
        solnid_2_meshid = man.dict( 
                              zip(self.data['solnid'], self.data['meshid']) )
        ### Build list of solnids that need to be transformed
        solns_2_run = man.list( [solnid for solnid, meshid in 
                            solnid_2_meshid.items() if meshid != target_id] )

        ### Check whether any Transforms were created
        if len( solns_2_run) == 0:
            logger.info('  --No conversions necessary' )
            return target_id

        ### Create all necessary transformations
        transforms = man.dict()
        for mesh_id, mesh in self.meshes.items():
            if mesh_id == target_id:
                continue
            if (mesh_id not in transforms):
                logger.debug(
                       " Making transformation from {:s} to {:s}".format(
                       mesh_id, target_id ) )
                trans = steam.interpolate.Transform(mesh, target_mesh)
                if k is None:
                    trans.inverse_distance()
                else:
                    trans.inverse_distance( k )
                transforms[mesh_id] = trans

        solnid_2_meshid = man.dict( 
                              zip(self.data['solnid'], self.data['meshid']) )
#        solns_2_run = man.list( [solnid for solnid, meshid in solnid_2_meshid
#                                 if meshid != target_id] )

        ### Apply transformations to each solution
        logger.debug(" Database.to_common_mesh:    Applying transformations")

        ### Create map-able function from convert_soln
        from functools import partial
        mapfunc = partial( _convert_soln, db=self, trans_dict=transforms, 
                           id_dict=solnid_2_meshid )

        ### Initialize progress_bar
        n_iter = 1
        total_iter = len(solns_2_run)
        with Pool( nproc ) as pool:
            for soln_tuple in pool.imap_unordered(mapfunc, solns_2_run):
                steam.util.progress_bar( n_iter, total_iter, 
                                         'Updating solutions:' )
                n_iter += 1

                solnid, new_soln = soln_tuple
                self.update_solnid( solnid, new_soln )

            logger.debug(" Database.to_common_mesh:    All solns updated")
        logger.debug(" Database.to_common_mesh:   Successfully completed")
        
        ### Return the meshid for the target_mesh
        return target_id

    def to_common_mesh_serial( self, target_mesh, inplace=True, k=None ):
        """ Given a mesh, convert all solutions in the database to that mesh.

        Unlike to_common_mesh, this method can only run on one processor.
        It therefore does not require the complications of the multiprocessing
        module and is much simpler.  When to_common_mesh is called with
        nproc=1, it calls this method.
        
        The transformation(s) built by this method use the 
        :obj:`steam.interpolate.Transform.inverse_distance` method.  Users can
        specify the number of nearest points to use with the "k" keyword 
        argument.

        Args:
            target_mesh (:obj:`~steam.mesh.Mesh`): Target mesh for conversion
            inplace (:obj:`bool`): Alter the database in place.  Defaults to True.
            k (:obj:`int`): Number of nearest points to use for inverse_distance transformation.  Defaults to the default value in :obj:`steam.interpolate.Transform.inverse_distance`
        """

        logger.info( 'Database.to_common_mesh_serial: ' + 
                     'Converting solutions to common mesh.')

        if not inplace:
            raise NotImplementedError( 
                              'Currently only modifies database in place.' )
        else:
            logger.debug(
                  'Database.to_common_mesh_serial modifying database in place')

        ### Get the mesh_id for the target mesh.  add_mesh checks if the new
        ### mesh is already stored in the database, so we can simply compare
        ### the ids 
        target_id = self.add_mesh( target_mesh )

        ### Build dict of solnids and meshids.
        solnid_2_meshid = dict( zip(self.data['solnid'], self.data['meshid']) )
        ### Build list of solnids that need to be transformed
        solns_2_run = [ solnid for solnid, meshid in 
                        solnid_2_meshid.items() if meshid != target_id ]

        ### Check whether any Transforms were created
        if len( solns_2_run) == 0:
            logger.info('  --No conversions necessary' )
            return target_id

        ### Create all necessary transformations
        transforms = {}
        for mesh_id, mesh in self.meshes.items():
            if (mesh_id not in transforms):
                logger.debug(
                       " Making transformation from {:s} to {:s}".format(
                       mesh_id, target_id ) )
                trans = steam.interpolate.Transform(mesh, target_mesh)
                if k is None:
                    trans.inverse_distance()
                else:
                    trans.inverse_distance( k )
                transforms[mesh_id] = trans

        ### Apply transformations to each solution
        logger.debug(" Database.to_common_mesh_serial:    Applying transformations")

        ### Create map-able function from convert_soln
        from functools import partial
        mapfunc = partial( _convert_soln, db=self, trans_dict=transforms, 
                           id_dict=solnid_2_meshid )

        ### Initialize progress_bar
        n_iter = 1
        total_iter = len(solns_2_run)
        for solnid in solns_2_run:
            steam.util.progress_bar( n_iter, total_iter, 
                                     'Updating solutions:' )
            meshid = solnid_2_meshid[ solnid ]
            transform = transforms[ meshid ]
            new_soln = transform.apply( self.get_soln( solnid ) )

            n_iter += 1
            self.update_solnid( solnid, new_soln )

        logger.debug(" Database.to_common_mesh_serial:  All solns updated")
        
        ### Return the meshid for the target_mesh
        return target_id

###################################
def _convert_soln( solnid, db, trans_dict, id_dict ):
    """ Given a solution ID and the necessary mappings between meshes
    and transformations, convert all solutions to the target mesh.

    NOTE:  This function is intended to be used EXCLUSIVELY by 
           Database.to_common_mesh.  It was originally defined within that
           method, but it must be a module-level function (as opposed to 
           a method or a function within a method) in order to be pickled
           and run in parallel with multiprocessing.Pool().

    Args:
        solnid (:obj:`str`): The ID for the solution to be updated.
        db (:obj:`~steam.database.Database`): The database in which solutions are being updated.
        trans_dict (:obj:`dict`): A dictionary of mesh transforms.  Keys are meshid strings and values are :obj:`~steam.interpolate.Transform` objects to project the solution onto the target grid.
        id_dict (:obj:`dict`): A dictionary with solnid strings as keys and the associated meshid strings as values.
    """

    old_meshid = id_dict[ solnid ]
    transform = trans_dict[ old_meshid ]
    new_soln = transform.apply( db.get_soln( solnid ) )

    return (solnid, new_soln)
###################################

class DBPoint():
    """Class definition for DBPoint

    A Database point encapsulates all of the data, mesh and
    solution variables from an STEAM database into a portable
    unit that can be passed between Databases or other classes
    and routines.

    It is comprised of a Pandas series, solution map (DataFrame),
    and a mesh object.
    """

    def __init__(self,
                 data=pd.Series(dtype=np.float64),
                 mesh=None,
                 soln=None,
                 ):
        """Constructor for DBPoint class."""

        #! If it's not a series, then convert it to one
        if (type(data) is not pd.Series):
            in_data = pd.Series(data)
        else:
            in_data = data

        # We deepcopy the data to allow users to pass in a subset of
        # their database and have a unique condition, but shared soln/mesh
        # - Examples include fake points to hold last
        self.data  = deepcopy(in_data)  # Series - Indep/Dep variables
        self.mesh  = mesh               # Mesh Object
        self.soln  = soln               # Dictionary of map DataFrames

        # The output from interpolation is a DBPoint and in that case,
        # it populates the interpolation stencil.  This is set in
        # database.interpolate()
        self._stencil = None

        # These are optional attributes that only show up if this point
        # was from a database.  They help us get back to that database
        # and the mesh/soln we represent.  These are set in
        # database.get_point().
        self._database = None
        self._meshid   = ""
        self._solnid   = ""
        self._index    = -1

        if (mesh is None and soln is not None):
            self.mesh = self.soln.mesh

    def __str__(self):
        """Return the print string of the DBPoint object."""

        string  = "\nDatabase Point:\n"
        string += " Data:\n"
        for key in self.data.keys():
            string += "   {:15} : {}\n".format(key,self.data[key])

        if self._stencil is not None:
            string += " Stencil:\n"
            string += self._stencil

        if self._database is not None:
            string += " Database Location:\n"
            string += "   Meshid : {}:\n".format(self._meshid)
            string += "   Solnid : {}:\n".format(self._solnid)
            string += "   Index  : {}:\n".format(self._index)

        import re

        string += " Soln:\n"
        if (self.soln is None):
            string +="   None\n";
        else:
            sol_str = str(self.soln)
            sol_str = re.sub( '^','   '  ,sol_str)
            sol_str = re.sub('\n','\n   ',sol_str)
            sol_str = re.sub('\n +$','\n',sol_str)
            #! Remove some lines we don't want
            sol_str = re.sub('^.*Object.*\n','',sol_str) 
            string += sol_str

        string += " Mesh:\n"
        if (self.mesh is None):
            string +="   None\n";
        else:
            grid_str = str(self.mesh)
            grid_str = re.sub( '^','   '  ,grid_str)
            grid_str = re.sub('\n','\n   ',grid_str)
            #! Remove some lines we don't want
            grid_str = re.sub('^.*Mesh Object.*\n','',grid_str) 
            string += grid_str

        return string
    
    def __repr__(self):
        return str(self)

    def set_stencil(self,stencil=""):
        """ Set the interpolation stencil string."""
        self._stencil = stencil

    def get_stencil(self):
        """ Get the interpolation stencil string."""
        return self._stencil

    def set_database(self,database=None):
        """ Set the source database."""
        self._database = database

    def get_database(self):
        """ Get the source database."""
        return self._database

    def set_index(self,index=-1):
        """ Set the source database index."""
        self._index = index

    def get_index(self):
        """ Get the source database index."""
        return self._index

    def set_meshid(self,meshid=""):
        """ Set the source database meshid."""
        self._meshid = meshid

    def get_meshid(self,meshid=""):
        """ Get the source database meshid."""
        return self._meshid

    def set_solnid(self,solnid=""):
        """ Set the source database solnid."""
        self._solnid = solnid

    def get_solnid(self,solnid=""):
        """ Set the source database solnid."""
        return self._solnid

    def update_soln(self, soln = None):
        """ Update the solution this point represents in the database it came 
        from.

        This keys off of the ._database and ._index attributes that were set 
        when this point was created.

        The ._solnid leaving this subroutine will be updated to reflect the 
        new solnid of the solution.

        Args:
            soln (:obj:`~/steam.solution.Solution`, optional): Solution with which to replace the solution currently associated with this DBPoint.  If this input is populated, the mesh associated with it will be used to update the DBPoint mesh for consistency.
        """

        ### Consistently set mesh
        if soln is not None:
            self.soln = soln
            self.mesh = soln.mesh

        if self._database is None:
            raise ValueError("No database associated with this DBPoint!")

        new_solnid = self._database.update_solnidx(self._index,self.soln)
        self.set_solnid(new_solnid)

    def update_mesh(self, mesh = None):
        """ Update the mesh this point represents in the database it came from.

        This keys off of the ._database and ._index attributes that were set when
        this point was created.

        The ._meshid leaving this subroutine will be updated to reflect the new meshid
        of the mesh.

        Args:
            mesh (:obj:`~/steam.mesh.Mesh`, optional): Mesh with which to replace the mesh currently associated with this DBPoint.  If this input is populated, the solution associated with it will have its mesh updated for consistency.

        """

        ### Consistently set mesh
        if mesh is not None:
            self.mesh = mesh
            self.soln.mesh = mesh

        if self._database is None:
            raise ValueError("No database associated with this DBPoint!")

        new_meshid = self._database.update_meshidx(self._index,self.mesh)
        self.set_meshid(new_meshid)


#def DBSolution(steam.solution.Solution):
#    """Class definition for STEAM DBSolution
#
#    This inherits from the STEAM Solution class.
#
#    This is designed to be an augmented solution that is aware of
#    its place in a database.  The purpose is to allow for users
#    to update and modify solutions in a database in a controlled
#    manner.
#    
#    For databases that are onDisk, this context aware solution
#    ensures that the data on disk is updated as necessary.  It
#    also ensures that solution hashes can be maintained by calls
#    to .commit().
#
#    When a solution is 'opened' then its open attribute is set.  This
#    ensures that when it is closed a commit is automatically called.
#    If a DBSolution is not opened and falls out of scope the changes
#    are lost.
#
#    Several methods create DBSolutions:
#        - When iterating through a database the provided solution is
#          a DBSolution.  It is not opened, so the user must manually
#          make a .commit() call to save changes.  Alternatively, they
#          could call the .open() attribute such that it falling out of
#          scope automatically calls a .commit().
#        - If the user calls a database's .open_soln method then an
#          opened DBSolution is returned and automatically calls
#          .commit() when it falls out of scope.
#        - More ...?
#
#    The mesh associated with the solution is not opened, regardless of
#    how the DBSolution was created.  It is included as a DBMesh and has
#    the same properties as the DBSolution.  It is never automatically
#    .commit()ed unless it has been .open()ed itself (in which case its
#    .close() method calls its .commit()) or unless it is specifically
#    called to .commit() (e.g. `dbsoln.mesh.commit()`).
#
#    Args:
#        active (:obj:`bool`): Should the solution be opened and updated on delete? [F]
#        mesh   (:obj:`~steam.database.DBMesh()`): The mesh associated with this solution.
#    """
#
#    def __init__(self,active=False,mesh=None):
#        """Constructor for DBSolution Class
#
#        """
#
#        raise NotImplementedError("More work to be done")
#
#        ### Run the Table class __init__ method
#        super().__init__(self,mesh=mesh)
#
#        self.solnid     = ""            # The ID in the database
#        self.database   = None          # The database I came from
#        self.opened     = active        # Was I opened?
#
#    def update():
#        """ Check DBSolution back into database.
#
#        This is going to recalculate the solution hash,
#        replace the data in the database, and close the
#        DBSolution.
#        """
#        
#        self.update_hash()
#        # Put it back into the database
#        self.database.update_soln(self.solnid,self.as_solution())
#        self.opened = False
#
#    def as_solution():
#        return steam.solution.Solution(copy=self)
#
