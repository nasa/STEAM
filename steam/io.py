""" Module for the Mesh/Solution IO

For read statements, return a mesh and soln (if availible and requested).

For write statements, write mesh and soln if requested.

"""
import logging
logger = logging.getLogger(__name__)
logger.debug("Top of file")

import pandas as pd
import numpy  as np
import steam.util

from code import interact
#interact( local = dict( globals(), **locals() ) )

################################################################################
#! ### Read Statements


def read_pickle(file,mesh=None,soln=None):
    """ Read mesh/solution from a compressed pickle file.
    
    This is most useful for writing to disk between various operations.

    Args:
        file (:obj:`string`): Path to file.
        mesh (:obj:`mesh`)      : Optional, mesh object to overwrite
        soln (:obj:`solution`)  : Optional, solution object to overwrite

    Returns:
        tuple( `mesh`,`soln`): mesh and solution (or `None` if either not read)
    """

    import pickle
    import zlib

    logger.info("Reading pickle file '{}'".format(file))

    f = open(file,"rb")
    payload = pickle.loads(zlib.decompress(f.read()))
    f.close()

    pmesh = payload['mesh']
    psoln = payload['soln']
    if (psoln is not None):
        psoln.mesh = pmesh # I reassociate since this was stripped during write

    if pmesh is not None:
        if mesh is None:
            mesh = steam.mesh.Mesh()
        mesh.copy(pmesh)

    if psoln is not None:
        if (soln is None):
            soln = steam.solution.Solution()
        soln.copy(psoln)

    return (mesh,soln)

def read_tecplot_ascii(file,ioMesh=True,ioSoln=True,mesh=None,soln=None,
                    toElement=True,collapseNodes=True,pointcloud=False):
    """ Read mesh/solution from a tecplot ASCII file.
    
    Args:
        file (:obj:`string`): Path to file.
        ioMesh (:obj:`bool`)  : Return mesh? [T]
        ioSoln (:obj:`bool`)  : Return solution (if exists)? [T]
        mesh (:obj:`mesh`)      : Optional, mesh object to overwrite
        soln (:obj:`solution`)  : Optional, solution object to overwrite
        toElement (:obj:`bool`) : Convert to element-based data? [T]
        collapseNodes (:obj:`bool`) : Collapse coincident nodes? [T]
        pointcloud (:obj:`bool') : File has point cloud data (True) vs Finite Element

    Returns:
        tuple( `mesh`,`soln`): mesh and solution (or `None` if either not read)
    """

    if steam.has_dataio:
        import pydata

        debug = False
    
        if (ioMesh is False and ioSoln is False):
            return (mesh,soln)
    
        if not pointcloud:
            # Read the tecplot file into a tecio object
            logger.info("Reading ASCII tecplot file '{}'".format(file))
            tec = pydata.data.PyData ()
            tec.read_file (file, 'tec')
    
            # If quads, then split to tris
            if tec.cells == 4: 
                logger.info("Splitting quad cells into tris")
                tec.split_quads_to_tris ()
    
            # Lets collapse the nodes just in case (may not need)
            if collapseNodes: 
                logger.info("Collapsing duplicate nodes (removing the common nodes)")
                tec.collapse_common_nodes ()
    
            # Populate the mesh object from the pydata object
            if ioMesh:
                if (mesh is None):
                    mesh = steam.mesh.Mesh()
                mesh.xyz_pt = pd.DataFrame (data = tec.nodes.transpose (), columns=['X', 'Y', 'Z'])
                mesh.conn   = pd.DataFrame (data = np.array (tec.connectivity, dtype=np.int32))
                mesh.comp   = pd.DataFrame (data = np.array (tec.components, dtype=np.int32))
    
                if (toElement):
                    mesh.get_xyz_el()
    
            # Populate the soln object from the pydata object
            if ioSoln:
                if (soln is None):
                    soln = steam.solution.Solution()
                soln.mesh = mesh
    
                # Determine variable locations (node, ele, mixed, none)
                if tec.soln_vars is not None:
                    if   len ([var for var in tec.soln_node if var is None]) == len (tec.soln_vars) and \
                         len ([var for var in tec.soln_ele  if var is None]) == 0: varloc = 'ele'
                    elif len ([var for var in tec.soln_ele  if var is None]) == len (tec.soln_vars) and \
                         len ([var for var in tec.soln_node if var is None]) == 0: varloc = 'node'
                    elif len ([var for var in tec.soln_node if var is None]) > 0 and \
                         len ([var for var in tec.soln_ele  if var is None]) > 0: varloc = 'mixed'
                    else: varloc = None
                else: varloc = None
    
                # Data at elements, leave it there
                if varloc == 'ele':
                    soln_dictionary = {var: tec.soln_ele [iv] for iv, var in enumerate (tec.soln_vars)}
                    soln.data  = pd.DataFrame.from_dict(soln_dictionary)
    
                # Data at nodes and toElement=True, move the data
                if varloc == 'node' and toElement:
                    soln_dictionary = {var: tec.soln_node [iv] for iv, var in enumerate (tec.soln_vars)}
                    soln.data  = pd.DataFrame.from_dict(soln_dictionary)
                    logger.info("Moving node based data to elements with steam.solution method 'node_to_element'")
                    soln.node_to_element()
    
                # Data at mixed and toElement=True, move only the data that needs to the element
                if varloc == 'mixed' and toElement:
                    logger.info("Moving node based data to elements with pydata.data.PyData method 'soln_to_elem'")
                    tec.soln_to_elem ()
                    soln_dictionary = {var: tec.soln_ele [iv] for iv, var in enumerate (tec.soln_vars)}
                    soln.data  = pd.DataFrame.from_dict(soln_dictionary)
    
                # Data at mixed and toElement=False, cannot do this, through warning
                if varloc == 'mixed' and not toElement:
                    logger.info("Cannot leave mixed (node based, element based) data mixed, must be moved to element with toElement=True option")
    
                # Data at node and toElement=False
                if varloc == 'node' and not toElement:
                    soln_dictionary = {var: tec.soln_node [iv] for iv, var in enumerate (tec.soln_vars)}
                    soln.data  = pd.DataFrame.from_dict(soln_dictionary)
        else:
            # Read the tecplot file into a tecio object
            logger.info("Reading ASCII tecplot point cloud file '{}'".format(file))
            logger.info("Method assumes the first three variables are the mesh variables.")
            tec = pydata.tecio.TecFile ()
            tec.read_file (file)

            if tec.zones is not None:
                if len (tec.zones) > 0:
                    # Populate the mesh object from the pydata object
                    if ioMesh:
                        if (mesh is None):
                            mesh = steam.mesh.Mesh()
                            df_mesh = None
                    if ioSoln:
                        if (soln is None):
                            soln = steam.solution.Solution()
                            df_soln = None

                    zone_dict = {}
                    for i, zone in enumerate(tec.zones):
                        if ioMesh:
                            comp_num = i + 1
                            df = pd.DataFrame (data = zone.nodes.transpose (), columns=['X', 'Y', 'Z'])
                            df['comp_num'] = comp_num
                            if df_mesh is not None: df_mesh = df_mesh.append (df, ignore_index=True)
                            else: df_mesh = df
                            ### Save mesh indices and zone titles for later
                            if zone.title not in zone_dict:
                                zone_dict[zone.title] = comp_num

                        if ioSoln and len (tec.vars) > 3:
                            soln_dictionary = {var: zone.soln [iv] for iv, var in enumerate (tec.vars [3:])}
                            df = pd.DataFrame.from_dict(soln_dictionary)
                            if df_soln is not None: df_soln = df_soln.append (df, ignore_index=True)
                            else: df_soln = df

                    if ioMesh:
                        mesh.xyz_pt = df_mesh[['X', 'Y', 'Z']]
                        mesh.to_point_cloud ()
                        ### Set the component names
                        mesh.comp = pd.DataFrame( df_mesh.comp_num.values )
                        for name in zone_dict.keys():
                            mesh.set_comp( name, str(zone_dict[name]) )
                        #interact( local = dict( globals(), **locals() ) )
                        
                    if ioSoln and len (tec.vars) > 3:
                        soln.data  = df_soln
                        soln.mesh = mesh
                            
        return (mesh,soln)
    else:
        logger.warn("Reading Tecplot ASCII files requires JSC's PyData module in PYTHONPATH")
        logger.warn("Failed to read '{}'".format(file))
        raise Exception("\n * Cannot read Tecplot file:\n" +
                          "   requires JSC's PyData module!"
                        )

def read_triq_ascii(file,ioMesh=True,ioSoln=True,mesh=None,soln=None,
                    toElement=True):
    """ Read ASCII tri or triq file.
    
    Args:
        file (:obj:`string`): Path to file.
        ioMesh (:obj:`bool`)  : Return mesh? [T]
        ioSoln (:obj:`bool`)  : Return solution (if exists)? [T]
        mesh (:obj:`mesh`)      : Optional, mesh object to overwrite
        soln (:obj:`solution`)  : Optional, solution object to overwrite
        toElement (:obj:`bool`) : Convert to element-based data? [T]

    Returns:
        tuple( `mesh`,`soln`): mesh and solution (or `None` if either not read)
    """

    debug = False

    if (ioMesh is False and ioSoln is False):
        return (mesh,soln)

    logger.info("Reading ASCII triq file '{}'".format(file))
    npts = 0
    nele = 0
    nq   = 0
    nodes = {'X' : [], 'Y' : [], 'Z' : []}
    elem  = {}
    comp  = {}
    qs    = {}
    f = open (file, 'r')

    #! First line is npt, ntri, (nq)
    for i in range(1):
        #! Read, trim, and split
        #rawline = ((f.readline()).strip()).split()
        #line  = rawline.strip()
        #split = line.split()
        split = ((f.readline()).strip()).split()

        npts = int(split[0])
        nele = int(split[1])
        if (len(split) > 2):
            nq   = int(split[2])
            # Set-up the dictionary
            for j in range(nq):
                num = j + 1
                qs["q{}".format(num)] = []
        if (debug):
            print(i,split)

    if (nq <= 0 and ioSoln):
        print(" - No q variables found, skipping solution")
        ioSoln=False

    #! Read in the nodes
    for i in range(npts):
        #! Read, trim, and split
        split = ((f.readline()).strip()).split()
        nodes['X'].append(float(split[0]))
        nodes['Y'].append(float(split[1]))
        nodes['Z'].append(float(split[2]))
        if (debug and i < 10):
            print(i,split)

    #! Read in the element connectivity
    for i in range(nele):
        #! Read, trim, and split
        split = ((f.readline()).strip()).split()
        #! Need to make it zero-based
        elem[i]    = [
                     int(split[0]) - 1,
                     int(split[1]) - 1,
                     int(split[2]) - 1
                     ]
        if (debug and i < 10):
            print(i,split,'->',elem[i])

    #! Read in the element components
    for i in range(nele):
        #! Read, trim, and split
        split = ((f.readline()).strip()).split()
        comp[i] = int(split[0])
        if (debug and i < 10):
            print(i,split,'->',comp[i])

    #! Read in the solution data
    for i in range(npts):
        #! Read, trim, and split
        split = ((f.readline()).strip()).split()
        #! In case the values span multiple lines, we need to
        #! keep track and append
        all_qs = split
        while (len(all_qs) < nq):
            split = ((f.readline()).strip()).split()
            all_qs = all_qs + split

        if (len(all_qs) != nq):
            print(nq,all_qs)
            print("* nq for vertex varies from global count!")
            raise IOError

        for j in range(nq):
            num = j + 1
            qs["q{}".format(num)].append(float(all_qs[j]))

    f.close()

    if ioMesh:
        if (mesh is None):
            mesh = steam.mesh.Mesh()
        mesh.xyz_pt = pd.DataFrame.from_dict (nodes)
        mesh.conn   = pd.DataFrame.from_dict(elem, 'index')
        mesh.comp   = pd.DataFrame.from_dict(comp, 'index')
        if (toElement):
            mesh.get_xyz_el()

    if ioSoln:
        if (soln is None):
            soln = steam.solution.Solution()
        soln.mesh  = mesh
        soln.data  = pd.DataFrame.from_dict(qs)
        soln.point = True
        if (toElement):
            soln.node_to_element()
        soln.init()

    return (mesh,soln)



def read_triq_uform(file,ioMesh=True,ioSoln=True,mesh=None,soln=None,
                    toElement=True,dp=True):
    """ Read fortran-unformated tri or triq file.

    This will check endianness, but it assumes single precision!
    
    Args:
        file (:obj:`string`): Path to file.
        ioMesh (:obj:`bool`)    : Return mesh? [T]
        ioSoln (:obj:`bool`)    : Return solution (if exists)? [T]
        mesh (:obj:`mesh`)      : Optional, mesh object to overwrite
        soln (:obj:`solution`)  : Optional, solution object to overwrite
        toElement (:obj:`bool`) : Convert to element-based data? [T]
        dp (:obj:`bool`)        : Double-precision floating point data? [T]

    Returns:
        tuple( `mesh`,`soln`): mesh and solution (or `None` if either not read)
    """

    from scipy.io import FortranFile
    import numpy as np

    debug = False

    if (ioMesh is False and ioSoln is False):
        return (mesh,soln)

    logger.info("Reading unformatted triq file '{}'".format(file))
    npts = 0
    nele = 0
    nq   = 0

    #! Floating precision
    dp_type  = 'f4'
    int_type = 'i4'
    int_sizes = [8,12]  # I have no idea
    if dp:
        dp_type = 'f8'
        int_type = 'i8'
        int_sizes = [16,24]  # I have no idea

    #! Determine the endianess
    #! - The first record is either two or three ints
    #!   so it's going to be size 12 bytes or 8 bytes.
    (be,le) = (False,False)
    with open(file, 'rb') as f:
        size = np.fromfile(f, dtype='>i4', count=1)
        if size in int_sizes :
            be = True
    with open(file, 'rb') as f:
        size = np.fromfile(f, dtype='<i4', count=1)
        if size in int_sizes :
            le = True

    if (be == le):
        raise FloatingPointError("triq does not register as uniquely big or little endian!")
        
    #! Set byte order
    bo = "="
    if (be):
        bo = ">" # Big
        logger.debug("Big Endian file")
    if (le):
        bo = "<" # Little
        logger.debug("Little Endian file")
    logger.debug("Byte Order: "+bo)

    f = FortranFile(file,'r',header_dtype=bo+'u4')
    dim = f.read_ints(bo+int_type)
    npts = int(dim[0])
    nele = int(dim[1])
    if (len(dim) > 2):
        nq   = int(dim[2])

    if (nq <= 0 and ioSoln):
        print(" - No q variables found, skipping solution")
        ioSoln=False

    #! Read in the nodes
    tmp = f.read_reals(bo+dp_type)
    xyz = np.reshape(tmp,(npts,3))

    #! Read in the element connectivity
    tmp = f.read_reals(bo+int_type)
    conn= np.reshape(tmp,(nele,3))

    #! Read in the element components
    comp = f.read_reals(bo+int_type)

    #! Read the q information if any
    if (nq > 0):
        tmp = f.read_reals(bo+dp_type)
        q   = np.reshape(tmp,(npts,nq))

    f.close()
    
    #! Swap to little endian if I just read big endian
    if (be):
        xyz  = xyz.byteswap().newbyteorder()
        conn = conn.byteswap().newbyteorder()
        comp = comp.byteswap().newbyteorder()
        if (nq > 0):
            q = q.byteswap().newbyteorder()

    # Make an array of the q names   
    qvars = []
    for j in range(nq):
        num = j + 1
        qvars.append("q{}".format(num))


    if ioMesh:
        if (mesh is None):
            mesh = steam.mesh.Mesh()
        mesh.xyz_pt = pd.DataFrame(xyz, columns=["X","Y","Z"])
        mesh.conn   = pd.DataFrame(conn) - 1
        mesh.comp   = pd.DataFrame(comp)
        if (toElement):
            mesh.get_xyz_el()

    if ioSoln:
        if (soln is None):
            soln = steam.solution.Solution()
        soln.mesh  = mesh
        soln.data  = pd.DataFrame(q,columns=qvars)
        soln.point = True
        if (toElement):
            soln.node_to_element()
        soln.init()

    return (mesh,soln)

################################################################################
#! ### Write Statements

def write_pickle(mesh=None,soln=None,name="output.pkl",ioMesh=True,):
    """ Write mesh/solution as a compressed pickle file.
    
    A mesh will only write nodes and connectivity.
    A solution will write the mesh and solution variables, by default..
    This is most useful for writing to disk between various operations.

    Args:
        mesh (:obj:`steam.mesh.Mesh()`): Mesh to write
        soln (:obj:`steam.solution.Solution()`): Solution to write
        name (:obj:`string`): The output file name
        ioMesh (:obj:`bool`): Have the solution write the mesh? [T]
    """

    if (mesh is None and soln is None):
        return

    logger.info("Writing pickle file '{}'".format(name))

    # I can pass in a mesh. If not, use the solutions.  If the
    # solution doesn't have one, then don't write one
    if (mesh is None):
        mesh = soln.mesh
    if (mesh is None and ioMesh):
        print(" - I do not have a mesh, so this is only solution data.")
    if (ioMesh is not True):
        mesh = None

    import pickle
    import zlib

    f = open(name,"wb")

    payload = dict()
    payload['mesh'] = mesh
    payload['soln'] = None
    if (soln is not None):
        payload['soln'] = steam.solution.Solution(copy=soln)
        payload['soln'].mesh = None

    f.write(zlib.compress(pickle.dumps(payload)))
    f.close()

    return

def write_triq_uform(mesh=None,soln=None,name="output.triq",dp=True):
    """ Write mesh/solution as a fortran-unformated triq file.
    
    A mesh will only write nodes and connectivity.
    A solution will write the mesh and solution variables.

    If the data is element-based, then convert to nodes first.

    Args:
        mesh (:obj:`steam.mesh.Mesh()`): Mesh to write
        soln (:obj:`steam.solution.Solution()`): Solution to write
        name (:obj:`string`): The output file name
        dp (:obj:`bool`)        : Double-precision point data? [T]
    """

    if (mesh is None and soln is None):
        return

    # I can pass in a mesh. If not, use the solutions.  If the
    # solution doesn't have one, then freak out.
    if (mesh is None):
        mesh = soln.mesh
    if (mesh is None):
        raise Exception("\n * Cannot write triq file:\n" +
                          "   No mesh in solution and no mesh passed in!"
                        )

    logger.info("Writing unformated triq file '{}'".format(name))

    from scipy.io import FortranFile

    dp_type  = 'float32'
    int_type = 'int32'
    if dp:
        dp_type = 'float64'
        int_type = 'int64'

    f = FortranFile(name,'w')
    #! Write the dims
    npt  = mesh.xyz_pt.shape[0]
    nele = mesh.conn.shape[0]
    if (soln is not None):
        ### Set the store_at flag
        soln.set_store_flag()

        #! First check to make sure that things are node based:
        if (soln.store_at.upper().find( 'ELEM' ) != 0 and
            soln.store_at.upper().find( 'NODE' ) != 0):
            raise Exception("Solution needs to be node or element!")
#        if (soln.element or not soln.point):
        if soln.store_at.upper().find( 'ELEM' ) == 0:
            logger.info(" - Writing element-based data to nodes")
            soln = soln.element_to_node(return_new=True)
        nq   = soln.data.shape[1]
        dim  = np.array([npt,nele,nq],dtype=int_type)
    else:
        nq   = 0
        dim  = np.array([npt,nele],dtype=int_type)

    f.write_record(np.array(dim                , dtype=int_type ))
    f.write_record(np.array(mesh.xyz_pt.values , dtype=dp_type))
    f.write_record(np.array(mesh.conn.values + 1, dtype=int_type ))
    f.write_record(np.array(mesh.comp.values   , dtype=int_type ))
    if (nq > 0):
        f.write_record(np.array(soln.data.values,dtype=dp_type))

    f.close()

    return

def write_tecplot_ascii(mesh=None,soln=None,name="output.dat"):
    """ Write mesh/solution as an ASCII triq file.
    
    A mesh will only write nodes and connectivity.
    A solution will write the mesh and solution variables.

    Args:
        mesh (:obj:`steam.mesh.Mesh()`): Mesh to write
        soln (:obj:`steam.solution.Solution()`): Solution to write
        name (:obj:`string`): The output file name
    """

    if (mesh is None and soln is None):
        raise ValueError('No mesh or solution given; Nothing is going to happen.' )
        return

    # Use the mesh provided
    # If not provided then use the mesh associated with the solution 
    # If still no mesh then raise an exception
    if (mesh is None):
        mesh = soln.mesh
    if (mesh is None):
        raise Exception("\n * Cannot write Tecplot file:\n" +
                          "   No mesh in solution and no mesh passed in!"
                        )
    try:
        mesh.is_point_cloud
    except AttributeError:
        mesh.is_point_cloud = False


    if steam.has_dataio:
        import pydata
        
        if not mesh.is_point_cloud:

            if soln is not None:
                if soln.store_at is not 'NODES' and soln.store_at is not 'ELEMENTS':
                    soln.set_store_flag ()
        
            logger.info("Writing ASCII tecplot file '{}'".format(name))

            # Create pydata object
            tec = pydata.data.PyData ()

            # Populate the pydata object with mesh information
            tec.nodes        = mesh.xyz_pt.values.transpose ()
            tec.connectivity = mesh.conn.values.tolist ()
            tec.components   = [comp for sub in mesh.comp.values.tolist() 
                                for comp in sub]
            tec.xyz_vars     = mesh.xyz_pt.columns.values.tolist ()

            # Populate the pydata object with the solution information 
            if (soln is not None):
                tec.soln_vars = soln.data.columns.values.tolist ()

                # node based data
                if (soln.store_at == 'NODES'): 
                    tec.soln_node = [var  for var in 
                                     soln.data.values.transpose()]
                    tec.soln_ele  = [None for var in tec.soln_vars]

                # element based data
                elif (soln.store_at == 'ELEMENTS'):
                    tec.soln_ele  = [var for var in 
                                     soln.data.values.transpose()]
                    tec.soln_node = [None for var in tec.soln_vars]

            tec.write_file (name, 'tec')

        else:
            tec = pydata.tecio.TecFile ()
            tec.vars = mesh.xyz_pt.columns.values.tolist()
            if soln is not None:
                tec.vars += soln.data.columns.values.tolist ()
            tec.xvar = 0
            tec.yvar = 1
            tec.zvar = 2
            tec.zones = []

            #interact( local = dict( globals(), **locals() ) )
            for comp in set(np.squeeze(mesh.comp.values)):
                pt_subset = mesh.xyz_pt.loc[mesh.get_comp(comp)]
                zone = pydata.tecio.TecZone ()
                zone.title = str(comp)
                zone.zonetype = "ORDERED"
                zone.datapacking = "BLOCK"
                zone.nodes = pt_subset.values.transpose ()
                zone.i = pt_subset.shape [0]
                zone.j = 1
                zone.k = 1
                zone.varlocation = ['node', 'node', 'node']
                if soln is not None:
                    soln_subset = soln.data.loc[mesh.get_comp(comp)]
                    zone.varlocation += ['node' for var in range (soln_subset.shape [1])]
                    zone.soln = [var for var in soln_subset.as_matrix ().transpose ()]
                tec.zones.append (zone)

            tec.write_file (name)
    
    else:
        logger.info("Writing Tecplot ASCII files requires JSC's PyData module")
        logger.info("Failed to create file: '{}'".format(name))
        raise Exception("\n * Cannot write Tecplot file:\n" +
                          "   requires JSC's PyData module!"
                        )

    return

def write_triq_ascii(mesh=None,soln=None,name="output.triq"):
    """ Write mesh/solution as an ASCII triq file.
    
    A mesh will only write nodes and connectivity.
    A solution will write the mesh and solution variables.

    If the data is element-based, then convert to nodes first.

    Args:
        mesh (:obj:`steam.mesh.Mesh()`): Mesh to write
        soln (:obj:`steam.solution.Solution()`): Solution to write
        name (:obj:`string`): The output file name
    """

    if (mesh is None and soln is None):
        return

    # I can pass in a mesh. If not, use the solutions.  If the
    # solution doesn't have one, then freak out.
    if (mesh is None):
        mesh = soln.mesh
    if (mesh is None):
        raise Exception("\n * Cannot write triq file:\n" +
                          "   No mesh in solution and no mesh passed in!"
                        )

    logger.info("Writing ASCII tecplot file '{}'".format(name))

    if (soln is not None):
        ### Set the store_at flag
        soln.set_store_flag()

        #! First check to make sure that things are node based:
        #if (soln.element == soln.point):
        if (soln.store_at.upper().find( 'ELEM' ) != 0 and
            soln.store_at.upper().find( 'NODE' ) != 0):
            raise Exception("Solution needs to be node or element!")
        #if (soln.element or not soln.point):
        if soln.store_at.upper().find( 'ELEM' ) == 0:
            logger.info(" - Writing element-based data to nodes")
            soln = soln.element_to_node(return_new=True)
        nq   = soln.data.shape[1]
    else:
        nq = 0

    npts = mesh.xyz_pt.shape[0]
    nele = mesh.conn.shape[0]

    with open (name, 'w') as f:
        f.write("{} {}".format( npts, nele) )
        if (soln is not None):
            f.write(" {}".format( nq ) )
        f.write("\n")

        arg = {'header':False,
               'index':False,
               'sep':' ',
               'float_format':'%.16g'
               }

        ### Write points
        mesh.xyz_pt.to_csv( f, **arg)

        ### Write elements -- Add 1 to all elements because 
        ### they're indexed off 1 instead of 0
        (mesh.conn + 1).to_csv( f, **arg)

        ### Write components
        mesh.comp.to_csv( f, **arg)

        if (soln is not None):
            ### Write Solution data
            soln.data.to_csv( f, **arg)

    return

def write_cdat(mesh=None,soln=None,name="output.cdat",
               header=True,fmt='%.12e'):
    """ Write a mesh or solution as a VTU file.
    
    A mesh will only write element X,Y,Z values if the data is element-based
    and point X,Y,Z values if the data is point-based.  If there is not solution,
    then a mesh that contains element centroids will be written as element-based.
    A solution will write the mesh and solution variables.

    Args:
        mesh (:obj:`steam.mesh.Mesh()`): Mesh to write
        soln (:obj:`steam.solution.Solution()`): Solution to write
        name (:obj:`string`): The output file name
        header (:obj:`bool`): Print the variable name header?
        fmt (:obj:`string`): Format string according to numpy.savetxt(). ['%.12e']
    """

    if (mesh is None and soln is None):
        return

    # I can pass in a mesh. If not, use the solutions.  If the
    # solution doesn't have one, then freak out.
    if (mesh is None):
        mesh = soln.mesh
    if (mesh is None):
        raise Exception("\n * Cannot write cdat:\n" +
                          "   No mesh in solution and no mesh passed in!"
                        )

    logger.info("Writing cdat file '{}'".format(name))

    #! We want to write the mesh variables (X,Y,Z)
    vars = list(mesh.xyz_pt.columns)
    data = mesh.xyz_pt.values
    if soln is None:
        if (mesh.xyz_el is not None):
            data = mesh.xyz_el.values
    else:
        if soln.store_at.upper().find( 'ELEM' ) == 0:
            data = mesh.xyz_el.values

    #! and also any of the variables stored in the solution
    if (soln is not None):
        vars.extend(soln.data.columns)
        data = np.hstack((data,soln.data.values))

    head = ''
    if (header):
        head = "        ".join(vars)

    np.savetxt(name,data,header=head,fmt=fmt)

    return

def write_vtu(mesh=None,soln=None,name="output.vtu",writeComp=True,):
    """ Write a mesh or solution as a VTU file.
    
    A mesh will only write nodes and connectivity.
    A solution will write the mesh and solution variables.

    Args:
        mesh (:obj:`steam.mesh.Mesh()`): Mesh to write
        soln (:obj:`steam.solution.Solution()`): Solution to write
        name (:obj:`string`): The output file name
        writeComp : Export components? [T]
    """

    if (mesh is None and soln is None):
        return

    # I can pass in a mesh. If not, use the solutions.  If the
    # solution doesn't have one, then freak out.
    if (mesh is None):
        mesh = soln.mesh
    if (mesh is None):
        raise Exception("\n * Cannot write VTU:\n" +
                          "   No mesh in solution and no mesh passed in!"
                        )

    logger.info("Writing VTU file '{}'".format(name))

    #! First write the mesh
    import xml.dom.minidom

    doc = xml.dom.minidom.Document()
    root_element = doc.createElementNS("VTK", "VTKFile")
    root_element.setAttribute("type", "UnstructuredGrid")
    root_element.setAttribute("version", "0.1")
    root_element.setAttribute("byte_order", "LittleEndian")
    doc.appendChild(root_element)

    # Unstructured mesh element
    unstructuredGrid = doc.createElementNS("VTK", "UnstructuredGrid")
    root_element.appendChild(unstructuredGrid)

    # Piece 0 (only one)
    piece = doc.createElementNS("VTK", "Piece")
    piece.setAttribute("NumberOfPoints", str(mesh.xyz_pt.shape[0]))
    piece.setAttribute("NumberOfCells" , str(mesh.conn.shape[0]))
    unstructuredGrid.appendChild(piece)

    ### Points ####
    points = doc.createElementNS("VTK", "Points")
    piece.appendChild(points)

    # Point location data
    point_coords = doc.createElementNS("VTK", "DataArray")
    point_coords.setAttribute("type", "Float64")
    point_coords.setAttribute("format", "ascii")
    point_coords.setAttribute("NumberOfComponents", "3")
    points.appendChild(point_coords)

    if True:
        string = "\n"
        xyz = mesh.xyz_pt.values
        for i in range(xyz.shape[0]):
            point = xyz[i,:]
            string += " {:f} {:f} {:f}\n".format(point[0],point[1],point[2])
        point_coords_data = doc.createTextNode(string)
    else:
        # Attempt at binary
        from base64 import b64encode as b4e
        xyz    = mesh.xyz_pt.values.copy(order='C')
        string = ""
        for i in range(xyz.shape[0]):
            point = xyz[i,:]
            for j in range(3):
                string += b4e(point[j]).decode('ascii')

        #! 32-bit integer of length first
        string = b4e('{:032b}'.format(len(string)).encode()).decode('ascii') + string
        point_coords_data = doc.createTextNode("\n"+string)

    point_coords.appendChild(point_coords_data)

    #### Cells ####
    cells = doc.createElementNS("VTK", "Cells")
    piece.appendChild(cells)

    # Cell locations
    cell_connectivity = doc.createElementNS("VTK", "DataArray")
    cell_connectivity.setAttribute("type", "Int32")
    cell_connectivity.setAttribute("Name", "connectivity")
    cell_connectivity.setAttribute("format", "ascii")        
    cells.appendChild(cell_connectivity)

    # Cell connectivity data (Assumes triangles!
    string = "\n"
    conn = mesh.conn.values
    for i in range(conn.shape[0]):
        for j in range(len(conn[i,:])):
            string += " {:d}".format(conn[i,j])
       #string += "\n"
    connectivity = doc.createTextNode(string)
    cell_connectivity.appendChild(connectivity)

    cell_offsets = doc.createElementNS("VTK", "DataArray")
    cell_offsets.setAttribute("type", "Int32")
    cell_offsets.setAttribute("Name", "offsets")
    cell_offsets.setAttribute("format", "ascii")                
    cells.appendChild(cell_offsets)

    string1 = "\n"
    string2 = "\n"
    if (mesh.is_point_cloud):
        # Points
        off   = 1
        itype = 1
    else:
        # Tris
        off   = 3
        itype = 5

    for i in range(conn.shape[0]):
        string1 += " {:d}\n".format((i+1)*off) # Offset
        string2 += " {:d}\n".format(  itype  ) # Type
    offsets = doc.createTextNode(string1)
    cell_offsets.appendChild(offsets)

    cell_types = doc.createElementNS("VTK", "DataArray")
    cell_types.setAttribute("type", "UInt8")
    cell_types.setAttribute("Name", "types")
    cell_types.setAttribute("format", "ascii")                
    cells.appendChild(cell_types)
    types = doc.createTextNode(string2)
    cell_types.appendChild(types)

    #! Write out the component numbers.
    comp = mesh.comp.values

    type = "CellData"

    out_data = doc.createElementNS("VTK", type)
    piece.appendChild(out_data)

    if (writeComp):
        #! Loop over each variable
        data_id = doc.createElementNS("VTK", "DataArray")
        data_id.setAttribute("type", "Float64")
        data_id.setAttribute("Name", "Component")
        data_id.setAttribute("format", "ascii")                
        out_data.appendChild(data_id)

        string = "\n"
        for i in range(comp.shape[0]):
                string += " {:d}\n".format(comp[i,0])
        var_data = doc.createTextNode(string)
        data_id.appendChild(var_data)

    if (soln != None):

        data = soln.data.values

        #! If there is data and it's point-based, then close
        #! the cell-based data (components). If it's cell-based,
        #! then just keep going.
        type = None
        if soln.store_at == 'NODES':
            type = "PointData"
            out_data = doc.createElementNS("VTK", type)
            piece.appendChild(out_data)

        #! Loop over each variable
        for n in range(len(data[0,:])):
            var = soln.vars()[n]

            data_id = doc.createElementNS("VTK", "DataArray")
            data_id.setAttribute("type", "Float64")
            data_id.setAttribute("Name", var)
            data_id.setAttribute("format", "ascii")                
            out_data.appendChild(data_id)

            string = "\n"
            for i in range(data.shape[0]):
                    string += " {:f}\n".format(data[i,n])
            var_data = doc.createTextNode(string)
            data_id.appendChild(var_data)

    # Write to file and exit
    outFile = open(name, 'w')
#        xml.dom.ext.PrettyPrint(doc, file)
    doc.writexml(outFile, newl='\n')
    outFile.close()
