""" Module for the Container class
"""

import logging
logger = logging.getLogger(__name__)
logger.debug("Top of file")

import pandas as pd
import os
import steam

class Container():
    """Class definition for STEAM Container.

    The container handles packing and unpacking items into HDF5
    files.  It will keep track of the contents and save it into
    the HDF5 file.

    """

    def __init__(self,hdf_path,mode=None,desc=""):
        """Constructor for Container class.

        Several modes are availible for use.  These mimic those availible
        from Pandas/PyTables.  By default, existing files are opened in 'r' mode
        and paths that do not exist are opened in 'w' mode:

           * ``'r'`` Read-only; no data can be modified.
           * ``'w'`` Write; a new file is created (an existing file with the same name would be deleted).
           * ``'a'`` Append; an existing file is opened for reading and writing, and if the file does not exist it is created.
           * ``'r+'`` It is similar to ``'a'``, but the file must already exist.

        Args:
            hdf_path (:obj:`string`): Path to HDF5 file (to load).
            mode (:obj:`string`, optional): State to open file in [r].
            desc (:obj:`string`, optional): Description of container.

        """

        self.description = desc                   # Optional description for container
        self.toc         = {}       # The table of contents  | name -> (toc_path,type)
                                    # toc_path is the location in /toc/ for the object
        self.obj         = {}       # The contained objects  | name -> object
        self.hdf5_path   = hdf_path                      # The path to the HDF5 file if loaded

        self.to_delete   = []                            # List of items to delete when next write happens
        
        # Check to see if the file exists
        exists = os.path.isfile(hdf_path)

        # Set the default mode
        if mode is None:
            if exists:
                mode = 'r'
            else:
                mode = 'w'
        self.hdf5_mode   = mode                          # Mode to open the file 

        if exists:
            if  self.hdf5_mode == "w":
                logger.warning("File will be overwritten by container write: {}".format(hdf_path))
            else:
                self.read_toc(self.hdf5_path)
        else:
            if 'r' in self.hdf5_mode:
                raise KeyError("HDF5 file does not exist: '{}'".format(self.hdf5_path))

        logger.info(str(self))

    def __enter__(self):
        """Return self for a context."""
        return self

    def __exit__(self,exc_type, exc_val, exc_tb):
        """Exit the context.  Save if not 'r' mode."""

        #print(exc_type, exc_val, exc_tb)

        if not self.hdf5_mode == "r":
            self.write()

        # Gut the container so that it cannot be written by default:
        self.hdf5_path = None
        self.hdf5_mode = 'r'
        self.to_delete   = []
        return

    def __str__(self):
        """Return the print string of the Container object."""

        string  = "Container:\n {}\n".format(self.description)
        string += " HDF Path: {}\n".format(self.hdf5_path)
        string += " HDF Mode: {}\n".format(self.hdf5_mode)
        string += " TOC:\n"
        for name in sorted(self.toc.keys()):
            (path,type) = self.toc[name]
            string += "  {:30s} ({} , {})\n".format(name,path,type)
        return string
   
    def __repr__(self):
        return "STEAM Container Object"

    def open_store(self,path,mode='r'):
        """Uses Pandas to open an HDF5 Store and return it.

        File is opened and sets a default compression level.
        Returns the HDFStore object.
    
        Args:
            path (:obj:`string`, optional): Path to HDF5 file.
            mode : {'a', 'w', 'r', 'r+'}, default 'r'

        Returns:
            :obj:`~pandas.io.pytables.HDFStore`: Pandas HDFStore.

        """

        return pd.HDFStore(path,complevel=9,complib='blosc',mode=mode)


    def close_store(self,store):
        """ Closes a Pandas HDFStore and returns `None`.

        This returns a None so that it can be called similar to the open_store::

            hdf5 = open_store('foo.h5')
            hdf5 = close_store(hdf5)

        Args:
            store (:obj:`~pandas.io.pytables.HDFStore`): Pandas HDFStore.

        Returns:
            :obj:`None`: To mimic usage of :meth:`~steam.container.Container.open_store`.
        """

        store.close()
        return None

    def write_toc(self,hdf5,items=None,root="/"):
        """Store the meta data and TOC in the root group of the HDFStore.
        
        Args:
            hdf5 (:obj:`~pandas.io.pytables.HDFStore`): HDF5 file
            items (:obj:`list`): List of items to write [default: all]
            root (:obj:`string`,optional): Path inside of h5 to store under (Unused).
        
        """
 
        # Get the name of every object unless some were passed in
        if items is None:
            items = list(self.toc.keys())

        #! The TOC used to be written as a single entry into the HDF5; a pickled dict.
        #! However, I didn't want to have to read it and update it every time that something
        #! was changed/deleted in the container.  So now they are independant nodes with
        #! unique numbers (at the time the object was added).  This lets us remove or
        #! update each entry atomistically without having to reconcile two dictionaries.
        #!  - If an item is deleted, then we delete that TOC node when we delete the file.
        #!    Since the TOC is updated after objects are deleted/written, any subsequently
        #!    new object that was given that toc path will be written to that path afterwards.
        #!  - If an item has changed, then the act of writing it will update the TOC in memory
        #!    and the TOC will be written overwriting the existing node.
        #!  - When creating a new HDF5 file or writing only a subset of the items in the
        #!    container's memory, only those items' TOC entries will be written.
        logger.info("Writing Container TOC")
        for name in items:
            (tpath,ttype) = self.toc[name]
            steam.util.write_hdf5_pnode(hdf5, root+"/table_of_contents/"+tpath, (name,ttype))

        # Save out meta data
        #!  Make objects that will store attributes
        info = pd.Series([0])
        iroot = root+'/info'
        hdf5.put(iroot,info)
        hdf5.get_storer(iroot).attrs.description = self.description


    def read_toc(self,hdf5_path,root="/"):
        """Read the meta data and TOC from the root group of the HDFStore.

        Args:
            hdf5_path (:obj:`string`): HDF5 file
            path (:obj:`string`,optional): Path inside of h5 to read under (Unused).
        """

        #! This is dumb, but we need to run through the low level options
        #! in order to see what's in the table_of_contents group and read it
        import tables
        from tables.nodes import filenode
        #! I need a path since I can't just use the pandas HDF5store
        #! Open a filenode at /toc
        h5file = tables.open_file(hdf5_path,'r')

        logger.info("Reading Container TOC")
        try:
            items = h5file.list_nodes("/table_of_contents")
        except:
            h5file.close()
            raise

        hdf5   = self.open_store(hdf5_path,mode='r')
        node_names = [i._v_pathname for i in items]
        import os
        for node_name in node_names:
            tpath = os.path.basename(node_name)
            try:
                (name,ttype) = steam.util.read_hdf5_pnode(hdf5,root+"/table_of_contents/"+tpath)
            except:
                h5file.close()
                self.close_store(hdf5)
                raise

            self.toc[name] = (tpath,ttype)

        # Save out meta data
        iroot = root+'/info'
        try:
            self.description = hdf5.get_storer(iroot).attrs.description
        except:
            logger.info("No container description found")

        h5file.close()
        self.close_store(hdf5)

    def add_to_toc(self,name,obj):
        """ Add an object to the toc, just pass the name of the object.
        
        This will create a unique identifier for it to use when writting to a container."""

        i = 0
        toc_path = 'foo'
        # Fancy - check the first in the duple in a dicy of tuples
        while (i == 0 or toc_path in [v[0] for (k,v) in self.toc.items()]):
            i += 1
            toc_path = "toc_{}".format(i)

        obj_type = self.get_container_type(obj)

        path = steam.util.string_to_hdf5_key(toc_path)
        self.toc[name] = (path,obj_type)

    def get_container_type(self,obj):
        """ Return the type of the object.

        This returns either the actual type or 'pickle' based on if the object
        has a 'write_hdf5' method.
        """

        #! If it does not have a 'write_hdf5' method, then just pickle it
        write_op = getattr(obj, "write_hdf5", None)
        if not callable(write_op):
            obj_type = "pickle"
        else:
            obj_type = type(obj)

        return obj_type

    @steam.util.timer(logger)
    def read(self,items=None,options=dict()):
        """Read the contents of the container into memory.

        Options should contain keys equal to the object names. Each
        value is passed to the object type's loader as a keyword
        argument 'options'.  Example for a database object named 'database'::

           db_opts = {'onDisk':True,'ioMesh':True,'ioSoln':False}
           options = {'database':db_opts}
        
        Args:
            items (:obj:`list`, optional): List of items to load.
            options (:obj:`dict`, optional): options to pass to each item.

        """
    
        if self.hdf5_mode == 'w':
            raise IOError("Cannot read HDF5 in write mode")

        print("\nReading HDF5 File '{}'".format(self.hdf5_path))

        try:
            hdf5        = self.open_store(self.hdf5_path,mode=self.hdf5_mode)
        except:
            raise

        # Get the name of every object
        if items is None:
            items = list(self.toc.keys())

    
        #! Now go ahead and load each item
        for obj_name in sorted(items):

            #! Check to make sure that each item doesn't overlap with
            #! one that is already defined.  If it is, then raise exception.
            #! This shouldn't ever really happen.
            try:
                #! If I can exec the thing, then it exists, so error off.
                self.obj[obj_name]
                print("Variable {} already exists!  Cannot load container.".format(obj_name))
                hdf5        = self.close_store(hdf5)
                return
            except KeyError:
                #! This means it wasn't defined.  Good.
                #print(" - {} variable not found".format(var_name))
                pass

            #! Get the object and toc information          
            try:
                (tpath,ttype) = self.toc[obj_name]
            except KeyError:
                print(" !! '{}' object not found in container TOC!".format(obj_name))
                hdf5        = self.close_store(hdf5)
                return

            #!  - Make the path HDF5-safe
            root = "/"+steam.util.string_to_hdf5_key(obj_name)
            logger.info("Loading object '{}' from container".format(obj_name))
            #print(" Loading "+obj_name)

            if (ttype == "pickle"):
                node_obj = steam.util.read_hdf5_pnode(hdf5,root,compress=True)
            else:
                #! Make the requested type and load it
                node_obj = ttype()
                if (obj_name in options):
                    node_obj.read_hdf5(hdf5,root=root,options=options[obj_name])
                else:
                    node_obj.read_hdf5(hdf5,root=root)

            #! Set the variable inside this container
            self.obj[obj_name] = node_obj

        hdf5        = self.close_store(hdf5)

        ### Check whether the container was created using the same version
        ### of STEAM as the one reading it
            ### Check if the container was created with the same STEAM version
        version_node = 'steam_version'
        if version_node in self.toc:
            try:
                written_version = self.obj[version_node]
            except KeyError:
                ### steam_version wasn't read out of the container, likely
                ### just because the user was selective about reading
                ### instead of reading all the contents of the container
                logger.info( f'Container version unknown.  This probably ' + 
                    f'does not require your attention, though.  Go on living' +
                    f'your life and enjoy your time using STEAM.' )
            else:
                if written_version != steam.__version__:
                    logger.warning( f'Container written with version ' + 
                        f'{written_version} and reading with version ' +
                        f'{steam.__version__}.  This could cause problems.' )
        else:
            ### No version number in the container -- it must be old
            logger.warning( 
                f'Container written with an old version of STEAM.' )
        
    def write(self,path=None,mode=None,items=None,options=dict()):
        """Save the contents of the container to disk.

        If any remove commands have been executed, then the items on disk will
        be removed prior to writing.  The file will be adjusted in size to 
        clear space, too.  If a new path is identified, then nothing will be 
        deleted.  Further, 

        In "w" mode, this will overwrite the existing file (the first time) and
        then it will be changed to "a" mode for subsequent syncs.  Otherwise, 
        each successive write would overwrite the previous.  And if the 
        selection of items has changed, then it would obliterate the previous 
        ones that were written.

        Options should contain keys equal to the object names. Each
        value is passed to the object type's loader as a keyword
        arguement 'options'.  Example for a database object named 'database':

        ```
        db_opts = {'onDisk':False,'WriteMesh':True,'WriteSoln':False} 
        options = {'database':db_opts} 
        ```
        
        Args:
            path (:obj:`string`,optional): Path to .h5 file to save, if different from initial.
            mode (:obj:`string`,optional): Defaults to current mode.  If you have a new path that points to an existing file, you must explicitly set to 'w','a', or 'r+'.
            items (:obj:`list`, optional): List of items to write.
            options (:obj:`dict`, optional): options to pass to each item.

        """

        ### Store the version of STEAM used in all containers
        version_node = 'steam_version'
        if version_node not in self.toc:
            self.add( version_node, steam.__version__ )

        fmode = mode

        #! This is the existing path
        if path is None:
            fpath = self.hdf5_path
            if mode is None:
                fmode = self.hdf5_mode

        #! This is a new file
        else:
            fpath = path
            fmode = mode
            if os.path.isfile(fpath) and not mode is 'w':
                logger.warning("File will be overwritten by container write: {}".format(hdf_path))
                raise ValueError("Must explicitly state mode='w' if writing to an existing file.")

        if fmode == 'r':
            raise IOError("Cannot write to HDF5 in read-only mode")

        logger.info("Writing container to : {}".format(fpath))
        logger.info("Writing container with mode : {}".format(fmode))
        if fmode == 'w':
            logging.info(" - Overwriting any existing file.")

        hdf5  = self.open_store(fpath,mode=fmode)

        # Get the name of every object
        items_list = items
        if items_list is None:
            items_list = list(self.toc.keys())

        #! Check to see if I need to delete any of these (only if I'm appending to myself
        resize  = False
        if not mode == 'w':
            # List that will contain things that still need work after this
            new_delete = []
            for entry in self.to_delete:
                (dname,dpath) = entry
                
                # If I specified a list of items, then I only want to work on them
                # If this isn't one of those, then save it for later
                if items is not None:
                    if dname not in items:
                        # If not, then I need to keep it for next time
                        new_delete.append(entry)
                        continue
                resize = True

                logging.info("Deleting {} from container".format("/"+dname))
                #! Delete the thing
                #!  - Make the path HDF5-safe
                hdf5.remove("/"+steam.util.string_to_hdf5_key(dname))
                #! Delete the toc entry
                toc_path = "/table_of_contents/"+dpath
                logging.info("Deleting {} from container".format(toc_path))
                hdf5.remove(toc_path)

            self.to_delete = new_delete

        #! Run through and make sure that we have everything
        for obj_name in sorted(items_list):

            #! If this thing was deleted, then move on
            if obj_name not in self.obj:
                continue
 
            #! Get the object and toc information          
            try:
                obj           = self.obj[obj_name]
                (tpath,ttype) = self.toc[obj_name]
            except KeyError:
                print(" !! '{}' object not found in container!".format(obj_name))
                hdf5        = self.close_store(hdf5)
                return
            except:
                hdf5        = self.close_store(hdf5)
                raise

        #! Run through the toc and write each entry
        for obj_name in sorted(items_list):

            #! If this thing was deleted, then move on
            if obj_name not in self.obj:
                continue
 
            #! Get the object and toc information          
            try:
                obj           = self.obj[obj_name]
                (tpath,ttype) = self.toc[obj_name]
            except:
                hdf5        = self.close_store(hdf5)
                raise

            #! The type may have changed since it was first added, so figure it out
            #! and update the toc:
            var_type       = self.get_container_type(obj)
            self.toc[obj_name] = (tpath,var_type) 

            #! The path in the container is just the name plus "/"
            #!  - Make the path HDF5-safe
            root = "/"+steam.util.string_to_hdf5_key(obj_name)
            logger.info("Writing object '{}' to container as {}".format(obj_name,root))
                
        
            #! If it does not have a 'write_hdf5' method, then just pickle it
            if (var_type == "pickle"):
                steam.util.write_hdf5_pnode(hdf5,root,obj,compress=True)
            else:
                #! Check to see if we need to pass this options
                if (obj_name in options):
                    obj.write_hdf5(hdf5,
                                              root=root,
                                              options=options[obj_name]
                                              )
                else:
                    obj.write_hdf5(hdf5,
                                              root=root
                                              )

        self.write_toc(hdf5,items=items_list)

        hdf5        = self.close_store(hdf5)

        if resize:
            steam.container.resize_hdf5(fpath)

    def add(self,node,node_obj):
        """ Adds item to HDF5 file and TOC.

        The object will be added to the HDF5 file.  This is a simple
        call to each object looking for the .hdf5_write() command.
        If it doesn't have an hdf5_write() command, then it's stored
        as a pickled object.

        The location in the HDF5 file where the object will be stored is
        the node.  The node is parsed to get the resulting root_dir and
        variable name for the item.  HDF5 files store nodes similar to
        files on a file system.  The leading '/' for the HDF5 root is
        assumed and should not be provided.

        Ex:  database/first -> store item 'first' in location '/database/'
             The item will be stored at self.objs['database/first']

        The HD5F name cannot contain a '-', '.', or begin with a non-letter, non-_.
        Those are replaced in the path with "_d5s_", "_p5s_", and "_h5s_",
        respectively.

        Args:
            node (:obj:`string`): Where to store the object.
            node_obj: The object to be stored in the container.
        """

        import os.path
        import re

        if node is "table_of_contents":
            print("Cannot name item '{}'".format(node))
            return

        path = node

        # Make sure it begins with the root location
        path = "/"+path

        # Make sure noone put two slashes in it
        path = path.replace("//","/")

        # The variable name is the entire path minus the leading "/"
        var_name = path.replace("/","",1)

        # Make the path HDF5-safe
        path = steam.util.string_to_hdf5_key(path)
 
        #print("Adding ",var_name,"at",path)

        if (var_name == ""):
            print(" Node cannot end in a slash :",path,"\n")
            return

        #! Check to see if the thing exists already.  You can't put something
        #! inside of another object or replace one that is already there without
        #! deleting it first.
        if var_name in self.toc:
            raise KeyError(" Container already contains item at {}!".format(var_name))
        for item in self.obj:
            pattern = re.match("^{}/".format(item),var_name)
            if pattern is not None:
                logging.error(" - Existing object: {}".format(item))
                logging.error(" - New      object: {}".format(var_name))
                raise KeyError("Cannot place item inside of another item in container")

        #! Make sure it doesn't conflict with anything that's already there
        try:
            #exec("self.obj['{}']".format(var_name))
            self.obj[var_name]
        except KeyError:
            #! This means it wasn't defined.  Good.
            pass
        else:
            raise ValueError("Variable {} already exists in container!".format(var_name))

        logging.info("Adding {} to container".format(var_name))

        #! Set the variable inside this container
        self.obj[var_name] = node_obj
        #! Get a TOC location
        self.add_to_toc(var_name,node_obj)

        return

    def rm(self,obj_name):
        """ Remove a object from the container.

        Ags:
            obj_name (:obj:`string`): The name of the object to remove.
        """

        #! Make sure this exists
        try:
            self.toc[obj_name]
        except KeyError:
            raise NameError("Object '{}' not found in container!".format(obj_name))
         
        #! Add it to the list of things to delete
        (tpath,ttype) = self.toc[obj_name]
        self.to_delete.append((obj_name,tpath))

        #! Remove it
        logging.info("Removing {} from container".format(obj_name))
        try:
            # It might not exist if we haven't .read() the container
            del self.obj[obj_name]
        except:
            pass
        self.toc.pop(obj_name)

def resize_hdf5(fname):
    """ This will resize the HDF5 file to remove unused space.

    After removing objects, the space in the container is not cleared.
    This will attempt to resize the file to reclaim that space.

    The approved way from HDF5 docs is to make a new file with only
    the things you want to keep.
    """

    from os         import rename,stat

    logging.info("Attempting to resize HDF5 file, '{}'".format(fname))
    logging.info(" - Start size: {:.3f}M".format(stat(fname).st_size/1024/1024))
    # Temporary outupt file to use
    tmp = 'tmp_resize.h5'
    # Resize command
    cmd = "ptrepack -o --chunkshape=auto --propindexes {} {}".format(fname, tmp)
    #print(cmd)
    # Run the command, check for errors
    (stdout,stderr,code) = steam.util.run_command(cmd)

    # Check for code
    if code != 0:
        logging.error("Command failed: cmd")
        logging.error("STDOUT: {}".format(stout))
        logging.error("STDERR: {}".format(sterr))
        raise Exception("Error resizing HDF5 file, '{}'".format(fname))
    else:
        # Move temp file back to original location
        rename(tmp,fname)
    logging.info(" - End   size: {:.3f}M".format(stat(fname).st_size/1024/1024))
