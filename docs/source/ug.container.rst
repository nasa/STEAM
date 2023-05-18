Containers
------------

STEAM is designed to save objects into containers.  A container is an HDF5 file that includes a table of contents as well as all objects themselves.  It is designed to be cross-platform to facilitate distributed development.  Each container is flexible and allows for the user to store data in a directory-like structure.

STEAM-specific objects have specialized read/write routines to efficiently store their data in the HDF5 file.  General python objects are stored using pickle.  You can read and write to the same HDF5 file, but there is not much consistency checking to ensure that the HDF5 file is cleaned of old information when overwriting.

HDF5 files are a "Hierarchial Data Format" and contain groups and nodes that are similar to a file system directory tree.  A node can exist nested into any number of groups.  The user can take advantage of this in order to organize the contets of a container to suite their needs.  When storing objects to a container, the complete path to the node is required and the object will be referenced by that path.  The path to the node would properly contain an initial backslash, imply the root of the HDF5, but this is left off for convenience when working with containers in STEAM.

Containers in STEAM have a dual existance.  Any open container and its contents exist in memory.  Only during a ``.write()`` operation are any modifications to the container written to disk.  Modifications include adding, modifying, moving, or removing objects.  For new container objects, an HDF5 file is not created until the first ``.write()`` is called.

Database objects can be used in an ``onDisk`` state.  When ``onDisk``, the database object is aware of its location in the HDF5 file and can read/write to that location independant of the container.  Provided the user does not manipulate an ``onDisk`` database while also manipulating that same node in the container, there should not be conflicts.

Here is a quick example that shows the creation a new container.  Several STEAM objects are include as is a Python string object.  To access the stored data, a special dictionary in the container (``.obj{}``) will include a reference to the object.  The key to each objects is the path that was specified with first adding the object.

.. literalinclude:: ../../examples/scripts/container.py
   :caption: Creating a new container
   :pyobject: create
   :start-after: CODE_START
   :end-before: CODE_END
   :dedent: 4
   :linenos:

.. literalinclude:: ../../examples/scripts/container.py
   :caption: Output
   :pyobject: create
   :start-after: OUTPUT
   :end-before: END
   :dedent: 4

In this example, the ``mesh`` object was added twice.  While the container is in memory, both objects store a reference to the original object.  Once it is written to disk, however, the object is written twice into two distinct locations.  This essentially copies the object when the container is next loaded from disk into memory.

Context management for containers automatically call the ``.write()`` method when exiting the context.  The above example could be simplified by replaceing container creating with the following:

.. literalinclude:: ../../examples/scripts/container.py
   :caption: Creating a new container with context manager
   :pyobject: create_context
   :start-after: CODE_START
   :end-before: CODE_END
   :dedent: 4
   :linenos:

When obtaining a new ``Container()`` object, a path must be used.  If this is a path to an existing file, then the file will be opened as a container in read-only mode.  If the path does not exist, then the path will be stored and used to write the container when ``.write()`` is next called.  An existing container can be open in read/write mode if the ``mode`` keyword argument is used.

Extension
~~~~~~~~~~~~

The container class can add or write any class that has the ``.write_hdf5(hdf5,root="/")`` and ``.read_hdf5(hdf5,root="/")`` methods.  Where ``hdf5`` is an HDFStore object and ``root`` is the path inside the HDF5 to use as the base directory or Group.  As new classes are added, appropriate methods should be included to maintain compatibility.  You can also add dictionaries for options that can be passed to the read/write methods in order to allow flexibility to the user.  Use the existing methods from other objects as an example.  

Future work can include support for ``.to_fs()`` and ``.from_fs()`` methods that will convert classes into file system directories and files that can be manipulated or re-encoded.  This might provide an easy way to import/export from the proprietary container format.

Examples
~~~~~~~~~~~
Files referenced in this section are collected in the `examples/` directory.  To not clutter 
the installation directory, copy the `examples/` directory and use it as your working directory 
for the following examples.

.. contents:: :local:

Reading a Container
^^^^^^^^^^^^^^^^^^^

.. literalinclude:: ../../examples/scripts/container.py
   :caption: Reading a container
   :pyobject: read_container
   :start-after: CODE_START
   :end-before: CODE_END
   :dedent: 4
   :linenos:

.. literalinclude:: ../../examples/scripts/container.py
   :caption: Output
   :pyobject: read_container
   :start-after: OUTPUT
   :end-before: END
   :dedent: 4

Viewing the Contents of a Container
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: ../../examples/scripts/container.py
   :caption: Viewing a Container
   :pyobject: view_contents
   :start-after: CODE_START
   :end-before: CODE_END
   :dedent: 4
   :linenos:

.. literalinclude:: ../../examples/scripts/container.py
   :caption: Output
   :pyobject: view_contents
   :start-after: OUTPUT
   :end-before: END
   :dedent: 4

Removing Objects from a Container
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: ../../examples/scripts/container.py
   :caption: Creating a new container
   :pyobject: remove_items
   :start-after: CODE_START
   :end-before: CODE_END
   :dedent: 4
   :linenos:

.. literalinclude:: ../../examples/scripts/container.py
   :caption: Output
   :pyobject: remove_items
   :start-after: OUTPUT
   :end-before: END
   :dedent: 4
