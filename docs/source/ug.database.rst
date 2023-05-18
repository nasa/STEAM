Database
----------
Database objects are derived from table objects, so manipulation of the data inside of a database uses identical methods as those used for tables.  Databases, however, include additional methods to handle grids and solutions, identifying dependant and independant variables, and interfacing with interpolation space objects.  Interpolation of scalar values and distributed maps (solutions) uses identical syntax and can be performed simultaneously.  Databases also provide methods for efficiently storing data, grids, and solutions into containers.

The simplest way to create a database is similar to creating a STEAM table file.  It can include as many commented lines at the top of the file (begining with a '#') followed by a header line and then the data.  The header line should not have a comment character and both it and the data must be white-space delimited.  Some of the columns in the file should correspond to independant paramers, others will be dependant parameters, and still more columns can be included for reference.  This is sufficient for a scalar database.

A 'solution database' is one that may or may not contain scalar data alongside the independant parameters, but also contains mesh/solution maps.  There are many possible formats for the mesh and solution data (see Mesh and Solution objects).  Each solution can have an individual mesh or they can share a common mesh.  If there is more than one mesh for the solution maps and the solution files do not contain their own mesh, then a mesh datafile must be provided (below).  If all solutions will share a common mesh and the solutions do not contain their mesh, then that mesh must be loaded prior to the solutions.  If all solutions will share a common mesh and the solution files do not contain their mesh, then the mesh must be loaded prior to loading the solutions (using the read_mesh() method).

To form this database, your input data file must contain the columns 'soln_path', 'soln_type', and (optionally) 'mesh_name':

* **soln_path:** The file path to the solution file.  This is relative to the input data file's directory or an absolute path.
* **soln_type:** The type of solution file.  See the Solution object for more details on possible types.
* **mesh_name:** The name of the mesh in the accomanying mesh file.  If there is only one mesh for all solution files, then this column can be excluded.  If the solution files contain their meshes (triq, for instance), then this column can be excluded.

If required, the mesh datafile must be provided and must have the columns: 'mesh_name', 'mesh_path', and 'mesh_type':

* **mesh_name:** A name for the mesh.  Used to identify the mesh in the input data file.
* **mesh_path:** The file path to the mesh file.  This can be relative to the current directory or an absolute path.
* **mesh_type:** The type of mesh file.  See the Mesh object for more details on possible types.

TODO: db.interpolate keyword 'cols' should be 'vars'?

Databases can be used as iterators.  When doing so, each case in the database is returned as a DBPoint inside of the loop context.  
CEbner: By reference?  Copy?  Cam you open/close each solution?  How do I update things (example below)?

Hashes
~~~~~~~~
CEbner: Describe, briefly what the hashes are; when and why we use them.  I think it's important to mention this to the user so 
that they can help diagnose their own issues.  If we tell them that the hash is the combination of these things, then 
they might be able to see where conflicts occure.

On Disk
~~~~~~~~~

For large databases, it might be desirable to store them on disk and not in memory.  This will reduce the speed at which the data can be read and written, but dramatically reduces the memory footprint.  To take advantage of this for a database. set the object's onDisk attribute to True.  All functionality available to a database in memory is also available to one that is on disk.

A database that is on disk must placed into a container and the container must have a location on disk available for use.  So the way to initialize a database as on disk is to do the following:


.. literalinclude:: ../../examples/scripts/database.py
   :caption: Create an database onDisk.
   :pyobject: onDiskExample
   :start-after: CODE_START
   :end-before: CODE_END
   :dedent: 4
   :linenos:

TODO: Make a method to have a database be onDisk and have it check to make sure that HDF5 file is set, exists, and can be opened.

Examples
~~~~~~~~~~~~
Files referenced in this section are collected in the `examples/` directory.  To not clutter 
the installation directory, copy the `examples/` directory and use it as your working directory 
for the following examples.

.. contents:: :local:

Read Scalar Database
^^^^^^^^^^^^^^^^^^^^^^^^

Reading the database file is straightforward:


.. literalinclude:: ../../examples/scripts/database.py
   :caption: Create an database onDisk.
   :pyobject: scalar_read
   :start-after: CODE_START
   :end-before: CODE_END
   :dedent: 4
   :linenos:


Once it is read, you can ensure that it has read the data file by looking at its data.  There are two ways.  
First, 'print(db)' yields a complete picture of it and its attributes - note that the DataFrame reports 
non-zero size.  The second way is to print the data directly:


.. literalinclude:: ../../examples/scripts/database.py
   :caption: Output
   :pyobject: scalar_read
   :start-after: OUTPUT
   :end-before: END
   :dedent: 4


For a database, the '.data' attribute is a Pandas DataFrame.  You can work with it directly to investigate
what variables exist, their ranges, or convert them to NumPy arrays:

.. literalinclude:: ../../examples/scripts/database.py
   :caption: Create an database onDisk.
   :pyobject: scalar_read
   :start-after: CODE_START2
   :end-before: CODE_END2
   :dedent: 4
   :linenos:


Perform Scalar Interpolation at Point
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Set Independant/Dependant Variables**
"""""""""""""""""""""""""""""""""""""""

Variables in the database are refered to by name and are case sensitive.  Independant varibles are used when creating interpolation spaces with which to do scalar/solution interpolation.  The dependant variabels are those variables that should be interpolated when performing interpolation.  If no dependant varibles are specified, then all non-independant variables are considered dependant.  Any entry in a 'SOLN' column is assumed to be a map of data and is always considered dependant.

To set the variables for this example, make 'MACH' and 'ALPHA' the independant parameters by providing a list of them.  They are stored as an attribute after performing a check to ensure that they exist in the data.  We will also set the dependant variables as "VAR1" and "VAR2".  This will ensure that any interpolated results do not contain interpolated 'RHO' and 'QBAR'.

Using 'print(db)' shows the state of the database after setting these variables.

.. literalinclude:: ../../examples/scripts/database.py
   :caption: Set database parameters
   :pyobject: scalar_set_vars
   :start-after: CODE_START
   :end-before: CODE_END
   :dedent: 4
   :linenos:

.. literalinclude:: ../../examples/scripts/database.py
   :caption: Output
   :pyobject: scalar_set_vars
   :start-after: OUTPUT
   :end-before: END
   :dedent: 4

**Create Interpolation Space**
"""""""""""""""""""""""""""""""""""""""

To perform interpolation, we first have to make an interpolation space.  These are detailed elsewhere.  For this simple database, we want to use all of our data in a single interpolation space.  Since we have specified the independant parameters, the default method call is appropriate.

.. literalinclude:: ../../examples/scripts/database.py
   :caption: Create interpolation space
   :pyobject: scalar_mk_space
   :start-after: CODE_START
   :end-before: CODE_END
   :dedent: 4
   :linenos:

.. literalinclude:: ../../examples/scripts/database.py
   :caption: Output
   :pyobject: scalar_mk_space
   :start-after: OUTPUT
   :end-before: END
   :dedent: 4

This informs us that it used scaling factors of 1.0 for both Mach and Alpha.  For some databases, you may want to scale the parameters to form a unit space or use custom scaling for each parameter.  You can also specify to perform interpolation in a modified space, "MACH^2" or "log(VEL)".  See the documentation for details on how to do this.

**Perform Interpolation**
"""""""""""""""""""""""""

Now that we have an interpolation space, we can interpolate inside of the 2-D space by specifying values of Mach and Alpha.  The interpolation method expects either a dataframe of points to interpolate to, a STEAM table, or a dictionary.  However input is supplied, it must contain values for all of the independant paramters.  The data is return as a list of Database Points (DBPoints objects).

For instance, let's interpolate to a Mach/Alpha point and provide an additional parameter, too.  Extra parameter just pass through the interpolation process and both them and the independent parameters are returned in the results.  This is to aid in workflows where you may have several paramters in your case list that are not necessary for the current interpolation step but might be necessary in subsequent operations.  By including it all in the input to the interpolation call, it remains availible in the output as a consolidated list.

.. literalinclude:: ../../examples/scripts/database.py
   :caption: Interpolate to point
   :pyobject: scalar_interpolate
   :start-after: CODE_START
   :end-before: CODE_END
   :dedent: 4
   :linenos:

.. literalinclude:: ../../examples/scripts/database.py
   :caption: Output
   :pyobject: scalar_interpolate
   :start-after: OUTPUT
   :end-before: END
   :dedent: 4

Notice that MACH, ALPHA, and FOO are all populated as are the values for the independant varibles, VAR1 and VAR2.  If we wanted to add RHO and QBAR to the output, we would have to either add them to the dependant paramter list (`db.set_dep(['RHO','QBAR','VAR1','VAR2'])`), clear the dependant paramter list (`db.set_dep([])` - since the default uses all parameters), or manually request them in the interpolate call:


.. literalinclude:: ../../examples/scripts/database.py
   :caption: Interpolate with additional variables
   :pyobject: scalar_interpolate
   :start-after: CODE2_START
   :end-before: CODE2_END
   :dedent: 4
   :linenos:

.. literalinclude:: ../../examples/scripts/database.py
   :caption: Output
   :pyobject: scalar_interpolate
   :start-after: OUTPUT2
   :end-before: END2
   :dedent: 4


For this simple problem, VAR1 is equal to QBAR+ALPHA and VAR2 is equal to QBAR-ALPHA.  That is a linear relationship
and we used linear interpolation, so we can confirm that the values are correct.

One important piece of diagnostics is the interpolation stencil that was used
to calculate the output data.  This can be shown by looking at the `.stencil`
attribute from the returned DBPoint.  It shows the target point and the
database points with their location in space and weights.  In this example,
two points were averaged (each have a weight of 0.5) to create data at the
requested target.

.. literalinclude:: ../../examples/scripts/database.py
   :caption: Viewing the interpolation stencil.
   :pyobject: scalar_interpolate
   :start-after: CODE3_START
   :end-before: CODE3_END
   :dedent: 4
   :linenos:

.. literalinclude:: ../../examples/scripts/database.py
   :caption: Output
   :pyobject: scalar_interpolate
   :start-after: OUTPUT3
   :end-before: END3
   :dedent: 4

Perform Solution Interpolation to Trajectory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Read Data**
"""""""""""""

In this example, we will be reading a solution database that includes a single mesh and multiple solutions that use that mesh.  The solutions all use the same mesh and do not contain it, so it is read prior to the solutions.  The grid is store in a triq file, and the solutions are stored in the cdat.



.. literalinclude:: ../../examples/scripts/database.py
   :caption: Create solution database
   :pyobject: solution_read
   :start-after: CODE_START
   :end-before: CODE_END
   :dedent: 4
   :linenos:

.. literalinclude:: ../../examples/scripts/database.py
   :caption: Output
   :pyobject: solution_read
   :start-after: OUTPUT
   :end-before: END
   :dedent: 4

This database has three independant parameters:  MACH, ALPHA, and BETA.  MACH has a range 0:5 and ALPHA and BETA have ranges -5:5.

**Create Interpolation Space and Interpolate**
""""""""""""""""""""""""""""""""""""""""""""""

This example uses the solution database to interpolate to a trajectory, represented by a table.  Each trajectory point will be populated with and scalar data in the database as well as interpolated solution maps.


.. literalinclude:: ../../examples/scripts/database.py
   :caption: Create solution database
   :pyobject: solution_mk_space
   :start-after: CODE_START
   :end-before: CODE_END
   :dedent: 4
   :linenos:

.. literalinclude:: ../../examples/scripts/database.py
   :caption: Output
   :pyobject: solution_mk_space
   :start-after: OUTPUT
   :end-before: END
   :dedent: 4

Each database point in the output list contains .data with the flightspace point and scalar variables, a .mesh with the database's common mesh, and a .soln that contains the solution. It also contains a .stencil string describing which points were used and the interpolation weights.  These solutions can be written out for later analysis.

.. literalinclude:: ../../examples/scripts/database.py
   :caption: Print resulting DBPoint
   :pyobject: solution_mk_space
   :start-after: CODE_START2
   :end-before: CODE_END2
   :dedent: 4
   :linenos:

.. literalinclude:: ../../examples/scripts/database.py
   :caption: Output
   :pyobject: solution_mk_space
   :start-after: OUTPUT2
   :end-before: END2
   :dedent: 4


The interpolate command returned a warning because the tenth point in the trajectory was outside of the database.  In this case, the interpolator does not extrapolate or hold-last, instead it returns blank data:


.. literalinclude:: ../../examples/scripts/database.py
   :caption: Print out of bounds data
   :pyobject: solution_mk_space
   :start-after: CODE_START3
   :end-before: CODE_END3
   :dedent: 4
   :linenos:

.. literalinclude:: ../../examples/scripts/database.py
   :caption: Output
   :pyobject: solution_mk_space
   :start-after: OUTPUT3
   :end-before: END3
   :dedent: 4

TODO: Have write_cdat include a comment or something that can encode the .data from the DBPoint.

**Write Results**
"""""""""""""""""

   for (i,point) in enumerate(output):
      steam.io.write_cdat(name="point_{}.cdat".format(i),soln=point.soln)

Replace or Modify a Mesh or Solution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Depending on the usage, a database can store data in memory or on disk.  When modifying or replacing meshes or solutions that are already in a database, you need to use to `update_` methods described below.  This is to enable to database to ensure that consistancy is maintained across the database and on disk (if relevant).  For read-only access, use the `get_` methods.  These methods return a copy of the data stored in the database.::

 db.get_mesh()
 db.get_soln()
 db.get_point()

To replace or modify the data in a database, you must call the relevant `.update_` method.  The solution or mesh object is replaced wholesale by the new
object pass to the method.  The database negotiates linking if the object already exists elsewhere in the database, maintains connectivity between solutions
and meshes, updates object hashes for fast comparisons, and maintains the HDF5 if this is an `onDisk` instance.

In a database, each soln is refered to by
a unique `solnid` and each mesh a unique `meshid`.  Multiple database points could use the same `solnid` if the solution data was identical at those points.
Frequently, multiple database points (or the entire database) will use the same `meshid`.  For this reason, there are two similar calls when updating
objects:
 * **update_meshid**  :  Replace the mesh at a given `meshid`.  This can affect multiple database points if multiple points point to the `meshid`.
 * **update_meshidx** :  Replace the mesh at a given database index.  This will only affect the specified database point.
 * **update_solnid**  :  Replace the solution at a given `solnid`.  This can affect multiple database points if multiple points point to the `solnid`.
 * **update_solnidx** :  Replace the solution at a given database index.  This will only affect the specified database point.

As an example, this replaces the solution at index 4 with one where the solution quantity, `q1`, is doubled:

.. literalinclude:: ../../examples/scripts/database.py
   :caption: Update solution at index 4
   :pyobject: modify_soln_mesh
   :start-after: CODE_START
   :end-before: CODE_END
   :dedent: 4
   :linenos:

If multiple database points used that solution, then all of them could be updated with a slightly different call:

.. literalinclude:: ../../examples/scripts/database.py
   :caption: Update all solutions using a given solnid
   :pyobject: modify_soln_mesh
   :start-after: CODE_START2
   :end-before: CODE_END2
   :dedent: 4
   :linenos:

To modify the mesh, similar calls can be used:

.. literalinclude:: ../../examples/scripts/database.py
   :caption: Update a database's mesh
   :pyobject: modify_soln_mesh
   :start-after: CODE_START3
   :end-before: CODE_END3
   :dedent: 4
   :linenos:

**NOTE:**  When updating a solution, the mesh used by that solution is also updated to use the mesh that the solution uses at the time of the update call.

**NOTE2:**  When modifying or replacing a mesh, the solutions are not interpolated to the new mesh.  It is assumed that the mesh has been only cosmetically altered (moved, rotated, mirrored) such that the connectivity between the mesh and solution is the same.  In the above example, each solution in the database that was associated with 'g0001' will remain associated with it and the object referencing will be updated.

When closing a mesh, the code ensures that the number of vertices and elements does not change.  If it did, then it would invalidate the solutions' mapping to it.  This is only a cosmetic check, it is possible for the user to replace one mesh with another, entirely different mesh that had the same dimension and 'trick' STEAM.

If you intend to transform solutions from one mesh to another by more than a cosmetic change, use the ``db.transform_mesh(new_mesh)`` method.  See arguments and usage in the API documentation.

Iterate through a Database
^^^^^^^^^^^^^^^^^^^^^^^^^^

CEbner: How can iterate through a database and update the solutions?  Maybe renaming all of the variables from "q1" to "cp"?
