Database Point
---------------

A database point is an object separate from a database.  It has two main uses:

1) To store the data returned during interpolation inside of a database.

2) To extract a case from a database or add a new case to a database.

Database points contains all variables, the mesh, and solution from a single index in a database.  When returned from
an interpolat cal, it also includes the interpolation stencil.  These are all contained in attributes of the DBPoint object::

 dbp = steam.database.DBPoint()
 dbp.data
 dbp.mesh
 dbp.soln
 dbp.stencil  # If from database interpolation


CEbner: This needs updated to conform the the new DBPoint usage.  Now, things are more object oriented and should use
the .{get,set}_{meshid,solnid,index,stencil} methods.  Also, points are copies of what was in the data, changing them does NOT update the
database.  If you use the .update_{mesh,soln} methods then it will replace what was in the database.  I think we should
have examples showing the disconnect and how to sync changes in a point to the source database.

Examples
~~~~~~~~~~~~~~~~~~~~
Files referenced in this section are collected in the `examples/` directory.  To not clutter 
the installation directory, copy the `examples/` directory and use it as your working directory 
for the following examples.

.. contents:: :local:

Extract a Point from a Database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For example, we can extract the first index from a solution database and observe that it has the mesh, soln, and freestream points associated with that index.


.. literalinclude:: ../../examples/scripts/dbpoint.py
   :caption: Extract point from database
   :pyobject: extract_point
   :start-after: CODE_START
   :end-before: CODE_END
   :dedent: 4
   :linenos:

.. literalinclude:: ../../examples/scripts/dbpoint.py
   :caption: Output
   :pyobject: extract_point
   :start-after: OUTPUT
   :end-before: END
   :dedent: 4

TODO: Have write_cdat include a comment or something that can encode the .data from the DBPoint.


Add a Point to a Database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let's expand the example database in the previous by holding-last in ALPHA.
To do that, we could take this point, modify its freestream parameters (in `.data`) and
then add it back to the database.  Continuing from the example above:

.. literalinclude:: ../../examples/scripts/dbpoint.py
   :caption: Add single DBPoint to database
   :pyobject: expand_alpha
   :start-after: CODE_START
   :end-before: CODE_END
   :dedent: 4
   :linenos:

.. literalinclude:: ../../examples/scripts/dbpoint.py
   :caption: Output
   :pyobject: expand_alpha
   :start-after: OUTPUT
   :end-before: END
   :dedent: 4

Or, to copy all of the cases with ALPHA set to -5 down to ALPHA -180, the following could be done:

.. literalinclude:: ../../examples/scripts/dbpoint.py
   :caption: Extend all ALPHA == -5 solutions in database
   :pyobject: expand_alpha
   :start-after: CODE_START2
   :end-before: CODE_END2
   :dedent: 4
   :linenos:

In this example, we duplicated several of the points in the database and their accompanying solution data.  
To be efficient, the STEAM database does not duplicate the solution data.  Instead, it merely tracks a 
second reference to the solution.  This is done by comparing the solution hash to existing hashes in the 
database.  You can confirm that by printing the database and observing that there are still only 58 
solutions but now 62 points in the DataFrame (`print(db)`).

.. literalinclude:: ../../examples/scripts/dbpoint.py
   :caption: Output
   :pyobject: expand_alpha
   :start-after: OUTPUT2
   :end-before: END2
   :dedent: 4

Modify a Point in a Database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For database points that are part of a database, changes to the mesh or solution in the database point
can be passed back to the originating database.  This performs the corresponding `.update_XXXXidx` method
calls from the database.

This example will start with a sample database containing a single grid.  Each point in the database will be updated to use a unique grid that has been translated in the x-direction.  The solutions will not be changed.  The mesh for case 0 remains the same because it is translated by 0 grid units.  Solutions could be updated with a similar call to `.update_soln()`.

.. literalinclude:: ../../examples/scripts/dbpoint.py
   :caption: Update the mesh for database points
   :pyobject: update_mesh
   :start-after: CODE_START
   :end-before: CODE_END
   :dedent: 4
   :linenos:

.. literalinclude:: ../../examples/scripts/dbpoint.py
   :caption: Relevant Output
   :pyobject: update_mesh
   :start-after: OUTPUT
   :end-before: END
   :dedent: 4
