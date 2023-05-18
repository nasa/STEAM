Data Tables
------------------

Data tables (or 'tables') are the underlying structure used to store database values and tabulated input to STEAM.  At their core, they are two-dimensional Pandas DataFrames.  The DataFrames have been wrapped by the STEAM table class in order to make manipulation more convenient, but Pandas operations will work on the table directly.

Below will follow some examples of usage based on the sample files in '/examples/table/'.  The files all have a common extension, '.cdat'.  In STEAM, the extension '.cdat' is used to denote a column-based ASCII data file.  It can include as many commented lines at the top of the file (beginning with a '#') followed by a header line and then the data.  The header line should not have a comment character and both it and the data must be white-space delimited.

TODO: Expand available comment characters, allow for non-white-space deliniation, allow for header to be commented for true compatibility with cdat library.

Examples
~~~~~~~~~~~

Files referenced in this section are collected in the `examples/` directory.  To not clutter 
the installation directory, copy the `examples/` directory and use it as your working directory 
for the following examples.

.. contents:: :local:

Load Table into STEAM
^^^^^^^^^^^^^^^^^^^^^^

Loading a table into steam is as easy as pointing a new object to the data file.

.. literalinclude:: ../../examples/scripts/table.py
   :caption: Creating a new Table
   :pyobject: load
   :start-after: CODE_START
   :end-before: CODE_END
   :dedent: 4
   :linenos:

.. literalinclude:: ../../examples/scripts/table.py
   :caption: Output
   :pyobject: load
   :start-after: OUTPUT
   :end-before: END
   :dedent: 4


Work with DataFrame Directly
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A table's 'data' attribute is a Pandas DataFrame and can be called directly:

.. literalinclude:: ../../examples/scripts/table.py
   :caption: Use Pandas DataFrame calls
   :pyobject: min_max
   :start-after: CODE_START
   :end-before: CODE_END
   :dedent: 4
   :linenos:

.. literalinclude:: ../../examples/scripts/table.py
   :caption: Outout
   :pyobject: min_max
   :start-after: OUTPUT
   :end-before: END
   :dedent: 4

Add/Remove Columns from Table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Some common operations have been turned into method of the table object in order to make the syntax more convenient.  To remove a variable from a table:


.. literalinclude:: ../../examples/scripts/table.py
   :caption: Remove Table Columns
   :pyobject: rm_column
   :start-after: CODE_START
   :end-before: CODE_END
   :dedent: 4
   :linenos:

.. literalinclude:: ../../examples/scripts/table.py
   :caption: Outout
   :pyobject: rm_column
   :start-after: OUTPUT
   :end-before: END
   :dedent: 4


To add a new column, you must have a two-dimensional array that has the same number of rows as your data.  Let's add two columns of data, one filled with 1.0 and the other with 2.0.


.. literalinclude:: ../../examples/scripts/table.py
   :caption: Add Column to Table
   :pyobject: add_column
   :start-after: CODE_START
   :end-before: CODE_END
   :dedent: 4
   :linenos:

.. literalinclude:: ../../examples/scripts/table.py
   :caption: Outout
   :pyobject: add_column
   :start-after: OUTPUT
   :end-before: END
   :dedent: 4


Write Table to ASCII File
^^^^^^^^^^^^^^^^^^^^^^^^^^
Write Table to HDF5 Container
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
