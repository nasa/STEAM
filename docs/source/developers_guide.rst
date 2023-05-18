Developer's Guide
==================

.. toctree::
   :maxdepth: 2

   developers_guide


Introduction
------------------

This is a quick guide detailing some developer conventions.

I've Created/Changed a Thing
----------------------------

Great!  Let's make sure that you've documented the change, added a unit test so that it remains functional in the future, and included any relevant examples:

1. Check to make sure that docstrings are correct.

   * Follow the conventions here:
   * Check with ``make docs`` to ensure it was parsed correctly.

2. Create a unit test for your method/object in ``tests/``

   * New objects can be give their own directory.  All new directory names must end in ``_test``.
   * New test scripts can be added inside of any directory.  New scripts must start with ``test``.  There is an example script to start from: ``tests/test_template.py``.
   * New methods can be added to any script.  New methods must start with ``test_``.
   * Test it: ``./tests/run_all.sh``.

3. Update/Create relevant examples and documentation:

   * Update any relevant portion of the User's Manual: ``docs/source/users_guide.rst`` and ``docs/source/ug.*.rst``.
   * You can add/modify the examples in the manual.  The actual example code resides in ``examples/scripts/*.py`` and any necessary files for the examples are also in ``examples/``.
   * For examples, use ``.. literalinclude`` blocks.  There are a number of examples, but the important things to keep in mind are:

      * All examples should be functions in a module inside of ``examples/scripts/*.py``.  It should execute without raising an exception (this is tested in ``tests/users_man_tests``).
      * Examples' code blocks are surrounded by ``CODE_START`` and ``CODE_END``.
      * Examples' output blocks are surrounded by ``OUTPUT="""`` and ``""" END``.  This ensures that it conforms to python syntax, does not contain '#' and can be cleanly parsed.

4. Check that all tests work after your modifications.

   * Ensure that you have all of the changes from develop : 

      * ``git checkout develop``
      * ``git pull``
      * ``git checkout -``
      * ``git merge --no-ff develop``
   * Run all of the tests and fix any errors: ``./tests/run_all.sh``.

5. Commit and merge into develop: 

   * ``git checkout develop``
   * ``git merge --no-ff NAME/FEATURE``
   * ``git push``

