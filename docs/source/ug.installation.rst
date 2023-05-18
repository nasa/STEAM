Installing STEAM
~~~~~~~~~~~~~~~~~
STEAM comes with a module file for use in JSC's Flight Sciences Lab (FSL) and
other systems that use 
`Lmod <https://lmod.readthedocs.io/en/latest/index.html>`_.  This module file 
is located in ``/path/to/steam/utils/module_set/steam``.  In order to use
the module file, the ``module_set`` directory must be added to the user's
module path with this command::

    module use /path/to/steam/utils/module_set

This line can be copied into the user's ``.bashrc`` or ``.cshrc`` file in
order to make STEAM available in all shells.

Once the STEAM module has been added to the module path, you can check its 
availability and load it::

    user$ module use /path/to/steam/utils/module_set
    user$ module avalilable steam

    ---------------- /path/to/steam/utils/module_set -----------------
       steam/default


    user$ module purge
    user$ module list
    No modules loaded

    user$ module load steam
    user$ module list

    Currently Loaded Modules:
      1) gcc/4.8.5       3) libmesh/1.4.1            5) anaconda/3-2018-12
      2) openmpi/3.1.4   4) mpi4py/3.0.1-python3.7   6) steam/default

The included module automatically loads the other modules necessary for using
STEAM.  It also runs the ``conda.sh`` script that enables activation of the
requisite Python environment.

In order to simplify the installation of the provided STEAM module, a simple
script is included in the module directory, ``add_module_path.py``.  When
run without any arguments, this script will output the correct module command
for the user to be able to load the STEAM module.::

    cebner$ ./add_module_path.py
    module use /aerolab/cebner/cebner/steam/utils/module_set

    cebner$ module use /aerolab/cebner/cebner/steam/utils/module_set
    cebner$ module available steam

    ----------------/aerolab/cebner/cebner/steam/utils/module_set -----------------
       steam/default


Additionally, the user can provide an optional argument of their .bashrc or 
.cshrc (or analogous file for whatever shell they use) and have the script 
append the ``module use`` command to that file.  This will not automatically
load the STEAM module, but it will make it available in all shells in the 
future so that the user may call ``module load steam`` from anywhere at any 
time.

Managing Python Environments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
STEAM needs to run in a specific Python 3.7.4 environment that contains various 
modules.  In order to standardize this environment, we use ``conda``, a Python
package manager.  Assuming successful execution of the module set up 
described above, ``conda`` should be in your path any time you load the 
provided ``steam`` module.

The first time that you use STEAM, you'll need to create the ``conda`` 
environment that you'll use every time you use STEAM.  This is done using the
``py37.yml`` file that comes with the STEAM repository with the following 
command::

    conda env update -f /path/to/steam/utils/py37.yml --prune

One common bug while building conda environments involves conda's inability to
install the necessary packages for the environment.  This generally happens
because the user does not have write privileges to the location to which 
conda is attempting to save the Python packages.  The workaround is to 
explicitly specify where a user wants their packages installed.

This can be done very simply by creating a file in the user's home directory
called ``.condarc`` and pasting the following text in it::

    pkgs_dirs:
      - /home/user/.conda/pkgs

Conda will then create a directory called ``.conda`` in the user's home
directory and all the packages required by conda will be downloaded and
stored here in the future.

Once the ``conda env update`` command above has been run successfully,
a new environment called ``pysteam37`` should have been created.  To activate 
the environment and run code inside it::

    conda activate pysteam37

Though the environment only needs to be created once, this ``conda activate`` 
command will need to be issued in every terminal that will use STEAM before 
you run the code.
