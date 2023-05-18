--[[
Module file for loading STEAM and all other necessary modules.  
The intention of this file to dramatically enhance portability
between STEAM users, at least while running on the FSL.
--]]

-- Load prerequisite modules
depends_on( "gcc/4.8.5" )
depends_on( "openmpi/3.1.4" )
depends_on( "libmesh/1.4.1" )
depends_on( "mpi4py/3.0.1-python3.10" )
depends_on( "anaconda/3-2022.05" )

--[[
Users will have their own installs of STEAM in different locations.
We need to find the Python executable directories based on the location of
this file and our knowledge of the STEAM directory structure.
--]]
this_file = myFileName()
this_dir, fname  = splitFileName( this_file )
-- Top STEAM dir is 3 levels up
relative_path = "utils/module_set/steam"
steam_dir = string.gsub( this_dir, relative_path, "" )
prepend_path	('PATH'      , steam_dir)
prepend_path	('PYTHONPATH', steam_dir)

--[[
In order to use conda, it must be activated with a script that is included
with it.
--]]
execute{ cmd=". /software/x86_64/anaconda/3-2022.05/etc/profile.d/conda.sh", modeA={"load"} }
