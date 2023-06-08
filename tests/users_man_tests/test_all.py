#! /bin/env python

import unittest as ut
import steam

class TestUserExamples( ut.TestCase ):

    @classmethod
    def setUpClass( cls ):
        """Create temporary files. This method gets called once, at the very
        beginning of testing.
        """

        # Copy the contents of "examples" to the current directory
        import os
        test_home = os.environ['STEAM_TEST_HOME']
        example_home = test_home+"/../examples/"
        import shutil
        for dir in os.listdir(example_home):
            if (not os.path.exists(dir)):
                shutil.copytree("{}/{}".format(example_home,dir),dir)

        # Copy all of the users guide test modules
        example_home = test_home+"/../examples/scripts"
        if (not os.path.isdir('modules')):
            shutil.copytree(example_home,'modules')

    def test_all_examples( self ):
        """Run all methods for all test modules in
        """
        import importlib
        import sys
        sys.path.insert(0,'modules')
        import glob
        module_files = glob.glob('modules/*.py')
        modules = [mod.replace("modules/","").replace(".py","") for mod in module_files]
        for module in modules:
            mod = importlib.import_module(module)
            #import code
            #code.interact(local=locals())
            func_list = [func for func in dir(mod) if callable(getattr(mod, func))]
            for func in func_list:
                print("Testing {:20s}: {}".format(module,func))
                call = getattr(mod,func)
                # Store the output somewhere else
                with open('{}.{}.out'.format(module,func), 'w') as f:
                    sys.stdout = f
                    call()  # Check for exception
                # return output to screen
                sys.stdout = sys.__stdout__

    @classmethod
    def tearDownClass( cls ):
        """This method gets called once, at the very end of testing in order
        to clean out temporary files.
        """

#! Run the tests
if __name__ == '__main__':

    ut.main()

