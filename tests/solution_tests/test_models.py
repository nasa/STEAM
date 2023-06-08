#! /bin/env python

import unittest as ut
import steam

class TestProject( ut.TestCase ):

    def test_scalar_to_comp( self ):
        """ Create dummy data, operate on it, and check.
        """
        mesh = steam.mesh.mk_cube( 1, 5 )
        soln = steam.solution.uniform_soln( mesh, node_elem_flag = 'NODE',
                                            values = [1.0,10.0,20.0])
        soln.node_to_element()

        steam.models.scalar_to_component(soln,1,'q1',1.0)
        steam.models.scalar_to_component(soln,2,'q1',2.0)
        steam.models.scalar_to_component(soln,3,'q1',3.0)
        steam.models.scalar_to_component(soln,4,'q1',4.0,'repl')
        steam.models.scalar_to_component(soln,5,'q1',5.0,'repl')
        steam.models.scalar_to_component(soln,6,'q1',6.0,'repl')

        steam.models.scalar_to_component(soln,1,'q2',1.0,"+")
        steam.models.scalar_to_component(soln,2,'q2',2.0,"add")
        steam.models.scalar_to_component(soln,3,'q2',3.0,"+")
        steam.models.scalar_to_component(soln,4,'q2',4.0,"*")
        steam.models.scalar_to_component(soln,5,'q2',5.0,"mult")
        steam.models.scalar_to_component(soln,6,'q2',6.0,"*")

        steam.models.scalar_to_component(soln,1,'q3',1.0,"-")
        steam.models.scalar_to_component(soln,2,'q3',2.0,"sub")
        steam.models.scalar_to_component(soln,3,'q3',3.0,"-")
        steam.models.scalar_to_component(soln,4,'q3',4.0,"/")
        steam.models.scalar_to_component(soln,5,'q3',5.0,"div")
        steam.models.scalar_to_component(soln,6,'q3',6.0,"/")

        for comp in [1,2,3,4,5,6]:
            for ele in mesh.get_comp(comp):
                #! q1 is always component number
                self.assertEqual(
                    soln.data['q1'][ele],
                    comp
                    )

                #! q2 is 10 + comp if comp <= 3 and 10 * comp for comp >3
                if (comp < 4):
                    self.assertEqual(
                        soln.data['q2'][ele],
                        10.0 + comp
                        )

                if (comp > 3):
                    self.assertEqual(
                        soln.data['q2'][ele],
                        10.0 * comp
                        )

                #! q3 is 20 - comp if comp <= 3 and 20 / comp for comp >3
                if (comp < 4):
                    self.assertEqual(
                        soln.data['q3'][ele],
                        20.0 - comp
                        )

                if (comp > 3):
                    self.assertEqual(
                        soln.data['q3'][ele],
                        20.0 / comp
                        )


#! Run the tests
if __name__ == '__main__':
    ut.main()

