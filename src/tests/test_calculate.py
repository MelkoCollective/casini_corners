'''
Created on Jun 17, 2013

@author: ghwatson
'''
import unittest
from nose.tools import eq_  # @UnresolvedImport

import calculate, maple  # @UnresolvedImport


class LatticeTest(unittest.TestCase):


    def testLatticeGen(self):
        
        #Lattice sizes
        sizes = [0,1,2,3]
        sets = {0:set([(0,0)]),
                1:set([(0,0),(0,1),(1,0),(1,1)]),
                2:set([(0,0),(0,1),(1,0),(1,1),(2,0),(0,2),(1,2),(2,1),(2,2)]),
                3:set([(0,0),(0,1),(1,0),(1,1),(2,0),(0,2),(2,1),(1,2),(2,2),(3,0),(0,3),(3,1),(1,3),(3,2),(2,3),(3,3)])
                }
        
        for size in sizes:
            lattice = calculate.generate_square_lattice(size)
            eq_(sets[size],set(tuple(tuple(x) for x in lattice)))


class CorrelationsTest(unittest.TestCase):
    #TODO: Currently limited to a small case.  In future, when program optimized, expand this!
    
    # The correlation matrices
    X = None
    P = None
    # ...and the product of them
    XP = None
    
    @classmethod
    def setUpClass(cls):
        
        # Startup Maple. #TODO: Replace with sympy when they fix bug (see src.main). This is ugly since relying on a hardcoded dir.
        maplelink = maple.MapleLink("/Library/Frameworks/Maple.framework/Versions/12/bin/maple -tu")
                
        # We just do for a particular case, for now, since slow.
        lattice = calculate.generate_square_lattice(3)
        
        cls.X,cls.P = calculate.calculate_correlations(lattice,maplelink,25, True)
        cls.XP = cls.X*cls.P
    
    def testSymmetric(self):
        eq_(True,self.X.is_symmetric())
        eq_(True,self.P.is_symmetric())
        
    def testReal(self):
        eq_(True, self.X == self.X.conjugate())
        eq_(True, self.P == self.P.conjugate())

    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()