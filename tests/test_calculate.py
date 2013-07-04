'''
Created on Jun 17, 2013

@author: ghwatson
'''
import unittest
import random

from nose.tools import eq_  # @UnresolvedImport
import sympy
import sympy.mpmath
import scipy as sp
from scipy import linalg

import maple  # @UnresolvedImport
from calculate import Calculate


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
            lattice = Calculate.square_lattice(size)
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
        if cls.X is None and cls.P is None:
            # Startup Maple. #TODO: Replace with sympy when they fix bug (see src.main). This is ugly since relying on a hardcoded dir.
            maplelink = maple.MapleLink("/Library/Frameworks/Maple.framework/Versions/12/bin/maple -tu")
                    
            # We just do for a particular case, for now, since slow.
            lattice = Calculate.square_lattice(3)
            
            cls.X,cls.P = Calculate.correlations(lattice,maplelink,20)
            cls.XP = cls.X*cls.P
    
    def testSymmetric(self):
        eq_(True,self.X.is_symmetric())
        eq_(True,self.P.is_symmetric())
        
    def testReal(self):
        eq_(True, self.X == self.X.conjugate())
        eq_(True, self.P == self.P.conjugate())
        
    # Test the precomputed ability.
    def testPrecomputed(self):
        pass
    
        
class Correlator_ijTest(unittest.TestCase):
    # See if we match values from maple here (use varying precisions?)
    # Test symmetry of i,j.
    # Test entropy calculation somehow....(use matlab)
    def testSymmetry(self):
        # Calculate for i,j and j,i. See if get same results.
        
        # Some random ijs.
        ijs = [[1,2],[2,3],[3,5],[1,5],[2,8]]
        
        maple_link = maple.MapleLink("/Library/Frameworks/Maple.framework/Versions/12/bin/maple -tu")
        precision = 20
        
        for ij in ijs:
            i = ij[0]
            j = ij[1]
            phi_corr, pi_corr = Calculate.correlators_ij(i, j, maple_link, precision)
            phi_corr2, pi_corr2 = Calculate.correlators_ij(j, i, maple_link, precision)
            eq_(phi_corr,phi_corr2)
            eq_(pi_corr,pi_corr2)
            
    def testValues(self):
#         # Test values against maple values.
#         maple_vals = []
#         
#         # Some random ijs.
#         ijs = [[1,2],[2,3],[3,5],[1,5],[2,8], [10,10]]
#         
#         maple_link = maple.MapleLink("/Library/Frameworks/Maple.framework/Versions/12/bin/maple -tu")
#         precision = 30
#         
#         for count, ij in enumerate(ijs):
#             i = ij[0]
#             j = ij[1]
#             phi_corr, pi_corr = Calculate.correlators_ij(i, j, maple_link, precision)
#             eq_(phi_corr,maple_vals[count][0])
#             eq_(pi_corr,maple_vals[count][])
        pass

class PrecisionTest(unittest.TestCase):
    # Try varying precisions, see if the precision increases accordingly.
    #TODO: Just add this in to each of the other test cases.
    # NOTE: this is too slow a thing to currently implement. Maybe in future
    # when calculations are quicker.
    pass
        
class CalculationsTest(unittest.TestCase):
    '''
    Test logical features that result from the calculations.
    '''
    
    def testPolygonScramble(self):
        # Here, rearrange polygon randomly and repeatedly get results. See if come out different.

        # Presets.
        precision = 20
        maple_link = maple.MapleLink("/Library/Frameworks/Maple.framework/Versions/12/bin/maple -tu") #TODO: hopefully remove this maple dependency...too hardwirey.
        L = 3
        
        
        polygon = Calculate.square_lattice(L)
        
        X,P = Calculate.correlations(polygon, maple_link, precision)

        # Get eigenvals
        with sympy.mpmath.extraprec(500):
            XP = X*P
        XPnum = sympy.matrix2numpy(XP)
        sqrt_eigs = sp.sqrt(linalg.eigvals(XPnum))
        sqrt_eigs_orig = sqrt_eigs.real
        
        # Scramble polygon.
        sp.random.shuffle(polygon)
        
        
        
        # Perform the integration again.        
        X,P = Calculate.correlations(polygon, maple_link, precision)

        # Get eigenvals
        with sympy.mpmath.extraprec(500):
            XP = X*P
        XPnum = sympy.matrix2numpy(XP)
        sqrt_eigs = sp.sqrt(linalg.eigvals(XPnum))
        sqrt_eigs_scrambled = sqrt_eigs.real
        
        
        
        # Compare.
        a = set(tuple([round(float(x),10) for x in sqrt_eigs_orig]))
        b = set(tuple([round(float(x),10) for x in sqrt_eigs_scrambled]))
        eq_(a,b)

    def testValues(self):    
        # Test a diagonalization against matlab values.
        pass
    
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()