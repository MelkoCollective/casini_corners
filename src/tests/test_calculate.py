'''
Created on Jun 17, 2013

@author: ghwatson
'''
import unittest
from nose.tools import eq_  # @UnresolvedImport

import maple  # @UnresolvedImport
from calculate import Calculate
import sympy
import sympy.mpmath
import scipy as sp
from scipy import linalg
import random



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
    
    
        
class Correlator_ijTest(unittest.TestCase):
    # See if we match values from maple here (use varying precisions?)
    # Test symmetry of i,j.
    # Test entropy calculation somehow....(use matlab)
    pass

class PrecisionTest(unittest.TestCase):
    # Try varying precisions, see if the precision increases accordingly.
    #TODO: Just add this in to each of the other test cases.
    pass
        
class CalculationsTest(unittest.TestCase):
    
    def polygonScrambleTest(self):
        # Here, rearrange polygon randomly and repeatedly get results. See if come out different.
        # Generate a poly.
        precision = 20
        maple_link = maple.MapleLink("/Library/Frameworks/Maple.framework/Versions/12/bin/maple -tu") #TODO: hopefully remove this maple dependency...too hardwirey.
        
        polygon = Calculate.square_lattice(5)
        # Get eigenvals.
        X,P = Calculate.correlations(polygon, maple_link, precision)

        # Get eigenvals
        with sympy.mpmath.extraprec(500):
            XP = X*P
    #         v = XP.eigenvals()
    #         sqrt_eigs = []
    #         for eig, mult in v.iteritems():
    #             sqrt_eigs.append([sqrt(sympy.re(sympy.N(eig,precision)))] * mult)
    #         sqrt_eigs = list(chain.from_iterable((sqrt_eigs)))
        
        #Scipy eigenvalues.
        #TODO: see if this works, given mpmath structure.
        XPnum = sympy.matrix2numpy(XP)
        speigs = linalg.eigvals(XPnum)
        sqrtspeigs = sp.sqrt(speigs)
        sqrt_eigs = sqrtspeigs
        
        # Check that the eigenvalues are well-defined.
        for eig in sqrt_eigs:
            if eig.real <= 0.5:
                raise ValueError("At least one of the eigenvalues of sqrt(XP) is below 0.5! \n eig = {0}".format(eig))
            if eig.imag != 0:
                raise ValueError("Warning: getting imaginary components in eigenvalues! \n imag = {0}".format(eig.imag))
        
        # Convert to float. Chop off imaginary component.
        sqrt_eigs_orig = sqrt_eigs.real
        
        
        # Scramble the indices and rebuild the poly.
        sp.random.shuffle(polygon)
        
        # Perform the integration again.        
        X,P = Calculate.correlations(polygon, maple_link, precision)

        # Get eigenvals
        with sympy.mpmath.extraprec(500):
            XP = X*P
    #         v = XP.eigenvals()
    #         sqrt_eigs = []
    #         for eig, mult in v.iteritems():
    #             sqrt_eigs.append([sqrt(sympy.re(sympy.N(eig,precision)))] * mult)
    #         sqrt_eigs = list(chain.from_iterable((sqrt_eigs)))
        
        #Scipy eigenvalues.
        #TODO: see if this works, given mpmath structure.
        XPnum = sympy.matrix2numpy(XP)
        speigs = linalg.eigvals(XPnum)
        sqrtspeigs = sp.sqrt(speigs)
        sqrt_eigs = sqrtspeigs
        
        # Check that the eigenvalues are well-defined.
        for eig in sqrt_eigs:
            if eig.real <= 0.5:
                raise ValueError("At least one of the eigenvalues of sqrt(XP) is below 0.5! \n eig = {0}".format(eig))
            if eig.imag != 0:
                raise ValueError("Warning: getting imaginary components in eigenvalues! \n imag = {0}".format(eig.imag))
        
        # Convert to float. Chop off imaginary component.
        sqrt_eigs_scrambled = sqrt_eigs.real
        a = set(tuple(sqrt_eigs_orig))
        b = set(tuple(sqrt_eigs_scrambled))
        eq_(a,b)

    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()