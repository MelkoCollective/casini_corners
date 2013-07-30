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
from scipy import linalg, spatial

import maple  # @UnresolvedImport
from calculate import Calculate

def getSignatureMatrix(matrix_out,matrix_it,shape):
    """
    Get the pattern of a matrix. Elements are cast via float.  User is
    responsible to remove bits subject to roundoff, if they exist.
    Ex: [20,80,50,
         30,20,40,
         40,80,20]
          
         gives 
          
         [0,1,2,
          3,0,4,
          4,1,0]
           
    :param matrix: The input matrix.
    :param shape: The shape of the matrix.
    """
     
    unique_elements = {}
    count = 0

     
    for idx_1d,el in enumerate(matrix_it):
        if el not in unique_elements.keys():
            unique_elements[float(el)] = count
            count+=1
        matrix_out[idx_1d/shape[1],idx_1d%shape[1]] = unique_elements[float(el)]
        

class LatticeTest(unittest.TestCase):


    def testLatticeGen(self):
        
        calc = Calculate()
        
        #Lattice sizes
        sizes = [0,1,2,3]
        sets = {0:set([(0,0)]),
                1:set([(0,0),(0,1),(1,0),(1,1)]),
                2:set([(0,0),(0,1),(1,0),(1,1),(2,0),(0,2),(1,2),(2,1),(2,2)]),
                3:set([(0,0),(0,1),(1,0),(1,1),(2,0),(0,2),(2,1),(1,2),(2,2),(3,0),(0,3),(3,1),(1,3),(3,2),(2,3),(3,3)])
                }
        
        for size in sizes:
            lattice = calc.square_lattice(size)
            eq_(sets[size],set(tuple(tuple(x) for x in lattice)))


class CorrelationsTest(unittest.TestCase):
    #TODO: Currently limited to a small case.  In future, when program optimized, expand this!
    
    maple_link = None
    # For L = 3...
    X = None
    P = None
    XP = None
    precision = None
    
    @classmethod
    def setUpClass(cls):
        cls.calc = Calculate()
        if cls.X is None and cls.P is None:
            # Startup Maple. #TODO: Replace with sympy when they fix bug (see src.main). This is ugly since relying on a hardcoded dir.
            cls.maple_link = maple.MapleLink("/opt/maple13/bin/maple -tu")
            cls.precision = 25
            
            lattice = cls.calc.square_lattice(5)
            cls.X,cls.P = cls.calc.correlations(lattice,cls.maple_link,cls.precision)
        else:
            # TODO: rremove this if this doesn't show up
            print "Bug: This function has been called twice"
    
    def setUp(self):
        self.calc = Calculate()
        
    def testSymmetric(self):
        eq_(True,self.X.is_symmetric())
        eq_(True,self.P.is_symmetric())
        
    def testReal(self):
        eq_(True, self.X == self.X.conjugate())
        eq_(True, self.P == self.P.conjugate())
        
    def testPrecision(self):
        """
        See if the precision increases with the variable.
        """
        polygon = self.calc.square_lattice(1)
        precisions = [20,30,40,50]
         
        # Find X,P for each precision
        corrs_collection = []
        for precision in precisions:
            X,P = self.calc.correlations(polygon, self.maple_link, precision)
            corrs_collection.append([X,P])
             
        for i in range(0, len(precisions)-1):
            # get X and P pairs for precisions[i] and precisions[i+1]
            pairX = sp.array([corrs_collection[i][0],corrs_collection[i+1][0]])
            pairP = sp.array([corrs_collection[i][1],corrs_collection[i+1][1]])
             
            # Check if increasing the precision variable actually does this.
            bound = 10**(-precisions[i])
             
            eq_(True, ( abs(pairX[0] - pairX[1]) < bound ).all())
            eq_(True, ( abs(pairP[0] - pairP[1]) < bound ).all())
        
    def testPrecomputed(self):
        """
        Test the precomputed correlations feature.
        """
        precision = 25

        lattice3 = self.calc.square_lattice(3)
        lattice4 = self.calc.square_lattice(4)
        
        X3,P3 = self.calc.correlations(lattice3, self.maple_link, precision, verbose=True)
        X4,P4 = self.calc.correlations(lattice4, self.maple_link, precision, verbose=True)
        
        # With precomputed.
        precomputed_correlations = {}
        X3_quick, P3_quick, precomputed_correlations = self.calc.correlations(lattice3, self.maple_link, precision, precomputed_correlations, verbose=True)
        X4_quick, P4_quick, precomputed_correlations = self.calc.correlations(lattice4, self.maple_link, precision, precomputed_correlations, verbose=True)

        eq_(True, X3 == X3_quick)
        eq_(True, P3 == P3_quick)
        eq_(True, X4 == X4_quick)
        eq_(True, P4 == P4_quick)
        
    def testMatrixSymmetries(self):
        """
        X and P should have the same pattern as the distance matrix defined by
        the lattice.
        """
        precision = 20
        polygon = self.calc.square_lattice(5)
        X,P = self.calc.correlations(polygon, self.maple_link, precision)
        
        # Round X and P down so we can see if elements are distinct or not.
        X = sympy.matrix2numpy(X)
        P = sympy.matrix2numpy(P)
        X = X.astype('float')
        P = P.astype('float')
        
        # Get the pattern of the distance matrix.
        D = spatial.distance.cdist(polygon,polygon)
        
        # The pattern of the distance matrix
        D_pat = sp.zeros(D.shape)
        getSignatureMatrix(D_pat,sp.nditer(D),D.shape)
        
        # Get the pattern of X and P.
        X_pat = sp.zeros(X.shape)
        P_pat = sp.zeros(P.shape)
        getSignatureMatrix(X_pat,sp.nditer(X),X.shape)
        getSignatureMatrix(P_pat,sp.nditer(P),P.shape)
        
        # Check if patterns match.
        eq_(False,(D_pat - X_pat).all())
        eq_(False,(D_pat - P_pat).all())
        
  
        
class Correlator_ijTest(unittest.TestCase):
    # See if we match values from maple here (use varying precisions?)
    maple_link = None
    calc = None
    
    @classmethod
    def setUpClass(cls):
        # Startup Maple. #TODO: Replace with sympy when they fix bug (see src.main). This is ugly since relying on a hardcoded dir.
        cls.maple_link = maple.MapleLink("/opt/maple13/bin/maple -tu")
        
    def setUp(self):
        self.calc = Calculate()
        
    def testPrecision(self):
        """
        See if precision increases with the variable "precision".
        
        Note: There is some redundancy in that this test acts in much the same 
        way as testPrecision in the CorrelationsTest TestCase, but acts on a 
        set of only a few [i,j].  This could be seen as a quicker alternative.
        """
        # TODO: Is this actually quicker?
        ijs = [[1,2],[3,5],[2,8]]

        precisions = [20,30,40,50]
        
        for [i,j] in ijs:
            # Find X,P for each precision
            phipi_collection = []
            for precision in precisions:
                phi,pi = self.calc.correlators_ij(i, j, self.maple_link, precision)
                phipi_collection.append([phi,pi])
                
            for i in range(0, len(precisions)-1):
                # get X and P pairs for precisions[i] and precisions[i+1]
                pairPhi = [phipi_collection[i][0],phipi_collection[i+1][0]]
                pairPi = [phipi_collection[i][1],phipi_collection[i+1][1]]
                
                # Check if increasing the precision variable actually does this.
                bound = 10**(-precisions[i])
                
                eq_(True, abs(pairPhi[0] - pairPhi[1]) < bound )
                eq_(True, abs(pairPi[0] - pairPi[1]) < bound ) 
        
            
    def testSymmetry(self):
        # Calculate for i,j and j,i. See if get same results.
        
        # Some random ijs.
        ijs = [[1,2],[3,5],[2,8]]
        
        precision = 20
        
        for ij in ijs:
            i = ij[0]
            j = ij[1]
            phi_corr, pi_corr = self.calc.correlators_ij(i, j, self.maple_link, precision)
            phi_corr2, pi_corr2 = self.calc.correlators_ij(j, i, self.maple_link, precision)
            eq_(phi_corr,phi_corr2)
            eq_(pi_corr,pi_corr2)
            
    def testValues(self):
        # Test values against maple values.
        maple_vals = [[0.03490480875099482936704802,], \
                      [0.01359884308448908943342273,], \
                      [0.009657881846921001497223915,]]
          
        # Some random ijs.
        ijs = [[1,2],[3,5],[2,8]]
          
        precision = 25
          
        for count, [i,j] in enumerate(ijs):
            phi_corr, pi_corr = self.calc.correlators_ij(i, j, self.maple_link, precision)
            eq_(phi_corr,maple_vals[count][0])
            eq_(pi_corr,maple_vals[count][1])
        
class CalculationsTest(unittest.TestCase):
    '''
    Test logical features that result from the calculations.
    '''
    
    def testPolygonScramble(self):
        # Here, rearrange polygon randomly and repeatedly get results. See if come out different.

        # Presets.
        precision = 20
        maple_link = maple.MapleLink("/opt/maple13/bin/maple -tu") #TODO: hopefully remove this maple dependency...too hardwirey.
        L = 3
        calc = Calculate()
        
        polygon = calc.square_lattice(L)
        
        X,P = calc.correlations(polygon, maple_link, precision)

        # Get eigenvals
        with sympy.mpmath.extraprec(500):
            XP = X*P
        XPnum = sympy.matrix2numpy(XP)
        sqrt_eigs = sp.sqrt(linalg.eigvals(XPnum))
        sqrt_eigs_orig = sqrt_eigs.real
        
        # Scramble polygon.
        sp.random.shuffle(polygon)
        
        # Perform the integration again.        
        X,P = calc.correlations(polygon, maple_link, precision)

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

    def testDiag2Chain(self):
        # Perform the diagonalization for a 2-site chain. Compare
        # an exact diagonalization against scipy.
        
        # First calculate X,P
        calc = Calculate()
        polygon = sp.array([[0,0],[2,0]]) 
        maple_link = maple.MapleLink("/opt/maple13/bin/maple -tu")
        precision = 25
        
        X,P = calc.correlations(polygon, maple_link, precision)
        
        entropy = calc.entropy(X, P, 1, precision,False)

        # Manual diag.
        # [ a, b,
        #   c, d]
        XP = X*P
        a = XP[0,0]
        b = XP[0,1]
        c = XP[1,0]
        d = XP[1,1]
        
        eig_plus = ((a+d)+sympy.mpmath.sqrt((a+d)**2-4*(a*d-b*c)))/2
        eig_minus = ((a+d)-sympy.mpmath.sqrt((a+d)**2-4*(a*d-b*c)))/2
        
        sqrt_eigs = [sympy.mpmath.sqrt(eig_plus),sympy.mpmath.sqrt(eig_minus)]
        
        S = 0
        for vk in sqrt_eigs:
            S += ((vk + 0.5)*sympy.mpmath.log(vk + 0.5) - (vk - 0.5)*sympy.mpmath.log(vk - 0.5))
            
        # Scipy operates in double, so we check equality up to a tolerance
        bound = 10**(-13)
        eq_(True, abs(S - entropy) < bound)
        
#     def testFit(self):
#         pass        
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()