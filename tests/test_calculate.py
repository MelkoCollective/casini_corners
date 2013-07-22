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
    
    maple_link = None
    # For L = 3...
    X = None
    P = None
    XP = None
    precision = None
    
    @classmethod
    def setUpClass(cls):
        if cls.X is None and cls.P is None:
            # Startup Maple. #TODO: Replace with sympy when they fix bug (see src.main). This is ugly since relying on a hardcoded dir.
            cls.maple_link = maple.MapleLink("/opt/maple13/bin/maple -tu")
            cls.precision = 25
                    
            
            lattice = Calculate.square_lattice(5)            
            cls.X,cls.P = Calculate.correlations(lattice,cls.maple_link,cls.precision)
            cls.XP = cls.X*cls.P
    
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
        polygon = Calculate.square_lattice(3)
        precisions = [20,30,40,50]
        tolerance = 2
        
        # Find X,P for each precision
        corrs_collection = []
        for precision in precisions:
            X,P = Calculate.correlations(polygon, self.maple_link, precision)
            corrs_collection.append([X,P])
            
        for i in range(0, len(precisions)-1):
            # get X and P pairs for precisions[i] and precisions[i+1]
            pairX = sp.array([corrs_collection[i][0],corrs_collection[i+1][0]])
            pairP = sp.array([corrs_collection[i][1],corrs_collection[i+1][1]])
            
            # Check if increasing the precision variable actually does this.
            digits_diff = precisions[i] - precisions[i+1]
            bound = 10^(digits_diff + tolerance)
            eq_(True, ( abs(pairX[0] - pairX[1]) < bound ).all())
            eq_(True, ( abs(pairP[0] - pairP[1]) < bound ).all())
        
    def testPrecomputed(self):
        """
        Test the precomputed correlations feature.
        """
        precision = 25
        lattice3 = Calculate.square_lattice(3)
        lattice4 = Calculate.square_lattice(4)
        
        X3,P3 = Calculate.correlations(lattice3, self.maple_link, precision, verbose=True)
        X4,P4 = Calculate.correlations(lattice4, self.maple_link, precision, verbose=True)
        
        # With precomputed.
        precomputed_correlations = {}
        X3_quick, P3_quick, precomputed_correlations = Calculate.correlations(lattice3, self.maple_link, precision, precomputed_correlations, verbose=True)
        X4_quick, P4_quick, precomputed_correlations = Calculate.correlations(lattice4, self.maple_link, precision, precomputed_correlations, verbose=True)

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
        polygon = Calculate.square_lattice(5)
        X,P = Calculate.correlations(polygon, self.maple_link, precision)
        
        # Round X and P down so we can see if elements are distinct or not.
        X = sympy.matrix2numpy(X)
        P = sympy.matrix2numpy(P)
        X = X.astype('float')
        P = P.astype('float')
        
        print X
        print ' '
        print P

        # Get the pattern of the distance matrix.
        D = spatial.distance.cdist(polygon,polygon)
        
        # The pattern of the distance matrix
        D_pat = sp.unique(D,return_inverse=True)[1]
        
        # Get the pattern of X and P.
        X_pat = sp.unique(X,return_inverse=True)[1]
        P_pat = sp.unique(P,return_inverse=True)[1]
        
        # Check if patterns match.
        print D_pat
        print ' '
        print X_pat
        eq_(False,(D_pat - X_pat).all())
        eq_(False,(D_pat - P_pat).all())
        
    @classmethod
    def getSignatureMatrix(matrix,shape):
        """
        Get the pattern of a matrix.
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
        
        for idx_1d,el in enumerate(matrix):
            if el not in unique_elements.keys():
                unique_elements[el] = count
                count+=1
            matrix[idx_1d/shape[1],idx_1d%shape[0]] = unique_elements[el]
        
        
class Correlator_ijTest(unittest.TestCase):
    # See if we match values from maple here (use varying precisions?)
    maple_link = None
    
    @classmethod
    def setUpClass(cls):
        # Startup Maple. #TODO: Replace with sympy when they fix bug (see src.main). This is ugly since relying on a hardcoded dir.
        cls.maple_link = maple.MapleLink("/opt/maple13/bin/maple -tu")
        
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
        tolerance = 2
        
        for [i,j] in ijs:
            # Find X,P for each precision
            phipi_collection = []
            for precision in precisions:
                phi,pi = Calculate.correlators_ij(i, j, self.maple_link, precision)
                phipi_collection.append([phi,pi])
                
            for i in range(0, precisions.size()-1):
                # get X and P pairs for precisions[i] and precisions[i+1]
                pairPhi = [phipi_collection[i][0],phipi_collection[i+1][0]]
                pairPi = [phipi_collection[i][1],phipi_collection[i+1][1]]
                
                # Check if increasing the precision variable actually does this.
                digits_diff = precision[i] - precision[i+1]
                bound = 10^(digits_diff + tolerance)
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
            phi_corr, pi_corr = Calculate.correlators_ij(i, j, self.maple_link, precision)
            phi_corr2, pi_corr2 = Calculate.correlators_ij(j, i, self.maple_link, precision)
            eq_(phi_corr,phi_corr2)
            eq_(pi_corr,pi_corr2)
            
    def testValues(self):
        # Test values against maple values.
#         maple_vals = [[0.03490480875099482936704802,], \
#                       [0.01359884308448908943342273,], \
#                       [0.009657881846921001497223915,], \
#                       [,]]
         
#         # Some random ijs.
#         ijs = [[1,2],[3,5],[2,8],[10,10]]
#          
#         precision = 25
#          
#         for count, [i,j] in enumerate(ijs):
#             phi_corr, pi_corr = Calculate.correlators_ij(i, j, self.maple_link, precision)
#             eq_(phi_corr,maple_vals[count][0])
#             eq_(pi_corr,maple_vals[count][1])
        pass
        
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
        
        # Have some matlab input here in comments to show what can be done.
        
        # Input the same matrix into sympy with mpf.
        
        # Diagonalize using scipy (with and without doubln down)
        
        # Check the eigenvalues against each other.
        pass
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()