'''
Created on May 31, 2013

@author: ghwatson
'''
from __future__ import division

from scipy import linalg, spatial
import scipy as sp
from sympy.mpmath import mpf, extraprec, cos, sqrt, log, pi, linspace, isinf, isnan
import sympy.mpmath
# TODO: The maple implementation is ugly, here.  Should pass namespace to the 
# maple object at initiation.
from maple import MapleLink, maple_EllipticK, maple_EllipticE, maple_EllipticF
import multiprocessing
import Queue

import re

class Calculate(object):
    '''
    Class containing calculations used in Casini's paper.
    '''
    
    _guard_start = 100 # Remembers the guard bits required to make the
                       # calculation work.
    
#     @staticmethod
    def correlations(self, polygon, maple_link, precision, correlators=None, verbose=False):
        
        sympy.mpmath.mp.dps = precision
        
        # Find all displacement vectors.
        dist_vects = [abs(p1-p2) for p1 in polygon for p2 in polygon]
        
        # The correlators depend symmetrically on i,j. To optimize calculations,
        # trim down to the effectively unique displacement vectors.
        dist_vects_ordered = [sp.array([p[1],p[0]]) if p[0] > p[1] else p for p in dist_vects]
        unique_ij = list(set(tuple(p) for p in dist_vects_ordered))
        
        # Allocate X and P.
        n = len(polygon)
        X = sp.zeros((n,n))
        P = sp.zeros((n,n))

        # Get all the unique correlators.
        for [i,j] in unique_ij:
            if correlators is not None and (i,j) in correlators.keys():
                phi_corr, pi_corr = correlators[(i,j)]
            else:
                phi_corr, pi_corr = self.correlators_ij(i,j,maple_link,precision)
                
                if verbose == True:
                    print "Calculated integrals for i,j = {0}:".format([i,j])
                    print "phi: {0}".format(phi_corr)
                    print "pi:  {0}".format(pi_corr)
                
            correlators[(i,j)] = [phi_corr, pi_corr]
            
        # Fill out X,P.
        for n1,p1 in enumerate(polygon):
            for n2,p2 in enumerate(polygon):
                [i,j] = abs(p1-p2)
                if i > j:
                    i,j = j,i
                X[n1,n2],P[n1,n2] = correlators[(i,j)]
            
        if correlators is not None:
            return X,P,correlators
        else:
            return X,P

    # A variant of the above function which will distribute the integrals
    # across many processes to utilize multiple cores. 
    def correlations_multicore(self, polygon, maple_dir, precision, correlators=None, verbose=False):
        
        sympy.mpmath.mp.dps = precision
        
        # Find all displacement vectors.
        dist_vects = [abs(p1-p2) for p1 in polygon for p2 in polygon]
        
        # The correlators depend symmetrically on i,j. To optimize calculations,
        # trim down to the effectively unique displacement vectors.
        dist_vects_ordered = [sp.array([p[1],p[0]]) if p[0] > p[1] else p for p in dist_vects]
        unique_ij = list(set(tuple(p) for p in dist_vects_ordered))
        
        # Allocate X and P.
        n = len(polygon)
        X = sp.zeros((n,n))
        P = sp.zeros((n,n))
        
        # Resources for allowing parallelization across multiple cores. 
        jobs = []
        q = multiprocessing.Queue()

        # Get all the unique correlators.
        for [i,j] in unique_ij:
            if correlators is not None and (i,j) in correlators.keys():
                phi_corr, pi_corr = correlators[(i,j)]

                correlators[(i,j)] = [phi_corr, pi_corr]
            else:
                def worker(q,i,j):
                    phi_corr, pi_corr = self.correlators_ij(i,j,maple_dir,precision)
                    
                    if verbose == True:
                        print "Calculated integrals for i,j = {0}:".format([i,j])
                        print "phi: {0}".format(phi_corr)
                        print "pi:  {0}".format(pi_corr)

                    q.put([(i,j),phi_corr, pi_corr])
                p = multiprocessing.Process(target=worker, args=(q,i,j,))
                p.start()
                jobs.append(p)

        # Complete any unfinished jobs.
        for j in jobs:
            j.join()
            
        # Get all the correlators from the finished processes.
        while True:
            try:
                (i,j),phi,pi = q.get(block = False)
                correlators[(i,j)] = [phi,pi]
            except Queue.Empty:
                break
            
        # Fill out X,P.
        for n1,p1 in enumerate(polygon):
            for n2,p2 in enumerate(polygon):
                [i,j] = abs(p1-p2)
                if i > j:
                    i,j = j,i
                X[n1,n2],P[n1,n2] = correlators[(i,j)]
            
        if correlators is not None:
            return X,P,correlators
        else:
            return X,P
        
#     @staticmethod
#     def correlators_ij(self, i, j, maple_link, precision):
    def correlators_ij(self, i, j, maple_dir, precision):
        '''
        Calculate the integral cos(ix)cos(jy)/sqrt(2(1-cos(x))+2(1-cos(y))) for x,y = -Pi..Pi.
        This is done by computing the inner integral symbolically, and the outer numerically.
        The integral is done on a quarter quadrant since it is symmetric.
        
        :param i: lattice index i
        :param j: lattice index j
        :param maple_link: the MapleLink class to communicate with Maple
        :param precision: the arithmetic precision
        '''
        # Set the mpmath precision.
        sympy.mpmath.mp.dps = precision
        
        maple_link = MapleLink(maple_dir,precision)
        
        # Quicker if the higher frequency occurs in numerically solved outer
        # integral.
        if i>j:
            i,j = j,i
        
        # Symbolically solve inner integral, using Maple:
#         phi_str = "cos({0}*x)/sqrt(2*(1-cos(x))+2*(1-cos(y)))".format(int(i))
#         phi_integ_str = "int({0},x=0..Pi) assuming y >= 0;".format(phi_str)
#         phi_integ_str = "simplify(int({0},x=0..Pi) assuming y >= 0);".format(phi_str)
#         pi_str = "cos({0}*x)*sqrt(2*(1-cos(x))+2*(1-cos(y)))".format(int(i))
#         pi_integ_str = "int({0},x=0..Pi) assuming y >= 0;".format(pi_str)
#         pi_integ_str = "simplify(int({0},x=0..Pi) assuming y >= 0);".format(pi_str)
#         inner_phi_str = maple_link.query(phi_integ_str)
#         inner_pi_str = maple_link.query(pi_integ_str)
#         
#         
#         def replacer(matchobj):
#             newstr = "mpf('{0}')".format(matchobj.group(0))
#             return newstr
#         inner_phi_str_old = re.sub(r"[0-9]*\.?[0-9]+",replacer,inner_phi_str)
#         inner_pi_str_old = re.sub(r"[0-9]*\.?[0-9]+",replacer,inner_pi_str)
#         
#         replace_dict = {'^':'**'}
#         for maple_entry,py_entry in replace_dict.iteritems():
#             inner_phi_str_old = inner_phi_str_old.replace(maple_entry,py_entry)
#             inner_pi_str_old = inner_pi_str_old.replace(maple_entry,py_entry)
 
        # Create function using maple output. #TODO: switch out the eval for a parser when possible. eval is dangerous.
#         def phi_inner_integral(y):
#             out = eval(inner_phi_str)
#             return out*cos(int(j)*y)
#         def pi_inner_integral(y):
#             out = eval(inner_pi_str)
#             return out*cos(int(j)*y)
        
        # New code here:
        phi_str = "cos({0}*x)/sqrt(2*(1-cos(x))+2*(1-cos(y)))".format(int(i))
#         phi_integ_str = "int({0},x=0..Pi) assuming y >= 0;".format(phi_str)
        phi_integ_str = "simplify(int({0},x=0..Pi) assuming y >= 0);".format(phi_str)
        pi_str = "cos({0}*x)*sqrt(2*(1-cos(x))+2*(1-cos(y)))".format(int(i))
#         pi_integ_str = "int({0},x=0..Pi) assuming y >= 0;".format(pi_str)
        pi_integ_str = "simplify(int({0},x=0..Pi) assuming y >= 0);".format(pi_str)

#         inner_phi_str = maple_link.query(phi_integ_str)
#         maple_link.parse(inner_phi_str)
#         phi_stack = maple_link.get_parsed_stack()

        inner_pi_str = maple_link.query(pi_integ_str)
        maple_link.parse(inner_pi_str)
        pi_stack = maple_link.get_parsed_stack()

        inner_phi_str = maple_link.query(phi_integ_str)
        maple_link.parse(inner_phi_str)
        phi_stack = maple_link.get_parsed_stack()
        
        def phi_inner_integral(y):
            vars = {'y':y}
            out = maple_link.eval_stack(phi_stack, vars)
            return out*cos(int(j)*y)
        def pi_inner_integral(y):
            vars = {'y':y}
            out = maple_link.eval_stack(pi_stack, vars)
            return out*cos(int(j)*y)
         
        # Perform the outer integrals.  
        # TODO: This blip of code is messy. There must be a more concise way 
        #       to achieve the implemented efficiency + readability here.
        def int_fail(integ):
            return (isnan(integ) or isinf(integ))
        for guard_bits in xrange(self._guard_start,1000,100):
            phi_integ = sympy.mpmath.quad(extraprec(guard_bits)(phi_inner_integral),[0,pi])
            if int_fail(phi_integ):
                continue
            pi_integ = sympy.mpmath.quad(extraprec(guard_bits)(pi_inner_integral),[0,pi])
            if int_fail(pi_integ):
                continue
            self._guard_start = guard_bits
            break
        else:
            raise ValueError("Hit guard bit tolerance of 1000!") #TODO: make this an input
        phi_integ *= mpf('1')/(2*pi**2)
        pi_integ *= mpf('1')/(2*pi**2)
        
        return phi_integ,pi_integ

    @staticmethod
    def entropy(X, P, n, precision, truncate, verbose=False):
        '''
        Using the correlator matrices, find the nth entropy.
        :param X: spatial correlator matrix (sympy)
        :param P: momentum correlator matrix (sympy)
        :param n: Renyi index
        :param precision: arithmetic precision 
        :param verbose: runtime output flag
        '''
#         with extraprec(500):
#             XP = X*P
        XP = sp.dot(X,P)
        
        # DEBUGGING: Saving the correlator matrices.
        prec = sympy.mpmath.mp.dps
        sp.savetxt("xmat.txt", X, fmt='%.{0}f'.format(prec),delimiter=",")
        sp.savetxt("pmat.txt",P, fmt='%.{0}f'.format(prec),delimiter=",")
                
        # Scipy eigenvalues
#         XPnum = sympy.matrix2numpy(XP)
#         sqrt_eigs = sp.sqrt(linalg.eigvals(XPnum.astype(sp.float32)))
        sqrt_eigs = sp.sqrt(linalg.eigvals(XP))
        
        # Check that the eigenvalues are well-defined.
        to_remove = []
        for i, eig in enumerate(sqrt_eigs):
            if eig.real <= 0.5:
                if truncate:
                    # Remove the eigenvalue
                    to_remove.append(i)
                    pass
                else:
                    raise ValueError("At least one of the eigenvalues of sqrt(XP) is below 0.5! \n eig = {0}".format(eig))
            if eig.imag != 0:
#                 raise ValueError("Warning: getting imaginary components in eigenvalues! \n imag = {0}".format(eig.imag))
                print "Warning: got an imaginary eigvalue component: " + str(eig.imag)
            
        sqrt_eigs = sp.delete(sqrt_eigs,to_remove)
       
        # Convert to float. Chop off imaginary component.
        sqrt_eigs = sqrt_eigs.real
        
        # Calculate entropy.
        S_n = 0
#         eps = 1.e-12
        eps = 0 # DEBUGGING
        if n == 1:
            for vk in sqrt_eigs:
#                 S_n += ((vk + 0.5)*log(vk + 0.5) - (vk - 0.5)*log(vk - 0.5 + eps))
                S_n += ((vk + 0.5)*log(vk + 0.5) - (vk - 0.5)*log(vk - 0.5))
        else:
            for vk in sqrt_eigs:
                S_n += log((vk + 0.5)**n - (vk - 0.5)**n)
            S_n *= 1./(n-1)
            
        if verbose==True:
            print "Calculated entropy of {0}".format(S_n)
            
        S_n = float(S_n)
        return S_n

    @staticmethod
    def square_lattice(L):
        '''
        Generates an array of points making up a square LxL lattice.
        
        :param L: the dimension of the square lattice.
        '''
        x = sp.linspace(0,L,L+1)
        coord_arrays = sp.meshgrid(x,x)
        polygon = (sp.dstack(coord_arrays))
        polygon = sp.reshape(polygon,((L+1)**2,2))    
        return polygon
