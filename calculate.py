'''
Created on May 31, 2013

@author: ghwatson
'''
from __future__ import division

from scipy import linalg
import scipy as sp
from sympy.mpmath import mpf, extraprec, cos, log, pi, isinf, isnan
import sympy.mpmath
from maple import MapleLink
import multiprocessing
import Queue

import bresenham

# class Calculate(object):
#     '''
#     Class containing calculations used in Casini's paper.
#     '''
    
GUARD_START = 100 # Remembers the guard bits required to make the
                       # calculation work.
    
#     @staticmethod
def correlations(polygon, maple_dir, precision, correlators=None, verbose=False):
    
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
            phi_corr, pi_corr = correlators_ij(i,j,maple_dir,precision)
            
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

# Used in the below function as a worker in multiprocessing code. This is the 
# kernel that gets parallelized.
def worker(data):
    # Unpack data.
    i,j,verbose,maple_dir,precision = data
    
    phi_corr, pi_corr = correlators_ij(i,j,maple_dir,precision)
    
    if verbose == True:
        print "Calculated integrals for i,j = {0}:".format([i,j])
        print "phi: {0}".format(phi_corr)
        print "pi:  {0}".format(pi_corr)

    return [i,j,phi_corr,pi_corr]
# A variant of the above correlations function which will distribute the integrals
# across many processes to utilize multiple cores. 
def correlations_multicore(polygon, maple_dir, precision, correlators=None, verbose=False):
    
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
    
    # Record which correlators need to be computed.
    pool_inputs = []
    for [i,j] in unique_ij:
        if correlators is not None and (i,j) in correlators.keys():
            phi_corr, pi_corr = correlators[(i,j)]

            correlators[(i,j)] = [phi_corr, pi_corr]
        else:
            # Aggregate jobs for the pool.
            pool_inputs.append([i,j,verbose,maple_dir,precision])

    # Calculate correlators using a pool of workers.
    pool_size = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=pool_size)
    pool_outputs = pool.map(worker,pool_inputs)
    pool.close()
    pool.join()
    for i,j,phi,pi in pool_outputs:
        correlators[(i,j)] = [phi,pi]
        
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
    
def correlators_ij(i, j, maple_dir, precision):
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
    
    global GUARD_START
    
    maple_link = MapleLink(maple_dir,precision)
    
    # Quicker if the higher frequency occurs in numerically solved outer
    # integral.
    if i>j:
        i,j = j,i
    
    phi_str = "cos({0}*x)/sqrt(2*(1-cos(x))+2*(1-cos(y)))".format(int(i))
    phi_integ_str = "simplify(int({0},x=0..Pi) assuming y >= 0);".format(phi_str)
    pi_str = "cos({0}*x)*sqrt(2*(1-cos(x))+2*(1-cos(y)))".format(int(i))
    pi_integ_str = "simplify(int({0},x=0..Pi) assuming y >= 0);".format(pi_str)

    inner_pi_str = maple_link.query(pi_integ_str)
    maple_link.parse(inner_pi_str)
    pi_stack = maple_link.get_parsed_stack()

    inner_phi_str = maple_link.query(phi_integ_str)
    maple_link.parse(inner_phi_str)
    phi_stack = maple_link.get_parsed_stack()
    
    def phi_inner_integral(y):
        vardict = {'y':y}
        out = maple_link.eval_stack(phi_stack, vardict)
        return out*cos(int(j)*y)
    def pi_inner_integral(y):
        vardict = {'y':y}
        out = maple_link.eval_stack(pi_stack, vardict)
        return out*cos(int(j)*y)
     
    # Perform the outer integrals.  
    # TODO: This blip of code is messy. There must be a more concise way 
    #       to achieve the implemented efficiency + readability here.
    def int_fail(integ):
        return (isnan(integ) or isinf(integ))
    for guard_bits in xrange(GUARD_START,1000,100):
        phi_integ = sympy.mpmath.quad(extraprec(guard_bits)(phi_inner_integral),[0,pi])
        if int_fail(phi_integ):
            continue
        pi_integ = sympy.mpmath.quad(extraprec(guard_bits)(pi_inner_integral),[0,pi])
        if int_fail(pi_integ):
            continue
        GUARD_START = guard_bits
        break
    else:
        raise ValueError("Hit guard bit tolerance of 1000!") #TODO: make this an input
    phi_integ *= mpf('1')/(2*pi**2)
    pi_integ *= mpf('1')/(2*pi**2)
    
    return phi_integ,pi_integ

# @staticmethod
def entropy(X, P, n, precision, truncate, verbose=False):
    '''
    Using the correlator matrices, find the nth entropy.
    :param X: spatial correlator matrix (sympy)
    :param P: momentum correlator matrix (sympy)
    :param n: Renyi index
    :param precision: arithmetic precision 
    :param verbose: runtime output flag
    '''
    XP = sp.dot(X,P)
    
    # DEBUGGING: Saving the correlator matrices.
#         prec = sympy.mpmath.mp.dps
#         sp.savetxt("xmat.txt", X, fmt='%.{0}f'.format(prec),delimiter=",")
#         sp.savetxt("pmat.txt",P, fmt='%.{0}f'.format(prec),delimiter=",")
            
    # Scipy eigenvalues
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
    if n == 1:
        for vk in sqrt_eigs:
            S_n += ((vk + 0.5)*log(vk + 0.5) - (vk - 0.5)*log(vk - 0.5))
    else:
        for vk in sqrt_eigs:
            S_n += log((vk + 0.5)**n - (vk - 0.5)**n)
        S_n *= 1./(n-1)
        
    if verbose==True:
        print "Calculated entropy of {0}".format(S_n)
        
    S_n = float(S_n)
    return S_n

# @staticmethod
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

# @staticmethod
def circle_lattice(L):
    coord_tuples = bresenham.raster_disk(0,0,L)
    unique_tuples = set(coord_tuples)

    coords = [sp.array(i) for i in unique_tuples]
    return coords
    