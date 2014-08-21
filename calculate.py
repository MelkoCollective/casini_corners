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

GUARD_START = 100 # Remembers the guard bits required to make the
                       # calculation work.
    
def correlations(polygon, maple_dir, precision, correlators=None, verbose=False):
    '''
    Given a point cloud defined by polygon, this function will find the X and P
    matrices of the phi and pi field correlations of free bosons on a square
    2D lattice.
    
    :param polygon: a list of points on the lattice ex: [[1,2],[5,8],[3,4]]
    :param maple_dir: directory to the commandline maple
    :param precision: precision with which to perform the integration
    :param correlators: precomputed correlators to avoid recomputation
    :param verbose: the user can turn standard output on and off.
    '''
    
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

def correlations_multicore(polygon, maple_dir, precision, correlators=None, verbose=False):
    '''
    A variant of the above correlations function which will distribute the 
    integrals across many subprocesses to utilize multiple cores.
    '''
    
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
    :param maple_dir: the directory to the commandline maple.
    :param precision: the arithmetic precision (for numerical integration)
    '''
    # Set the mpmath precision.
    sympy.mpmath.mp.dps = precision
    
    global GUARD_START
    
    maple_link = MapleLink(maple_dir,precision)
    
    # Quicker if the higher frequency occurs in numerically solved outer
    # integral.
    if i>j:
        i,j = j,i
        
    # Perform the inner integrals.
    
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
    
    def int_fail(integ):
        return (isnan(integ) or isinf(integ))
    # This loop will increase the bits of precision in the integral calculations
    # if they don't converge (low enough precision results in INF as a result).
    # The global is used to remember the precision that works.
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
        raise ValueError("Hit guard bit tolerance of 1000!") #TODO: make this an input if needed
    # Prefactors
    phi_integ *= mpf('1')/(2*pi**2)
    pi_integ *= mpf('1')/(2*pi**2)
    
    return phi_integ,pi_integ

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
            else:
                raise ValueError("At least one of the eigenvalues of sqrt(XP) is below 0.5! \n eig = {0}".format(eig))
        if eig.imag != 0:
            print "Warning: got an imaginary eigvalue component: " + str(eig.imag)
        
    sqrt_eigs = sp.delete(sqrt_eigs,to_remove)
   
    # Chop off imaginary component.
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

def circle_lattice(L):
    '''
    Generates an array of points making up a rasterized disc with integer
    radius L.
    
    :param L: the radius of the disk
    '''
    coord_tuples = bresenham.raster_disk(0,0,L)
    unique_tuples = set(coord_tuples)

    coords = [sp.array(i) for i in unique_tuples]
    return coords
    