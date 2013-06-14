'''
Created on May 31, 2013

@author: ghwatson
'''
from __future__ import division

from scipy import integrate, linalg, spatial
import scipy as sp

import sympy
from sympy.mpmath import cos, sqrt, log, pi
import sympy.mpmath

from tool.box import maple_EllipticK, maple_EllipticE, maple_EllipticF, flatten

def calculate_entropy(polygon, n, maple_link):
    '''
    Calculates the nth Renyi entropy for the polygonal set, as Casini does
    numerically in his paper.
    
    :param polygon: a Nx2 scipy.array of points composing the polygon region.
    :param n: the Renyi index.
    :param maple_link: link to Maple's command line, using MapleLink class.
    '''
    # Compute distance matrix of all possible pairs.
    D = spatial.distance.cdist(polygon,polygon)
    
    # Get the unique distances making up D, as well as id's of D giving these
    # and inverse id's to rebuild D from the unique distances.
    [unique,idx_1d,inv_idx] = sp.unique(D,True,True)
    
    # Convert 1d indices to 2d indices
    idx_2d = sp.unravel_index(idx_1d,D.shape)
    idx_2d = sp.dstack((idx_2d[0],idx_2d[1]))[0]
    
    # Unflatten inv_idx
    inv_idx = sp.reshape(inv_idx,D.shape)
      
    # Calculate all distinct correlator values for V.
    unique_phi_correlations = sympy.mpmath.zeros(unique.shape[0],1) 
    unique_pi_correlations = sympy.mpmath.zeros(unique.shape[0],1)
    
    for idx_1d, [r,r_prime] in enumerate(idx_2d):
        
        # Take difference in coordinates of the pairs to get i,j (i.e. set origin
        # at r for calculating this correlator).
        [i,j] = polygon[r_prime] - polygon[r]
                        
        # Symbolically solve inner integral, using Maple:
        phi_str = "cos({0}*x)/sqrt(2*(1-cos(x))+2*(1-cos(y)))".format(i)
        phi_integ_str = "int({0},x=0..Pi) assuming y >= 0;".format(phi_str)
        pi_str = "cos({0}*x)*sqrt(2*(1-cos(x))+2*(1-cos(y)))".format(i)
        pi_integ_str = "int({0},x=0..Pi) assuming y >= 0;".format(pi_str)
        inner_phi_str = maple_link.query(phi_integ_str)
        inner_pi_str = maple_link.query(pi_integ_str)

        # Create function using maple output. #TODO: switch out the eval for a parser when possible. eval is dangerous.
        def phi_inner_integral(y):
            out = eval(inner_phi_str)
            return out*cos(j*y)
        def pi_inner_integral(y):
            out = eval(inner_pi_str)
            return out*cos(j*y)
            
        # Perform the outer integrals.
        sympy.mpmath.mp.dps = 35
        phi_integ = sympy.mpmath.quad(phi_inner_integral,[0,pi])
        pi_integ = sympy.mpmath.quad(pi_inner_integral,[0,pi])
                
        # Save.
        unique_pi_correlations[idx_1d] = pi_integ*(1./(2*pi**2))
        unique_phi_correlations[idx_1d] = phi_integ*(1./(2*pi**2))
        
    # Populate matrix elements. 
    X = sympy.zeros(inv_idx.shape[0],inv_idx.shape[1])
    P = sympy.zeros(inv_idx.shape[0],inv_idx.shape[1])
    for i in xrange(inv_idx.shape[0]):
        for j in xrange(inv_idx.shape[1]):
            X[i,j] = unique_phi_correlations[inv_idx[i,j]]
            P[i,j] = unique_pi_correlations[inv_idx[i,j]]
            
    # Get the eigenvalues of sqrt(XP)
    XP = X*P
    v = XP.eigenvals()
    sqrt_eigs = []
    for eig, mult in v.iteritems():
        sqrt_eigs.append([sqrt(sympy.re(sympy.N(eig,20)))] * mult)
    sqrt_eigs = list(flatten(sqrt_eigs))
        
    # Calculate entropy.
    S_n = 0
    if n == 1:
        for vk in sqrt_eigs:
            S_n += ((vk + 0.5)*log(vk + 0.5) - (vk - 0.5)*log(vk - 0.5))
    else:
        for vk in sqrt_eigs:
            S_n += log((vk + 0.5)**n - (vk - 0.5)**n)
        S_n *= 1./(n-1)
        
    return S_n

def generate_square_lattice(L):
    '''
    Generates an array of points making up a square LxL lattice.
    
    :param L: the dimension of the square lattice.
    '''
    x = sp.linspace(0,L,L+1)
    coord_arrays = sp.meshgrid(x,x)
    polygon = (sp.dstack(coord_arrays))
    polygon = sp.reshape(polygon,((L+1)**2,2))    
    return polygon


# Old Scipy version. Didn't work. ---------------------------------------------



def calculate_entropy_old(polygon,n, maple_link):
    '''
    Calculates the nth Renyi entropy for the polygonal set, as Casini does
    numerically in his paper.
     
    :param polygon: a Nx2 scipy.array of points composing the polygon region.
    :param n: the Renyi index.
    '''
    #TODO: For small clusters, does it matter to even do this, or easier to just recompute the integrals? Run some speed tests, if this takes off.
    # Compute distance matrix of all possible pairs.
    D = spatial.distance.cdist(polygon,polygon)
     
    # Get the unique distances making up D, as well as id's of D giving these
    # and inverse id's to rebuild D from the unique distances.
    [unique,idx_1d,inv_idx] = sp.unique(D,True,True)
     
    # Convert 1d indices to 2d indices
    idx_2d = sp.unravel_index(idx_1d,D.shape)
    idx_2d = sp.dstack((idx_2d[0],idx_2d[1]))[0]
     
    # Unflatten inv_idx
    inv_idx = sp.reshape(inv_idx,D.shape)
         
    # Calculate all distinct correlator values for V.
    unique_phi_correlations = sp.zeros(unique.shape)
    unique_pi_correlations = sp.zeros(unique.shape)
    for idx_1d, [r,r_prime] in enumerate(idx_2d):
         
        # Take difference in coordinates of the pairs to get i,j (i.e. set origin
        # at r_prime for calculating this correlator).
        [i,j] = polygon[r_prime] - polygon[r]
                 
        
        # Symbolically solve inner integral, using Maple:
        phi_str = "cos({0}*x)/sqrt(2*(1-cos(x))+2*(1-cos(y)))".format(i)
        phi_integ_str = "int({0},x=0..Pi) assuming y >= 0;".format(phi_str)
        pi_str = "cos({0}*x)*sqrt(2*(1-cos(x))+2*(1-cos(y)))".format(i)
        pi_integ_str = "int({0},x=0..Pi) assuming y >= 0;".format(pi_str)
        inner_phi_str = maple_link.query(phi_integ_str)
        inner_pi_str = maple_link.query(pi_integ_str)
 
        # Create Python function from MapleLink
        def phi_inner_integral(y):
            out = eval(inner_phi_str)
            return out
        def pi_inner_integral(y):
            out = eval(inner_pi_str)
            return out
             
        # Perform outer integral:
        phi_integ = integrate.quad(phi_inner_integral,0,pi,epsabs=1.49e-13,epsrel=1.49e-13,weight='cos',wvar=j)
        pi_integ = integrate.quad(pi_inner_integral,0,pi,epsabs=1.49e-13,epsrel=1.49e-13,weight='cos',wvar=j)
                 
        unique_phi_correlations[idx_1d] = phi_integ[0]*(1./(2*pi**2))
        unique_pi_correlations[idx_1d] = pi_integ[0]*(1./(2*pi**2))
         
    # Populate matrix elements.
    X = unique_phi_correlations[inv_idx]
    P = unique_pi_correlations[inv_idx]
     
    # Diagonalize.
    XP = sp.dot(X,P)
    v = linalg.eigvalsh(XP)
    sqrt_v = sp.sqrt(v)
     
    # Calculate entropy.
    S_n = 0
    if n == 1:
        for vk in sqrt_v:
            S_n += ((vk + 0.5)*log(vk + 0.5) - (vk - 0.5)*log(vk - 0.5))
    else:
        for vk in sqrt_v:
            S_n += log((vk + 0.5)**n - (vk - 0.5)**n)
        S_n *= 1./(n-1)
         
    return S_n