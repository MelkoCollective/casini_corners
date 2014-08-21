#!/usr/bin/python
'''
Created on May 27, 2013

This script can calculate the coefficients in scaling laws of geometries on a 
free boson square lattice. In its present state, it can calculate the corner
coefficient in the square and the gamma in the circle but it is simple to add
more geometries and scaling laws.

It uses the method outlined in the appendix of Casini and Huerta's paper:
    http://arxiv.org/abs/hep-th/0606256
We also acknowledge that we have corresponded with Casini to get some tips
that have been helpful in overcoming some obstacles.

@author: ghwatson
'''
# TODO: Replace maplelink with sympy once they have elliptic integrals fully worked in.
#      Pull request here: 
#            https://github.com/sympy/sympy/pull/2165
#      Discussion on adding elliptic integrals here: 
#            http://colabti.org/irclogger/irclogger_log/sympy?date=2013-06-09#l151
#      My original query here:
#            http://stackoverflow.com/questions/16991997/decomposing-integral-into-elliptic-integrals/17013218?noredirect=1#comment24619933_17013218

import sys

from scipy import optimize, zeros, log, arange

import calculate as calc
import pickle

def main():
    geometry_set = sys.argv[1] # choices = ['SQUARE','CIRCLE']
    precision = int(sys.argv[2])    # precision = 15 or higher recommended
    n = int(sys.argv[3])            # the Renyi entropy index
    maple_dir = sys.argv[4]         # directory of maple executable
    fit_to = int(sys.argv[5])       # bound on the geometry characteristic length
    load_pickle = bool(int(sys.argv[6])) # loads old correlators (to avoid recomputing)
    save_pickle = bool(int(sys.argv[7])) # saves new correlators
    
    # The geometry is ultimately defined by a function which can generate the
    # lattice points, and the function to which to fit.
    # p0 is just a starting point for the curve_fit optimization routine.
    # This if-elif clause sets up these quantities.
    if geometry_set == 'SQUARE':
        calc_pointcloud = calc.square_lattice
        def func_to_fit(L,c0,c1,c2,c3,s_n):
            # Note that the coefficient's names are not the same as Casini's notation.
            return c0 + c1*(L) + c2*(1./(L)) + c3*(1./(L)**2) - 4*s_n*log(L)
        p0 = [1,1,1,1,0.05]
    elif geometry_set == 'CIRCLE':
        calc_pointcloud = calc.circle_lattice
        def func_to_fit(L,s_n,gamma):
            return gamma + s_n*L
        p0 = [1,0.05]
    
    # Storage to avoid repeated computation of ij correlations. If desired,
    # can load in pickled correlators from a previous run.
    saved_correlations = {}
    if load_pickle is True:
        pkl_file = open('corrs.pkl', 'rb')
        saved_correlations = pickle.load(pkl_file)
        pkl_file.close()
        
    # Print presets
    print "START-------------------------------------"
    print "------------------------------------------"
    print "PARAMETERS:"
    print 'geometry: {0}'.format(geometry_set)
    print "precision: {0}".format(precision)
    print "Renyi index: {0}".format(n)
    print "fit_to: {0}".format(fit_to)
    print "------------------------------------------"

    # Get the entropy for many different sizes of square lattice.
    sizes = arange(1,fit_to)

    entropies = zeros(sizes.shape)
    for count,L in enumerate(sizes):
        print "---------------Working on lattice size L={0}...".format(L)
        polygon = calc_pointcloud(L)
#         X,P,saved_correlations = calc.correlations(polygon,maple_dir,precision,saved_correlations,True)
        X,P,saved_correlations = calc.correlations_multicore(polygon,maple_dir,precision,saved_correlations,True)
        entropies[count] = calc.entropy(X,P,n,precision, True, True)
        
        # Save the 'saved' correlations to file for later use.
        if save_pickle is True:
            output = open('corrs.pkl','wb')
            pickle.dump(saved_correlations,output,-1)
            output.close()
        
    # Fitting.
    print "Performing fit..."
    popt,_ = optimize.curve_fit(func_to_fit,sizes,entropies,p0)

    # Get the corner coefficient (square) or gamma (circle)
    s_n = popt[-1]
    print "Coefficient for {0} is: ".format(geometry_set) 
    print s_n


if __name__ == '__main__':    
    main()
